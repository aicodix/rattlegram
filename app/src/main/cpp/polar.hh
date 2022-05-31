/*
CA-SCL polar coding for COFDMTV

Copyright 2022 Ahmet Inan <inan@aicodix.de>
*/

#pragma once

#include <cmath>
#include <iostream>
#include <algorithm>

#include "crc.hh"
#include "bitman.hh"
#include "polar_tables.hh"
#include "polar_helper.hh"
#include "polar_encoder.hh"
#include "polar_list_decoder.hh"

template<typename code_type>
class PolarEncoder {
	static const int code_order = 11;
	static const int data_bits = 1360;
	static const int crc_bits = data_bits + 32;
	CODE::CRC<uint32_t> crc;
	CODE::PolarSysEnc<code_type> encode;
	int8_t mesg[crc_bits];

	static int nrz(bool bit) {
		return 1 - 2 * bit;
	}

public:
	PolarEncoder() : crc(0x8F6E37A0) {}

	void operator()(code_type *code, const uint8_t *message) {
		for (int i = 0; i < data_bits; ++i)
			mesg[i] = nrz(CODE::get_le_bit(message, i));
		crc.reset();
		for (int i = 0; i < data_bits / 8; ++i)
			crc(message[i]);
		for (int i = 0; i < 32; ++i)
			mesg[i + data_bits] = nrz((crc() >> i) & 1);
		encode(code, mesg, frozen_2048_1392, code_order);
	}
};

template<typename code_type>
class PolarDecoder {
#ifdef __AVX2__
	typedef SIMD<code_type, 32 / sizeof(code_type)> mesg_type;
#else
	typedef SIMD<code_type, 16 / sizeof(code_type)> mesg_type;
#endif
	typedef typename CODE::PolarHelper<mesg_type>::PATH metric_type;
	static const int code_order = 11;
	static const int code_len = 1 << code_order;
	static const int data_bits = 1360;
	static const int crc_bits = data_bits + 32;
	CODE::CRC<uint32_t> crc;
	CODE::PolarEncoder<mesg_type> encode;
	CODE::PolarListDecoder<mesg_type, code_order> decode;
	mesg_type mesg[crc_bits], mess[code_len];

	void systematic() {
		encode(mess, mesg, frozen_2048_1392, code_order);
		for (int i = 0, j = 0; i < code_len && j < crc_bits; ++i)
			if (!((frozen_2048_1392[i / 32] >> (i % 32)) & 1))
				mesg[j++] = mess[i];
	}

public:
	PolarDecoder() : crc(0x8F6E37A0) {}

	bool operator()(uint8_t *message, const code_type *code) {
		metric_type metric[mesg_type::SIZE];
		decode(metric, mesg, code, frozen_2048_1392, code_order);
		systematic();
		int order[mesg_type::SIZE];
		for (int k = 0; k < mesg_type::SIZE; ++k)
			order[k] = k;
		std::sort(order, order + mesg_type::SIZE, [metric](int a, int b) { return metric[a] < metric[b]; });
		int best = -1;
		for (int k = 0; k < mesg_type::SIZE; ++k) {
			crc.reset();
			for (int i = 0; i < crc_bits; ++i)
				crc(mesg[i].v[order[k]] < 0);
			if (crc() == 0) {
				best = order[k];
				break;
			}
		}
		if (best < 0)
			return false;
		for (int i = 0; i < data_bits; ++i) {
			bool decoded = mesg[i].v[best] < 0;
			CODE::set_le_bit(message, i, decoded);
		}
		return true;
	}
};
