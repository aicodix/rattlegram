/*
CA-SCL polar coding for COFDMTV

Copyright 2022 Ahmet Inan <inan@aicodix.de>
*/

#pragma once

#include <cmath>
#include <iostream>

#include "crc.hh"
#include "bitman.hh"
#include "polar_tables.hh"
#include "polar_helper.hh"
#include "polar_encoder.hh"
#include "polar_list_decoder.hh"

template<typename code_type>
class PolarEncoder {
	static const int code_order = 11;
	static const int max_bits = 1360 + 32;
	CODE::CRC<uint32_t> crc;
	CODE::PolarSysEnc<code_type> encode;
	int8_t mesg[max_bits];

	static int nrz(bool bit) {
		return 1 - 2 * bit;
	}

public:
	PolarEncoder() : crc(0x8F6E37A0) {}

	void operator()(code_type *code, const uint8_t *message, const uint32_t *frozen_bits, int data_bits) {
		for (int i = 0; i < data_bits; ++i)
			mesg[i] = nrz(CODE::get_le_bit(message, i));
		crc.reset();
		for (int i = 0; i < data_bits / 8; ++i)
			crc(message[i]);
		for (int i = 0; i < 32; ++i)
			mesg[i + data_bits] = nrz((crc() >> i) & 1);
		encode(code, mesg, frozen_bits, code_order);
	}
};

template<typename code_type>
class PolarDecoder {
#ifdef __AVX2__
	typedef SIMD<code_type, 32 / sizeof(code_type)> mesg_type;
#else
	typedef SIMD<code_type, 16 / sizeof(code_type)> mesg_type;
#endif
	static const int code_order = 11;
	static const int code_len = 1 << code_order;
	static const int max_bits = 1360 + 32;
	CODE::CRC<uint32_t> crc;
	CODE::PolarEncoder<mesg_type> encode;
	CODE::PolarListDecoder<mesg_type, code_order> decode;
	mesg_type mesg[max_bits], mess[code_len];

	void systematic(const uint32_t *frozen_bits, int crc_bits) {
		encode(mess, mesg, frozen_bits, code_order);
		for (int i = 0, j = 0; i < code_len && j < crc_bits; ++i)
			if (!((frozen_bits[i / 32] >> (i % 32)) & 1))
				mesg[j++] = mess[i];
	}

public:
	PolarDecoder() : crc(0x8F6E37A0) {}

	int operator()(uint8_t *message, const code_type *code, const uint32_t *frozen_bits, int data_bits) {
		int crc_bits = data_bits + 32;
		decode(nullptr, mesg, code, frozen_bits, code_order);
		systematic(frozen_bits, crc_bits);
		int best = -1;
		for (int k = 0; k < mesg_type::SIZE; ++k) {
			crc.reset();
			for (int i = 0; i < crc_bits; ++i)
				crc(mesg[i].v[k] < 0);
			if (crc() == 0) {
				best = k;
				break;
			}
		}
		if (best < 0)
			return -1;
		int flips = 0;
		for (int i = 0, j = 0; i < data_bits; ++i, ++j) {
			while ((frozen_bits[j / 32] >> (j % 32)) & 1)
				++j;
			bool received = code[j] < 0;
			bool decoded = mesg[i].v[best] < 0;
			flips += received != decoded;
			CODE::set_le_bit(message, i, decoded);
		}
		return flips;
	}
};
