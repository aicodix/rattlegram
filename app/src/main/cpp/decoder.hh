/*
Decoder for COFDMTV

Copyright 2022 Ahmet Inan <inan@aicodix.de>
*/

#pragma once

#include <cmath>
#include <iostream>
#include <algorithm>

namespace DSP { using std::abs; using std::min; using std::cos; using std::sin; }

#include "schmidl_cox.hh"
#include "bip_buffer.hh"
#include "theil_sen.hh"
#include "xorshift.hh"
#include "complex.hh"
#include "hilbert.hh"
#include "blockdc.hh"
#include "bitman.hh"
#include "phasor.hh"
#include "const.hh"
#include "polar.hh"
#include "fft.hh"
#include "mls.hh"
#include "crc.hh"
#include "osd.hh"
#include "psk.hh"

#define STATUS_OKAY 0
#define STATUS_FAIL 1
#define STATUS_SYNC 2
#define STATUS_DONE 3
#define STATUS_HEAP 4
#define STATUS_NOPE 5

struct DecoderInterface {
	virtual int process(const int16_t *, int) = 0;

	virtual void cached(float *, int32_t *, int8_t *) = 0;

	virtual bool fetch(uint8_t *) = 0;

	virtual int rate() = 0;

	virtual ~DecoderInterface() = default;
};

template<int RATE>
class Decoder : public DecoderInterface {
	typedef DSP::Complex<float> cmplx;
	typedef DSP::Const<float> Const;
	typedef int8_t code_type;
	static const int code_order = 11;
	static const int mod_bits = 2;
	static const int code_len = 1 << code_order;
	static const int symbol_count = 4;
	static const int symbol_length = (1280 * RATE) / 8000;
	static const int guard_length = symbol_length / 8;
	static const int extended_length = symbol_length + guard_length;
	static const int filter_length = (((21 * RATE) / 8000) & ~3) | 1;
	static const int data_bits = 1360;
	static const int cor_seq_len = 127;
	static const int cor_seq_off = 1 - cor_seq_len;
	static const int cor_seq_poly = 0b10001001;
	static const int pre_seq_len = 255;
	static const int pre_seq_off = -pre_seq_len / 2;
	static const int pre_seq_poly = 0b100101011;
	static const int pay_car_cnt = 256;
	static const int pay_car_off = -pay_car_cnt / 2;
	static const int buffer_length = 4 * extended_length;
	static const int search_position = extended_length;
	DSP::FastFourierTransform<symbol_length, cmplx, -1> fwd;
	SchmidlCox<float, cmplx, search_position, symbol_length / 2, guard_length> correlator;
	DSP::BlockDC<float, float> block_dc;
	DSP::Hilbert<cmplx, filter_length> hilbert;
	DSP::BipBuffer<cmplx, buffer_length> buffer;
	DSP::TheilSenEstimator<float, pay_car_cnt> tse;
	DSP::Phasor<cmplx> osc;
	CODE::CRC<uint16_t> crc;
	CODE::OrderedStatisticsDecoder<255, 71, 2> osd;
	PolarDecoder<code_type> polar;
	cmplx temp[extended_length], freq[symbol_length], prev[pay_car_cnt], cons[pay_car_cnt];
	float index[pay_car_cnt]{}, phase[pay_car_cnt]{};
	code_type code[code_len];
	int8_t generator[255 * 71];
	int8_t soft[pre_seq_len];
	uint8_t data[(pre_seq_len + 7) / 8];
	int symbol_number = symbol_count;
	int symbol_position = search_position + extended_length;
	int cached_mode = 0;
	uint64_t cached_call = 0;

	static int bin(int carrier) {
		return (carrier + symbol_length) % symbol_length;
	}

	static int nrz(bool bit) {
		return 1 - 2 * bit;
	}

	static cmplx mod_map(code_type *b) {
		return PhaseShiftKeying<4, cmplx, code_type>::map(b);
	}

	static void mod_hard(code_type *b, cmplx c) {
		PhaseShiftKeying<4, cmplx, code_type>::hard(b, c);
	}

	static void mod_soft(code_type *b, cmplx c, float precision) {
		PhaseShiftKeying<4, cmplx, code_type>::soft(b, c, precision);
	}

	static void base37(int8_t *str, uint64_t val, int len) {
		for (int i = len - 1; i >= 0; --i, val /= 37)
			str[i] = " 0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"[val % 37];
	}

	static cmplx demod_or_erase(cmplx curr, cmplx prev) {
		if (norm(prev) <= 0)
			return 0;
		cmplx cons = curr / prev;
		if (norm(cons) > 4)
			return 0;
		return cons;
	}

	const cmplx *corSeq() {
		CODE::MLS seq(cor_seq_poly);
		for (int i = 0; i < symbol_length / 2; ++i)
			freq[i] = 0;
		for (int i = 0; i < cor_seq_len; ++i)
			freq[(i + cor_seq_off / 2 + symbol_length / 2) % (symbol_length / 2)] = nrz(seq());
		return freq;
	}

	cmplx analytic(float real) {
		return hilbert(block_dc(real));
	}

	const cmplx *next_sample(const int16_t *samples, int channel, int i) {
		switch (channel) {
			case 1:
				return buffer(analytic(samples[2 * i] / 32768.f));
			case 2:
				return buffer(analytic(samples[2 * i + 1] / 32768.f));
			case 3:
				return buffer(analytic(((int) samples[2 * i] + (int) samples[2 * i + 1]) / 65536.f));
			case 4:
				return buffer(cmplx(samples[2 * i], samples[2 * i + 1]) / 32768.f);
		}
		return buffer(analytic(samples[i] / 32768.f));
	}

	void compensate() {
		int count = 0;
		for (int i = 0; i < pay_car_cnt; ++i) {
			cmplx con = cons[i];
			if (con.real() != 0 && con.imag() != 0) {
				code_type tmp[mod_bits];
				mod_hard(tmp, con);
				index[count] = i + pay_car_off;
				phase[count] = arg(con * conj(mod_map(tmp)));
				++count;
			}
		}
		tse.compute(index, phase, count);
		for (int i = 0; i < pay_car_cnt; ++i)
			cons[i] *= DSP::polar<float>(1, -tse(i + pay_car_off));
	}

	float precision() {
		float sp = 0, np = 0;
		for (int i = 0; i < pay_car_cnt; ++i) {
			code_type tmp[mod_bits];
			mod_hard(tmp, cons[i]);
			cmplx hard = mod_map(tmp);
			cmplx error = cons[i] - hard;
			sp += norm(hard);
			np += norm(error);
		}
		return sp / np;
	}

	void demap() {
		float pre = precision();
		for (int i = 0; i < pay_car_cnt; ++i)
			mod_soft(code + mod_bits * (symbol_number * pay_car_cnt + i), cons[i], pre);
	}

	int preamble(const cmplx *buf) {
		DSP::Phasor<cmplx> nco;
		nco.omega(-correlator.cfo_rad);
		for (int i = 0; i < symbol_length; ++i)
			temp[i] = buf[correlator.symbol_pos + extended_length + i] * nco();
		fwd(freq, temp);
		CODE::MLS seq(pre_seq_poly);
		for (int i = 0; i < pre_seq_len; ++i)
			freq[bin(i + pre_seq_off)] *= nrz(seq());
		for (int i = 0; i < pre_seq_len; ++i)
			PhaseShiftKeying<2, cmplx, int8_t>::soft(soft + i, demod_or_erase(freq[bin(i + pre_seq_off)], freq[bin(i - 1 + pre_seq_off)]), 32);
		if (!osd(data, soft, generator))
			return STATUS_FAIL;
		uint64_t md = 0;
		for (int i = 0; i < 55; ++i)
			md |= (uint64_t) CODE::get_be_bit(data, i) << i;
		uint16_t cs = 0;
		for (int i = 0; i < 16; ++i)
			cs |= (uint16_t) CODE::get_be_bit(data, i + 55) << i;
		crc.reset();
		if (crc(md << 9) != cs)
			return STATUS_FAIL;
		cached_mode = md & 255;
		cached_call = md >> 8;
		if (cached_mode != 14)
			return STATUS_NOPE;
		if (cached_call == 0 || cached_call >= 129961739795077L) {
			cached_call = 0;
			return STATUS_NOPE;
		}
		return STATUS_OKAY;
	}

public:
	Decoder() : correlator(corSeq()), crc(0xA8F4) {
		CODE::BoseChaudhuriHocquenghemGenerator<255, 71>::matrix(generator, true, {
			0b100011101, 0b101110111, 0b111110011, 0b101101001,
			0b110111101, 0b111100111, 0b100101011, 0b111010111,
			0b000010011, 0b101100101, 0b110001011, 0b101100011,
			0b100011011, 0b100111111, 0b110001101, 0b100101101,
			0b101011111, 0b111111001, 0b111000011, 0b100111001,
			0b110101001, 0b000011111, 0b110000111, 0b110110001});
		block_dc.samples(2 * extended_length);
		osc.omega(-2000, RATE);
	}

	int rate() final {
		return RATE;
	}

	void cached(float *cfo, int32_t *mode, int8_t *call) final {
		*cfo = correlator.cfo_rad * (RATE / Const::TwoPi());
		*mode = cached_mode;
		base37(call, cached_call, 9);
	}

	bool fetch(uint8_t *payload) final {
		bool result = polar(payload, code);
		CODE::Xorshift32 scrambler;
		for (int i = 0; i < data_bits / 8; ++i)
			payload[i] ^= scrambler();
		return result;
	}

	int process(const int16_t *audio_buffer, int channel_select) final {
		int status = STATUS_OKAY;
		const cmplx *buf;
		for (int i = 0; i < extended_length; ++i) {
			buf = next_sample(audio_buffer, channel_select, i);
			if (correlator(buf)) {
				status = preamble(buf);
				if (status == STATUS_OKAY) {
					osc.omega(-correlator.cfo_rad);
					symbol_position = correlator.symbol_pos + i;
					symbol_number = 0;
					status = STATUS_SYNC;
				}
			}
		}
		if (status == STATUS_SYNC || symbol_number < symbol_count) {
			for (int i = 0; i < extended_length; ++i)
				temp[i] = buf[symbol_position + i] * osc();
			fwd(freq, temp);
			if (status != STATUS_SYNC) {
				for (int i = 0; i < pay_car_cnt; ++i)
					cons[i] = demod_or_erase(freq[bin(i + pay_car_off)], prev[i]);
				compensate();
				demap();
				if (++symbol_number == symbol_count)
					status = STATUS_DONE;
			}
			for (int i = 0; i < pay_car_cnt; ++i)
				prev[i] = freq[bin(i + pay_car_off)];
		}
		return status;
	}
};
