/*
Encoder for COFDMTV

Copyright 2022 Ahmet Inan <inan@aicodix.de>
*/

#pragma once

#include <cmath>
#include <iostream>
#include <algorithm>
#include "bose_chaudhuri_hocquenghem_encoder.hh"
#include "base37_bitmap.hh"
#include "xorshift.hh"
#include "complex.hh"
#include "bitman.hh"
#include "polar.hh"
#include "utils.hh"
#include "const.hh"
#include "papr.hh"
#include "fft.hh"
#include "mls.hh"
#include "crc.hh"
#include "psk.hh"
#include "qam.hh"
#include "wav.hh"

int log2_int(int n) {
    if (n == 0) {
        std::cout << "Number is 0" << std::endl;
        return 0;
    }
    int log = 0;
    while (n >>= 1) ++log;
    return log;
}

struct EncoderInterface {
	virtual void configure(const uint8_t *, const int8_t *, int, int, bool) = 0;

	virtual bool produce(int16_t *, int) = 0;

	virtual bool produce_write(int) = 0;

	virtual int rate() = 0;

	virtual ~EncoderInterface() = default;
};

template<int RATE, int SYMBOL_MAPPING>
class Encoder : public EncoderInterface {
	typedef DSP::Complex<float> cmplx;
	typedef DSP::Const<float> Const;
	typedef int8_t code_type;
	static const int code_order = 11;
	static const int map = SYMBOL_MAPPING;
	int mod_bits = log2_int(SYMBOL_MAPPING);
	static const int symbol_count = (SYMBOL_MAPPING == 2)? 8 : 4;
	static const int code_len = 1 << code_order;
	static const int symbol_length = (1280 * RATE) / 8000;
	static const int guard_length = symbol_length / 8;
	static const int extended_length = symbol_length + guard_length;
	static const int max_bits = 1360;
	static const int cor_seq_len = 127;
	static const int cor_seq_off = 1 - cor_seq_len;
	static const int cor_seq_poly = 0b10001001;
	static const int pre_seq_len = 255;
	static const int pre_seq_off = -pre_seq_len / 2;
	static const int pre_seq_poly = 0b100101011;
	static const int pay_car_cnt = 256;
	static const int pay_car_off = -pay_car_cnt / 2;
	static const int fancy_off = -(8 * 9 * 3) / 2;
	static const int noise_poly = 0b100101010001;
	DSP::WritePCM<float> *pcm;
	DSP::FastFourierTransform<symbol_length, cmplx, 1> bwd;
	CODE::CRC<uint16_t> crc;
	CODE::BoseChaudhuriHocquenghemEncoder<255, 71> bch;
	CODE::MLS noise_seq;
	ImprovePAPR<cmplx, symbol_length, RATE <= 16000 ? 4 : 1> improve_papr;
	PolarEncoder<code_type> polar;
	cmplx temp[extended_length], freq[symbol_length], prev[pay_car_cnt], guard[guard_length];
	uint8_t mesg[max_bits / 8], call[9];
	code_type code[code_len];
	uint64_t meta_data;
	int operation_mode = 0;
	int carrier_offset = 0;
	int symbol_number = symbol_count;
	int count_down = 0;
	int fancy_line = 0;
	int noise_count = 0;

	static uint8_t base37_map(int8_t c) {
		if (c >= '0' && c <= '9')
			return c - '0' + 1;
		if (c >= 'a' && c <= 'z')
			return c - 'a' + 11;
		if (c >= 'A' && c <= 'Z')
			return c - 'A' + 11;
		return 0;
	}

	static uint64_t base37(const int8_t *str) {
		uint64_t acc = 0;
		for (char c = *str++; c; c = *str++)
			acc = 37 * acc + base37_map(c);
		return acc;
	}

	static int nrz(bool bit) {
		return 1 - 2 * bit;
	}

	static cmplx mod_map(code_type *b) {
		switch (SYMBOL_MAPPING) {
			case 2:
				return PhaseShiftKeying<2, cmplx, code_type>::map(b);
				break;
			case 4:
				return PhaseShiftKeying<4, cmplx, code_type>::map(b);
				break;
			case 8:
				return PhaseShiftKeying<8, cmplx, code_type>::map(b);
				break;
			case 16:
				return QAM<16, cmplx, code_type>::map(b);
				break;
			default:
				return PhaseShiftKeying<4, cmplx, code_type>::map(b);
				break;
		}
	}

	int bin(int carrier) {
		return (carrier + carrier_offset + symbol_length) % symbol_length;
	}

	void schmidl_cox() {
		CODE::MLS seq(cor_seq_poly);
		float factor = std::sqrt(float(2 * symbol_length) / cor_seq_len);
		for (int i = 0; i < symbol_length; ++i)
			freq[i] = 0;
		freq[bin(cor_seq_off - 2)] = factor;
		for (int i = 0; i < cor_seq_len; ++i)
			freq[bin(2 * i + cor_seq_off)] = nrz(seq());
		for (int i = 0; i < cor_seq_len; ++i)
			freq[bin(2 * i + cor_seq_off)] *= freq[bin(2 * (i - 1) + cor_seq_off)];
		transform(false);
	}

	void preamble() {
		uint8_t data[9] = {0}, parity[23] = {0};
		for (int i = 0; i < 55; ++i)
			CODE::set_be_bit(data, i, (meta_data >> i) & 1);
		crc.reset();
		uint16_t cs = crc(meta_data << 9);
		for (int i = 0; i < 16; ++i)
			CODE::set_be_bit(data, i + 55, (cs >> i) & 1);
		bch(data, parity);
		CODE::MLS seq(pre_seq_poly);
		float factor = std::sqrt(float(symbol_length) / pre_seq_len);
		for (int i = 0; i < symbol_length; ++i)
			freq[i] = 0;
		freq[bin(pre_seq_off - 1)] = factor;
		for (int i = 0; i < 71; ++i)
			freq[bin(i + pre_seq_off)] = nrz(CODE::get_be_bit(data, i));
		for (int i = 71; i < pre_seq_len; ++i)
			freq[bin(i + pre_seq_off)] = nrz(CODE::get_be_bit(parity, i - 71));
		for (int i = 0; i < pre_seq_len; ++i)
			freq[bin(i + pre_seq_off)] *= freq[bin(i - 1 + pre_seq_off)];
		for (int i = 0; i < pre_seq_len; ++i)
			freq[bin(i + pre_seq_off)] *= nrz(seq());
		for (int i = 0; i < pay_car_cnt; ++i)
			prev[i] = freq[bin(i + pay_car_off)];
		transform();
	}

	void fancy_symbol() {
		int active_carriers = 1;
		for (int j = 0; j < 9; ++j)
			for (int i = 0; i < 8; ++i)
				active_carriers += (base37_bitmap[call[j] + 37 * fancy_line] >> i) & 1;
		float factor = std::sqrt(float(symbol_length) / active_carriers);
		for (int i = 0; i < symbol_length; ++i)
			freq[i] = 0;
		for (int j = 0; j < 9; ++j)
			for (int i = 0; i < 8; ++i)
				if (base37_bitmap[call[j] + 37 * fancy_line] & (1 << (7 - i)))
					freq[bin((8 * j + i) * 3 + fancy_off)] = factor * nrz(noise_seq());
		transform(false);
	}

	void noise_symbol() {
		float factor = std::sqrt(symbol_length / float(pay_car_cnt));
		for (int i = 0; i < symbol_length; ++i)
			freq[i] = 0;
		for (int i = 0; i < pay_car_cnt; ++i)
			freq[bin(i + pay_car_off)] = factor * cmplx(nrz(noise_seq()), nrz(noise_seq()));
		transform(false);
	}

	void payload_symbol() {
		for (int i = 0; i < symbol_length; ++i)
			freq[i] = 0;
		for (int i = 0; i < pay_car_cnt; ++i)
			freq[bin(i + pay_car_off)] = prev[i] *= mod_map(code + mod_bits * (pay_car_cnt * symbol_number + i));
		transform();
	}

	void silence() {
		for (int i = 0; i < symbol_length; ++i)
			temp[i] = 0;
	}

	void transform(bool papr_reduction = true) {
		if (papr_reduction && RATE <= 16000)
			improve_papr(freq);
		bwd(temp, freq);
		for (int i = 0; i < symbol_length; ++i)
			temp[i] /= std::sqrt(float(8 * symbol_length));
	}

	void next_sample(int16_t *samples, cmplx signal, int channel, int i) {
		switch (channel) {
			case 1:
				samples[2 * i] = std::clamp<float>(std::nearbyint(32767 * signal.real()), -32768, 32767);
				samples[2 * i + 1] = 0;
				break;
			case 2:
				samples[2 * i] = 0;
				samples[2 * i + 1] = std::clamp<float>(std::nearbyint(32767 * signal.real()), -32768, 32767);
				break;
			case 4:
				samples[2 * i] = std::clamp<float>(std::nearbyint(32767 * signal.real()), -32768, 32767);
				samples[2 * i + 1] = std::clamp<float>(std::nearbyint(32767 * signal.imag()), -32768, 32767);
				break;
			default:
				samples[i] = std::clamp<float>(std::nearbyint(32767 * signal.real()), -32768, 32767);
		}
	}

public:
	Encoder(DSP::WritePCM<float> *pcm) : pcm(pcm), noise_seq(noise_poly), crc(0xA8F4), bch({
		0b100011101, 0b101110111, 0b111110011, 0b101101001,
		0b110111101, 0b111100111, 0b100101011, 0b111010111,
		0b000010011, 0b101100101, 0b110001011, 0b101100011,
		0b100011011, 0b100111111, 0b110001101, 0b100101101,
		0b101011111, 0b111111001, 0b111000011, 0b100111001,
		0b110101001, 0b000011111, 0b110000111, 0b110110001}) {}

	int rate() final {
		return RATE;
	}



	bool produce(int16_t *audio_buffer, int channel_select) final {
		bool data_symbol = false;
		switch (count_down) {
			case 5:
				if (noise_count) {
					--noise_count;
					noise_symbol();
					break;
				}
				--count_down;
			case 4:
				schmidl_cox();
				data_symbol = true;
				--count_down;
				break;
			case 3:
				preamble();
				data_symbol = true;
				--count_down;
				if (!operation_mode)
					--count_down;
				break;
			case 2:
				payload_symbol();
				data_symbol = true;
				if (++symbol_number == symbol_count)
					--count_down;
				break;
			case 1:
				if (fancy_line) {
					--fancy_line;
					fancy_symbol();
					break;
				}
				silence();
				--count_down;
				break;
			default:
				for (int i = 0; i < extended_length; ++i)
					next_sample(audio_buffer, 0, channel_select, i);
				return false;
		}
		for (int i = 0; i < guard_length; ++i) {
			float x = i / float(guard_length - 1);
			float ratio(0.5);
			if (data_symbol)
				x = std::min(x, ratio) / ratio;
			float y = 0.5f * (1 - std::cos(DSP::Const<float>::Pi() * x));
			cmplx sum = DSP::lerp(guard[i], temp[i + symbol_length - guard_length], y);
			next_sample(audio_buffer, sum, channel_select, i);
		}
		for (int i = 0; i < guard_length; ++i)
			guard[i] = temp[i];
		for (int i = 0; i < symbol_length; ++i)
			next_sample(audio_buffer, temp[i], channel_select, i + guard_length);
		return true;
	}

	bool produce_write(int channel) {
		bool data_symbol = false;
		switch (count_down) {
			case 5:
				if (noise_count) {
					--noise_count;
					noise_symbol();
					break;
				}
				--count_down;
			case 4:
				schmidl_cox();
				data_symbol = true;
				--count_down;
				break;
			case 3:
				preamble();
				data_symbol = true;
				--count_down;
				if (!operation_mode)
					--count_down;
				break;
			case 2:
				payload_symbol();
				data_symbol = true;
				if (++symbol_number == symbol_count)
					--count_down;
				break;
			case 1:
				if (fancy_line) {
					--fancy_line;
					fancy_symbol();
					break;
				}
				silence();
				--count_down;
				break;
			default:
				return false;
		}
		for (int i = 0; i < guard_length; ++i) {
			float x = i / float(guard_length - 1);
			float ratio(0.5);
			if (data_symbol)
				x = std::min(x, ratio) / ratio;
			float y = 0.5f * (1 - std::cos(DSP::Const<float>::Pi() * x));
			cmplx sum = DSP::lerp(guard[i], temp[i + symbol_length - guard_length], y);
			float sum_write[2];
			switch (channel) {
				case 1:
					sum_write[0] = sum.real(); 
					sum_write[1] = 0.0;
					pcm->write(sum_write, 1, 2);
					break;
				case 2:
					sum_write[0] = 0.0; 
					sum_write[1] = sum.real();
					pcm->write(sum_write, 1, 2);
					break;
				case 3:
					sum_write[0] = sum.real(); 
					sum_write[1] = sum.imag();
					pcm->write(sum_write, 1, 2);
					break;
				default:
					sum_write[0] = sum.real();
					pcm->write(sum_write, 1, 1);
			}
		}
		for (int i = 0; i < guard_length; ++i)
			guard[i] = temp[i];
		for (int i = 0; i < symbol_length; ++i) {
			float temp_write[2];
			switch (channel) {
				case 1:
					temp_write[0] = temp[i].real(); 
					temp_write[1] = 0.0;
					pcm->write(temp_write, 1, 2);
					break;
				case 2:
					temp_write[0] = 0.0; 
					temp_write[1] = temp[i].real();
					pcm->write(temp_write, 1, 2);
					break;
				case 3:
					temp_write[0] = temp[i].real(); 
					temp_write[1] = temp[i].imag();
					pcm->write(temp_write, 1, 2);
					break;
				default:
					temp_write[0] = temp[i].real();
					pcm->write(temp_write, 1, 1);
			}
		}
			
		return true;
	}

	void configure(const uint8_t *payload, const int8_t *call_sign, int carrier_frequency, int noise_symbols, bool fancy_header) final {
		int len = 0;
		while (len <= 128 && payload[len])
			++len;
		if (!len)
			operation_mode = 0;
		else if (len <= 85)
			operation_mode = 16;
		else if (len <= 128)
			operation_mode = 15;
		else
			operation_mode = 14;
		carrier_offset = (carrier_frequency * symbol_length) / RATE;
		meta_data = (base37(call_sign) << 8) | operation_mode;
		for (int i = 0; i < 9; ++i)
			call[i] = 0;
		for (int i = 0; i < 9 && call_sign[i]; ++i)
			call[i] = base37_map(call_sign[i]);
		symbol_number = 0;
		count_down = 5;
		fancy_line = 11 * fancy_header;
		noise_count = noise_symbols;
		for (int i = 0; i < guard_length; ++i)
			guard[i] = 0;
		const uint32_t *frozen_bits;
		int data_bits;
		switch (operation_mode) {
			case 14:
				data_bits = 1360;
				frozen_bits = frozen_2048_1392;
				break;
			case 15:
				data_bits = 1024;
				frozen_bits = frozen_2048_1056;
				break;
			case 16:
				data_bits = 680;
				frozen_bits = frozen_2048_712;
				break;
			default:
				return;
		}
		CODE::Xorshift32 scrambler;
		for (int i = 0; i < data_bits / 8; ++i)
			mesg[i] = payload[i] ^ scrambler();
		polar(code, mesg, frozen_bits, data_bits);
	}
};
