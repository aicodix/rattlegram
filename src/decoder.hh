/*
Decoder for COFDMTV

Copyright 2022 Ahmet Inan <inan@aicodix.de>
*/

#pragma once

#include <cmath>
#include <iostream>
#include <algorithm>
#include <cassert>

namespace DSP { using std::abs; using std::min; using std::cos; using std::sin; }

#include "schmidl_cox.hh"
#include "bip_buffer.hh"
#include "theil_sen.hh"
#include "xorshift.hh"
#include "decibel.hh"
#include "complex.hh"
#include "hilbert.hh"
#include "blockdc.hh"
#include "filter.hh"
#include "window.hh"
#include "coeffs.hh"
#include "bitman.hh"
#include "phasor.hh"
#include "const.hh"
#include "image.hh"
#include "polar.hh"
#include "fft.hh"
#include "mls.hh"
#include "crc.hh"
#include "osd.hh"
#include "psk.hh"
#include "qam.hh"

#define STATUS_OKAY 0
#define STATUS_FAIL 1
#define STATUS_SYNC 2
#define STATUS_DONE 3
#define STATUS_HEAP 4
#define STATUS_NOPE 5
#define STATUS_PING 6

int log2_int(int n) {
    if (n == 0) {
        std::cout << "Number is 0" << std::endl;
        return 0;
    }
    int log = 0;
    while (n >>= 1) ++log;
    return log;
}

struct DecoderInterface {
	virtual bool feed(const int16_t *, int, int) = 0;

	virtual int process() = 0;

	virtual void spectrum(uint32_t *, uint32_t *, int) = 0;

	virtual void staged(float *, int32_t *, uint8_t *) = 0;

	virtual int fetch(uint8_t *) = 0;

	virtual int rate() = 0;

	virtual ~DecoderInterface() = default;
};

template<int RATE, int SYMBOL_MAPPING>
class Decoder : public DecoderInterface {
	typedef DSP::Complex<float> cmplx;
	typedef DSP::Const<float> Const;
	typedef int8_t code_type;
	static const int spectrum_width = 360, spectrum_height = 128;
	static const int spectrogram_width = 360, spectrogram_height = 128;
	static const int code_order = 11;
	static const int map = SYMBOL_MAPPING;
	int mod_bits = log2_int(SYMBOL_MAPPING);
	static const int symbol_count = (SYMBOL_MAPPING == 2)? 8 : 4;
	static const int code_len = 1 << code_order;
	static const int symbol_length = (1280 * RATE) / 8000;
	static const int guard_length = symbol_length / 8;
	static const int extended_length = symbol_length + guard_length;
	static const int filter_length = (((33 * RATE) / 8000) & ~3) | 1;
	static const int stft_length = extended_length / 2;
	static const int window_length = 2 * stft_length;
	static const int dB_min = -96, dB_max = 0;
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
	DSP::FastFourierTransform<stft_length, cmplx, -1> stft;
	SchmidlCox<float, cmplx, search_position, symbol_length / 2, guard_length> correlator;
	DSP::BlockDC<float, float> block_dc;
	DSP::Hilbert<cmplx, filter_length> hilbert;
	DSP::BipBuffer<cmplx, buffer_length> buffer;
	DSP::TheilSenEstimator<float, pay_car_cnt> tse;
	DSP::Phasor<cmplx> osc;
	DSP::Hann<float> hann;
	DSP::LowPass2<float> lowpass;
	DSP::Coeffs<window_length, float, true> window;
	CODE::CRC<uint16_t> crc;
	CODE::OrderedStatisticsDecoder<255, 71, 2> osd;
	PolarDecoder<code_type> polar;
	cmplx temp[extended_length], freq[symbol_length], prev[pay_car_cnt], cons[pay_car_cnt];
	float power[spectrum_width]{}, index[pay_car_cnt]{}, phase[pay_car_cnt]{};
	code_type code[code_len];
	int8_t generator[255 * 71];
	int8_t soft[pre_seq_len];
	uint8_t data[(pre_seq_len + 7) / 8];
	int symbol_number = symbol_count;
	int symbol_position = search_position + extended_length;
	int stored_position = 0;
	int staged_position = 0;
	int staged_mode = 0;
	int operation_mode = 0;
	int accumulated = 0;
	float stored_cfo_rad = 0;
	float staged_cfo_rad = 0;
	uint64_t staged_call = 0;
	bool stored_check = false;
	bool staged_check = false;
	const cmplx *buf;

	static uint32_t argb(float a, float r, float g, float b) {
		a = std::clamp<float>(a, 0, 1);
		r = std::clamp<float>(r, 0, 1);
		g = std::clamp<float>(g, 0, 1);
		b = std::clamp<float>(b, 0, 1);
		r *= a;
		g *= a;
		b *= a;
		r = std::sqrt(r);
		g = std::sqrt(g);
		b = std::sqrt(b);
		int A = (int) std::nearbyint(255 * a);
		int R = (int) std::nearbyint(255 * r);
		int G = (int) std::nearbyint(255 * g);
		int B = (int) std::nearbyint(255 * b);
		return (A << 24) | (R << 16) | (G << 8) | (B << 0);
	}

	static uint32_t rainbow(float v) {
		v = std::clamp<float>(v, 0, 1);
		float t = 4 * v - 2;
		return argb(4 * v, t, 1 - std::abs(t), -t);
	}

	static int bin(int carrier) {
		return (carrier + symbol_length) % symbol_length;
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

	static void mod_hard(code_type *b, cmplx c) {
		switch (SYMBOL_MAPPING) {
			case 2:
				return PhaseShiftKeying<2, cmplx, code_type>::hard(b, c);
				break;
			case 4:
				return PhaseShiftKeying<4, cmplx, code_type>::hard(b, c);
				break;
			case 8:
				return PhaseShiftKeying<8, cmplx, code_type>::hard(b, c);
				break;
			case 16:
				return QAM<16, cmplx, code_type>::hard(b, c);
				break;
			default:
				return PhaseShiftKeying<4, cmplx, code_type>::hard(b, c);
				break;
		}
	}

	static void mod_soft(code_type *b, cmplx c, float precision) {
		switch (SYMBOL_MAPPING) {
			case 2:
				return PhaseShiftKeying<2, cmplx, code_type>::soft(b, c, precision);
				break;
			case 4:
				return PhaseShiftKeying<4, cmplx, code_type>::soft(b, c, precision);
				break;
			case 8:
				return PhaseShiftKeying<8, cmplx, code_type>::soft(b, c, precision);
				break;
			case 16:
				return QAM<16, cmplx, code_type>::soft(b, c, precision);
				break;
			default:
				return PhaseShiftKeying<4, cmplx, code_type>::soft(b, c, precision);
				break;
		}
	}

	static void base37(uint8_t *str, uint64_t val, int len) {
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

	cmplx convert(const int16_t *samples, int channel, int i) {
		switch (channel) {
			case 1:
				return analytic(samples[2 * i] / 32768.f);
			case 2:
				return analytic(samples[2 * i + 1] / 32768.f);
			case 3:
				return analytic(((int) samples[2 * i] + (int) samples[2 * i + 1]) / 65536.f);
			case 4:
				return cmplx(samples[2 * i], samples[2 * i + 1]) / 32768.f;
		}
		return analytic(samples[i] / 32768.f);
	}

	void update_spectrum(uint32_t *pixels, uint32_t tint) {
		Image<uint32_t, spectrum_width, spectrum_height> img(pixels);
		img.fill(0);
		auto pos = [this, img](int i) {
			return (int) std::nearbyint((1 - power[i]) * (img.height - 1));
		};
		tint |= 0xff000000;
		for (int i = 1, j = pos(0), k; i < img.width; ++i, j = k)
			img.line(i - 1, j, i, k = pos(i), tint);
	}

	void update_spectrogram(uint32_t *pixels) {
		std::memmove(pixels + spectrogram_width, pixels, sizeof(uint32_t) * spectrogram_width * (spectrogram_height - 1));
		for (int i = 0; i < spectrogram_width; ++i)
			pixels[i] = rainbow(power[i]);
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

	int preamble() {
		DSP::Phasor<cmplx> nco;
		nco.omega(-staged_cfo_rad);
		for (int i = 0; i < symbol_length; ++i)
			temp[i] = buf[staged_position + i] * nco();
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
		staged_mode = md & 255;
		staged_call = md >> 8;
		if (staged_mode && (staged_mode < 14 || staged_mode > 16))
			return STATUS_NOPE;
		if (staged_call == 0 || staged_call >= 129961739795077L) {
			staged_call = 0;
			return STATUS_NOPE;
		}
		if (!staged_mode)
			return STATUS_PING;
		return STATUS_OKAY;
	}

public:
	Decoder() : correlator(corSeq()), crc(0xA8F4), lowpass(1, symbol_length), window(&hann, &lowpass) {
		CODE::BoseChaudhuriHocquenghemGenerator<255, 71>::matrix(generator, true, {
			0b100011101, 0b101110111, 0b111110011, 0b101101001,
			0b110111101, 0b111100111, 0b100101011, 0b111010111,
			0b000010011, 0b101100101, 0b110001011, 0b101100011,
			0b100011011, 0b100111111, 0b110001101, 0b100101101,
			0b101011111, 0b111111001, 0b111000011, 0b100111001,
			0b110101001, 0b000011111, 0b110000111, 0b110110001});
		block_dc.samples(filter_length);
		osc.omega(-2000, RATE);
	}

	int rate() final {
		return RATE;
	}

	void staged(float *cfo, int32_t *mode, uint8_t *call) final {
		*cfo = staged_cfo_rad * (RATE / Const::TwoPi());
		*mode = staged_mode;
		base37(call, staged_call, 9);
	}

	int fetch(uint8_t *payload) final {
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
				return -1;
		}
		int result = polar(payload, code, frozen_bits, data_bits);
		CODE::Xorshift32 scrambler;
		for (int i = 0; i < data_bits / 8; ++i)
			payload[i] ^= scrambler();
		for (int i = data_bits / 8; i < 170; ++i)
			payload[i] = 0;
		return result;
	}

	bool feed(const int16_t *audio_buffer, int sample_count, int channel_select) final {
		assert(sample_count <= extended_length);
		for (int i = 0; i < sample_count; ++i) {
			if (correlator(buffer(convert(audio_buffer, channel_select, i)))) {
				stored_cfo_rad = correlator.cfo_rad;
				stored_position = correlator.symbol_pos + accumulated;
				stored_check = true;
			}
			if (++accumulated == extended_length)
				buf = buffer();
		}
		if (accumulated >= extended_length) {
			accumulated -= extended_length;
			if (stored_check) {
				staged_cfo_rad = stored_cfo_rad;
				staged_position = stored_position;
				staged_check = true;
				stored_check = false;
			}
			return true;
		}
		return false;
	}

	int process() final {
		int status = STATUS_OKAY;
		if (staged_check) {
			staged_check = false;
			status = preamble();
			if (status == STATUS_OKAY) {
				operation_mode = staged_mode;
				osc.omega(-staged_cfo_rad);
				symbol_position = staged_position;
				symbol_number = -1;
				status = STATUS_SYNC;
			}
		}
		if (symbol_number < symbol_count) {
			for (int i = 0; i < extended_length; ++i)
				temp[i] = buf[symbol_position + i] * osc();
			fwd(freq, temp);
			if (symbol_number >= 0) {
				for (int i = 0; i < pay_car_cnt; ++i)
					cons[i] = demod_or_erase(freq[bin(i + pay_car_off)], prev[i]);
				compensate();
				demap();
			}
			if (++symbol_number == symbol_count)
				status = STATUS_DONE;
			for (int i = 0; i < pay_car_cnt; ++i)
				prev[i] = freq[bin(i + pay_car_off)];
		}
		return status;
	}

	void spectrum(uint32_t *spectrum_pixels, uint32_t *spectrogram_pixels, int spectrum_tint) final {
		for (int j = 0; j < 2; ++j) {
			for (int i = 0; i < stft_length; ++i)
				temp[i] = 0;
			for (int i = 0; i < window_length; ++i)
				temp[i % stft_length] += window[i] * buf[buffer_length - window_length + stft_length * (j - 1) + i];
			stft(freq, temp);
			for (int i = 0; i < spectrum_width; ++i)
				power[i] = std::clamp<float>((DSP::decibel(norm(freq[i])) - dB_min) / (dB_max - dB_min), 0, 1);
			update_spectrogram(spectrogram_pixels);
		}
		update_spectrum(spectrum_pixels, spectrum_tint);
	}
};
