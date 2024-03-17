/*
Schmidl & Cox correlator

Copyright 2021 Ahmet Inan <inan@aicodix.de>
*/

#pragma once

#include "fft.hh"
#include "sma.hh"
#include "phasor.hh"
#include "trigger.hh"

template<typename value, typename cmplx, int search_pos, int symbol_len, int guard_len>
class SchmidlCox {
	typedef DSP::Const<value> Const;
	static const int match_len = guard_len | 1;
	static const int match_del = (match_len - 1) / 2;
	DSP::FastFourierTransform<symbol_len, cmplx, -1> fwd;
	DSP::FastFourierTransform<symbol_len, cmplx, 1> bwd;
	DSP::SMA4<cmplx, value, symbol_len, false> cor;
	DSP::SMA4<value, value, 2 * symbol_len, false> pwr;
	DSP::SMA4<value, value, match_len, false> match;
	DSP::Delay<value, match_del> align;
	DSP::SchmittTrigger<value> threshold;
	DSP::FallingEdgeTrigger falling;
	cmplx tmp0[symbol_len], tmp1[symbol_len];
	cmplx kern[symbol_len];
	value timing_max = 0;
	value phase_max = 0;
	int index_max = 0;

	static int bin(int carrier) {
		return (carrier + symbol_len) % symbol_len;
	}

	static cmplx demod_or_erase(cmplx curr, cmplx prev, value pwr) {
		if (norm(curr) > pwr && norm(prev) > pwr) {
			cmplx cons = curr / prev;
			if (norm(cons) < value(4))
				return cons;
		}
		return 0;
	}

public:
	int symbol_pos = 0;
	value cfo_rad = 0;
	value frac_cfo = 0;

	SchmidlCox(const cmplx *sequence) : threshold(value(0.17 * match_len), value(0.19 * match_len)) {
		fwd(kern, sequence);
		for (int i = 0; i < symbol_len; ++i)
			kern[i] = conj(kern[i]) / value(symbol_len);
	}

	bool operator()(const cmplx *samples) {
		cmplx P = cor(samples[search_pos + symbol_len] * conj(samples[search_pos + 2 * symbol_len]));
		value R = value(0.5) * pwr(norm(samples[search_pos + 2 * symbol_len]));
		value min_R = 0.00001 * symbol_len;
		R = std::max(R, min_R);
		value timing = match(norm(P) / (R * R));
		value phase = align(arg(P));

		bool collect = threshold(timing);
		bool process = falling(collect);

		if (!collect && !process)
			return false;

		if (timing_max < timing) {
			timing_max = timing;
			phase_max = phase;
			index_max = match_del;
		} else if (index_max < symbol_len + guard_len + match_del) {
			++index_max;
		}

		if (!process)
			return false;

		frac_cfo = phase_max / value(symbol_len);

		DSP::Phasor<cmplx> osc;
		osc.omega(frac_cfo);
		symbol_pos = search_pos - index_max;
		index_max = 0;
		timing_max = 0;
		for (int i = 0; i < symbol_len; ++i)
			tmp1[i] = samples[i + symbol_pos + symbol_len] * osc();
		fwd(tmp0, tmp1);
		value min_pwr = 0;
		for (int i = 0; i < symbol_len; ++i)
			min_pwr += norm(tmp0[i]);
		min_pwr /= symbol_len;
		for (int i = 0; i < symbol_len; ++i)
			tmp1[i] = demod_or_erase(tmp0[i], tmp0[bin(i - 1)], min_pwr);
		fwd(tmp0, tmp1);
		for (int i = 0; i < symbol_len; ++i)
			tmp0[i] *= kern[i];
		bwd(tmp1, tmp0);

		int shift = 0;
		value peak = 0;
		value next = 0;
		for (int i = 0; i < symbol_len; ++i) {
			value power = norm(tmp1[i]);
			if (power > peak) {
				next = peak;
				peak = power;
				shift = i;
			} else if (power > next) {
				next = power;
			}
		}
		if (peak <= next * 4)
			return false;

		int pos_err = std::nearbyint(arg(tmp1[shift]) * symbol_len / Const::TwoPi());
		if (abs(pos_err) > guard_len / 2)
			return false;
		symbol_pos -= pos_err;

		cfo_rad = shift * (Const::TwoPi() / symbol_len) - frac_cfo;
		if (cfo_rad >= Const::Pi())
			cfo_rad -= Const::TwoPi();
		return true;
	}
};
