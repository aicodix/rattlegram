/*
Improve peak-to-average power ratio

Copyright 2021 Ahmet Inan <inan@aicodix.de>
*/

#pragma once

#include "fft.hh"

template<typename cmplx, int size, int fact>
struct ImprovePAPR {
	typedef typename cmplx::value_type value;
	DSP::FastFourierTransform<fact * size, cmplx, -1> fwd;
	DSP::FastFourierTransform<fact * size, cmplx, 1> bwd;
	cmplx temp[fact * size], over[fact * size];
	bool used[size];

	void operator()(cmplx *freq) {
		for (int i = 0; i < size; ++i)
			used[i] = freq[i].real() || freq[i].imag();
		for (int i = 0; i < size / 2; ++i)
			over[i] = freq[i];
		for (int i = size / 2; i < fact * size - size / 2; ++i)
			over[i] = 0;
		for (int i = size / 2; i < size; ++i)
			over[size * (fact - 1) + i] = freq[i];
		bwd(temp, over);
		value factor = 1 / std::sqrt(value(fact * size));
		for (int i = 0; i < fact * size; ++i)
			temp[i] *= factor;
		for (int i = 0; i < fact * size; ++i) {
			value pwr = norm(temp[i]);
			if (pwr > 1)
				temp[i] /= std::sqrt(pwr);
		}
		fwd(over, temp);
		for (int i = 0; i < size / 2; ++i)
			if (used[i])
				freq[i] = factor * over[i];
		for (int i = size / 2; i < size; ++i)
			if (used[i])
				freq[i] = factor * over[size * (fact - 1) + i];
	}
};

template<typename cmplx, int size>
struct ImprovePAPR<cmplx, size, 1> {
	typedef typename cmplx::value_type value;
	DSP::FastFourierTransform<size, cmplx, -1> fwd;
	DSP::FastFourierTransform<size, cmplx, 1> bwd;
	cmplx temp[size];
	bool used[size];

	void operator()(cmplx *freq) {
		for (int i = 0; i < size; ++i)
			used[i] = freq[i].real() || freq[i].imag();
		bwd(temp, freq);
		value factor = 1 / std::sqrt(value(size));
		for (int i = 0; i < size; ++i)
			temp[i] *= factor;
		for (int i = 0; i < size; ++i) {
			value pwr = norm(temp[i]);
			if (pwr > 1)
				temp[i] /= std::sqrt(pwr);
		}
		fwd(freq, temp);
		for (int i = 0; i < size; ++i)
			if (used[i])
				freq[i] *= factor;
			else
				freq[i] = 0;
	}
};
