/*
Discrete Hilbert transformation

Copyright 2020 Ahmet Inan <inan@aicodix.de>
*/

#pragma once

#include "window.hh"

namespace DSP {

template <typename TYPE, int TAPS>
class Hilbert
{
	static_assert((TAPS-1) % 4 == 0, "TAPS-1 not divisible by four");
	typedef TYPE complex_type;
	typedef typename TYPE::value_type value_type;
	value_type real[TAPS];
	value_type imco[(TAPS-1)/4];
	value_type reco;
public:
	Hilbert(value_type a = value_type(2))
	{
		Kaiser<value_type> win(a);
		reco = win((TAPS-1)/2, TAPS);
		for (int i = 0; i < (TAPS-1)/4; ++i)
			imco[i] = win((2*i+1)+(TAPS-1)/2, TAPS) * 2 / ((2*i+1) * Const<value_type>::Pi());
		for (int i = 0; i < TAPS; ++i)
			real[i] = 0;
	}
	complex_type operator()(value_type input)
	{
		value_type re = reco * real[(TAPS-1)/2];
		value_type im = imco[0] * (real[(TAPS-1)/2-1] - real[(TAPS-1)/2+1]);
		for (int i = 1; i < (TAPS-1)/4; ++i)
			im += imco[i] * (real[(TAPS-1)/2-(2*i+1)] - real[(TAPS-1)/2+(2*i+1)]);
		for (int i = 0; i < TAPS-1; ++i)
			real[i] = real[i+1];
		real[TAPS-1] = input;
		return complex_type(re, im);
	}
};

}

