/*
Coefficient array helper

Copyright 2018 Ahmet Inan <inan@aicodix.de>
*/

#pragma once

namespace DSP {

template <int TAPS, typename TYPE, bool NORMALIZE = false>
class Coeffs
{
	TYPE w[TAPS];
public:
	template <typename FUNC>
	Coeffs(FUNC *func)
	{
		for (int n = 0; n < TAPS; ++n)
			w[n] = (*func)(n, TAPS);
		if (NORMALIZE)
			normalize();
	}
	template <typename FUNC0, typename FUNC1>
	Coeffs(FUNC0 *func0, FUNC1 *func1)
	{
		for (int n = 0; n < TAPS; ++n)
			w[n] = (*func0)(n, TAPS) * (*func1)(n, TAPS);
		if (NORMALIZE)
			normalize();
	}
	void normalize(TYPE divisor = 1)
	{
		TYPE sum(0);
		for (int n = 0; n < TAPS; ++n)
			sum += w[n];
		for (int n = 0; n < TAPS; ++n)
			w[n] /= divisor * abs(sum);
	}
	inline TYPE operator () (int n) { return n >= 0 && n < TAPS ? w[n] : 0; }
	inline operator const TYPE * () const { return w; }
};

}

