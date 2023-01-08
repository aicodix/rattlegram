/*
Some window functions

Copyright 2018 Ahmet Inan <inan@aicodix.de>
*/

#pragma once

#include "const.hh"
#include "kahan.hh"
#include "utils.hh"
#include "unit_circle.hh"

namespace DSP {

template <typename TYPE>
struct Rect
{
	TYPE operator () (int n, int N) const { return n >= 0 && n < N ? 1 : 0; }
};

template <typename TYPE>
struct Hann
{
	TYPE operator () (int n, int N) const
	{
		return TYPE(0.5) * (TYPE(1) - UnitCircle<TYPE>::cos(n, N - 1));
	}
};

template <typename TYPE>
struct Hamming
{
	TYPE operator () (int n, int N) const
	{
		return TYPE(0.54) - TYPE(0.46) * UnitCircle<TYPE>::cos(n, N - 1);
	}
};

template <typename TYPE>
struct Lanczos
{
	TYPE operator () (int n, int N) const
	{
#if 0
		return sinc(TYPE(2 * n) / TYPE(N - 1) - TYPE(1));
#else
		return 2*n == N-1 ? TYPE(1) :
			UnitCircle<TYPE>::sin(2*n-(N-1), 2*(N-1)) / (Const<TYPE>::Pi()*TYPE(2*n-(N-1))/TYPE(N-1));
#endif
	}
};

template <typename TYPE>
class Blackman
{
	TYPE a0, a1, a2;
public:
	Blackman(TYPE a0, TYPE a1, TYPE a2) : a0(a0), a1(a1), a2(a2) {}
	Blackman(TYPE a) : Blackman((TYPE(1) - a) / TYPE(2), TYPE(0.5), a / TYPE(2)) {}
	// "exact Blackman"
	Blackman() : Blackman(TYPE(7938) / TYPE(18608), TYPE(9240) / TYPE(18608), TYPE(1430) / TYPE(18608)) {}
	TYPE operator () (int n, int N) const
	{
		return a0 - a1 * UnitCircle<TYPE>::cos(n, N-1) + a2 * UnitCircle<TYPE>::cos((2*n)%(N-1), N-1);
	}
};

template <typename TYPE>
class Gauss
{
	TYPE o;
public:
	Gauss(TYPE o) : o(o) {}
	TYPE operator () (int n, int N) const
	{
		return exp(- TYPE(0.5) * pow((TYPE(n) - TYPE(N - 1) / TYPE(2)) / (o * TYPE(N - 1) / TYPE(2)), TYPE(2)));
	}
};

template <typename TYPE>
class Kaiser
{
	TYPE a;
	/*
	i0() implements the zero-th order modified Bessel function of the first kind:
	https://en.wikipedia.org/wiki/Bessel_function#Modified_Bessel_functions:_I%CE%B1,_K%CE%B1
	$I_\alpha(x) = i^{-\alpha} J_\alpha(ix) = \sum_{m=0}^\infty \frac{1}{m!\, \Gamma(m+\alpha+1)}\left(\frac{x}{2}\right)^{2m+\alpha}$
	$I_0(x) = J_0(ix) = \sum_{m=0}^\infty \frac{1}{m!\, \Gamma(m+1)}\left(\frac{x}{2}\right)^{2m} = \sum_{m=0}^\infty \left(\frac{x^m}{2^m\,m!}\right)^{2}$
	We obviously can't use the factorial here, so let's get rid of it:
	$= 1 + \left(\frac{x}{2 \cdot 1}\right)^2 + \left(\frac{x}{2 \cdot 1}\cdot \frac{x}{2 \cdot 2}\right)^2 + \left(\frac{x}{2 \cdot 1}\cdot \frac{x}{2 \cdot 2}\cdot \frac{x}{2 \cdot 3}\right)^2 + .. = 1 + \sum_{m=1}^\infty \left(\prod_{n=1}^m \frac{x}{2n}\right)^2$
	*/
	static TYPE i0(TYPE x)
	{
		Kahan<TYPE> sum(1.0);
		TYPE val = 1.0;
		// converges for -3*Pi:3*Pi in less than:
		// float: 25 iterations
		// double: 35 iterations
		for (int n = 1; n < 35; ++n) {
			val *= x / TYPE(2 * n);
			if (sum.same(val * val))
				return sum();
		}
		return sum();
	}
	static TYPE sqr(TYPE x)
	{
		return x * x;
	}
public:
	Kaiser(TYPE a) : a(a) {}
	TYPE operator () (int n, int N) const
	{
		return i0(Const<TYPE>::Pi() * a * sqrt(TYPE(1) - sqr(TYPE(2 * n) / TYPE(N - 1) - TYPE(1)))) / i0(Const<TYPE>::Pi() * a);
	}
};

}

