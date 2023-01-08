/*
Some finite impulse response filter functions

Copyright 2018 Ahmet Inan <inan@aicodix.de>
*/

#pragma once

#include "const.hh"
#include "utils.hh"
#include "unit_circle.hh"

namespace DSP {

template <typename TYPE>
class LowPass
{
	TYPE f;
public:
	LowPass(TYPE cutoff) : f(TYPE(2) * cutoff) {}
	TYPE operator () (int n, int N) const
	{
		TYPE x = TYPE(n) - TYPE(0.5) * TYPE(N - 1);
		return f * sinc(f * x);
	}
};

template <typename TYPE>
class LowPass2
{
	int num, den;
	TYPE fac;
public:
	LowPass2(int num, int den) : num(num), den(den), fac(TYPE(2*num)/TYPE(den)) {}
	TYPE operator () (int n, int N) const
	{
		int twox = 2 * n - (N - 1);
		return !twox ? fac : fac *
			UnitCircle<TYPE>::sin((twox * num) % (2 * den), 2 * den) /
			(Const<TYPE>::HalfPi() * fac * TYPE(twox));
	}
};

template <typename TYPE>
class HighPass
{
	TYPE f;
public:
	HighPass(TYPE cutoff) : f(TYPE(2) * cutoff) {}
	TYPE operator () (int n, int N) const
	{
		TYPE x = TYPE(n) - TYPE(0.5) * TYPE(N - 1);
		// if (N%1) return delta(x) - f * sinc(f * x);
		return sinc(x) - f * sinc(f * x);

	}
};

template <typename TYPE>
class HighPass2
{
	int num, den;
	TYPE fac;
public:
	HighPass2(int num, int den) : num(num), den(den), fac(TYPE(2*num)/TYPE(den)) {}
	TYPE operator () (int n, int N) const
	{
		int twox = 2 * n - (N - 1);
		return !twox ? TYPE(1) - fac :
			UnitCircle<TYPE>::sin(twox % 4, 4) / (Const<TYPE>::HalfPi() * TYPE(twox))
			- fac * UnitCircle<TYPE>::sin((twox * num) % (2 * den), 2 * den) /
			(Const<TYPE>::HalfPi() * fac * TYPE(twox));
	}
};

template <typename TYPE>
class BandPass
{
	TYPE f0, f1;
public:
	BandPass(TYPE cutoff0, TYPE cutoff1) :
		f0(TYPE(2) * cutoff0), f1(TYPE(2) * cutoff1) {}
	TYPE operator () (int n, int N) const
	{
		TYPE x = TYPE(n) - TYPE(0.5) * TYPE(N - 1);
		return f1 * sinc(f1 * x) - f0 * sinc(f0 * x);
	}
};

template <typename TYPE>
struct HilbertTransform
{
	TYPE operator () (int n, int N) const
	{
		if (N&1) {
			int x = n - (N - 1) / 2;
			return x&1 ? TYPE(2) / (Const<TYPE>::Pi() * TYPE(x)) : TYPE(0);
		} else {
			TYPE x = TYPE(n) - TYPE(0.5) * TYPE(N - 1);
			return TYPE(1) / (Const<TYPE>::Pi() * x);
		}
	}
};

}

