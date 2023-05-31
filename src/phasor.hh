/*
Numerically controlled oscillator

Copyright 2019 Ahmet Inan <inan@aicodix.de>
*/

#include "unit_circle.hh"

#pragma once

namespace DSP {

template <typename TYPE>
class Phasor
{
	typedef TYPE complex_type;
	typedef typename complex_type::value_type value_type;
	complex_type prev, delta;
public:
	constexpr Phasor() : prev(1, 0), delta(1, 0)
	{
	}
	void omega(int n, int N)
	{
		delta = complex_type(
			UnitCircle<value_type>::cos(n, N),
			UnitCircle<value_type>::sin(n, N));
	}
	void omega(value_type v)
	{
		delta = complex_type(cos(v), sin(v));
	}
	void freq(value_type v)
	{
		omega(Const<value_type>::TwoPi() * v);
	}
	void reset()
	{
		prev = complex_type(1, 0);
	}
	complex_type operator()()
	{
		complex_type tmp = prev;
		prev *= delta;
		prev /= abs(prev);
		return tmp;
	}
};

}
