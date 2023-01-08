/*
Some little helpers

Copyright 2018 Ahmet Inan <inan@aicodix.de>
*/

#pragma once

#include "const.hh"

namespace DSP {

template <typename TYPE>
int signum(TYPE v)
{
	return (v > TYPE(0)) - (v < TYPE(0));
}

template <typename AB, typename X>
AB lerp(AB a, AB b, X x)
{
	return (X(1) - x) * a + x * b;
}

template <typename TYPE>
TYPE clamp(TYPE x, TYPE a, TYPE b)
{
	return x < a ? a : x > b ? b : x;
}

template <typename TYPE>
TYPE normal_pdf(TYPE x, TYPE m, TYPE s)
{
	return exp(-pow((x - m) / s, TYPE(2)) / TYPE(2)) / (Const<TYPE>::SqrtTwoPi() * s);
}

template <typename TYPE>
TYPE sinc(TYPE x)
{
	return TYPE(0) == x ? TYPE(1) : sin(Const<TYPE>::Pi() * x) / (Const<TYPE>::Pi() * x);
}

template <typename TYPE>
TYPE delta(TYPE x)
{
	return TYPE(0) == x ? TYPE(1) : TYPE(0);
}

}

