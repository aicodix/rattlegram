/*
Kahan summation algorithm

Copyright 2018 Ahmet Inan <inan@aicodix.de>
*/

#pragma once

namespace DSP {

template <typename T>
class Kahan
{
	T high, low;
public:
	Kahan() : high(0), low(0) {}
	Kahan(T init) : high(init), low(0) {}
	Kahan(const Kahan &a) : high(a.high), low(a.low) {}
#if __clang__
	[[clang::optnone]]
#elif __GNUC__
	[[gnu::optimize("no-associative-math")]]
#else
#error unsupported compiler
#endif
	T operator ()(T input)
	{
		T tmp = input - low;
		T sum = high + tmp;
		low = (sum - high) - tmp;
		return high = sum;
	}
	bool operator == (const Kahan &a)
	{
		return low == a.low && high == a.high;
	}
	bool same(T input)
	{
		Kahan tmp(*this);
		(*this)(input);
		return tmp == *this;
	}
	T operator ()() { return high; }
};

template <class I, class T>
T kahan_sum(I begin, I end, T init)
{
	Kahan<T> kahan(init);
	while (begin != end)
		kahan(*begin++);
	return kahan();
}

}

