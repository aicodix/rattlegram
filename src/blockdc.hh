/*
DC Blocker

Copyright 2019 Ahmet Inan <inan@aicodix.de>
*/

#pragma once

namespace DSP {

template <typename TYPE, typename VALUE>
class BlockDC
{
	TYPE x1, y1;
	VALUE a, b;
public:
	constexpr BlockDC() : x1(0), y1(0), a(0), b(0.5)
	{
	}
	void samples(int s)
	{
		a = VALUE(s - 1) / VALUE(s);
		b = (VALUE(1) + a) / VALUE(2);
	}
	TYPE operator()(TYPE x0)
	{
		TYPE y0 = b * (x0 - x1) + a * y1;
		x1 = x0; y1 = y0;
		return y0;
	}
};

}

