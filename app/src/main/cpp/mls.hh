/*
Maximum length sequence

Copyright 2020 Ahmet Inan <inan@aicodix.de>
*/

#pragma once

namespace CODE {

class MLS
{
	static int hibit(unsigned n)
	{
		n |= n >> 1;
		n |= n >> 2;
		n |= n >> 4;
		n |= n >> 8;
		n |= n >> 16;
		return n ^ (n >> 1);
	}
	int poly, test, reg;
public:
	MLS(int poly = 0b100000000000000001001, int reg = 1) : poly(poly), test(hibit(poly)>>1), reg(reg) {}
	void reset(int r = 1)
	{
		reg = r;
	}
	bool bad(int r = 1)
	{
		reg = r;
		int len = hibit(poly) - 1;
		for (int i = 1; i < len; ++i) {
			(*this)();
			if (reg == r)
				return true;
		}
		(*this)();
		return reg != r;
	}
	bool operator()()
	{
		bool fb = reg & test;
		reg <<= 1;
		reg ^= fb * poly;
		return fb;
	}
};

}

