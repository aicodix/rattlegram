/*
Digital delay line

Copyright 2020 Ahmet Inan <inan@aicodix.de>
*/

#pragma once

namespace DSP {

template <typename TYPE, int NUM>
class Delay
{
	TYPE buf[NUM];
	int pos;
public:
	Delay(TYPE init = 0) : pos(0)
	{
		for (int i = 0; i < NUM; ++i)
			buf[i] = init;
	}
	TYPE operator () (TYPE input)
	{
		TYPE tmp = buf[pos];
		buf[pos] = input;
		if (++pos >= NUM)
			pos = 0;
		return tmp;
	}
};

}

