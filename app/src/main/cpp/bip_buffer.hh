/*
Bip buffer

Copyright 2020 Ahmet Inan <inan@aicodix.de>
*/

#pragma once

namespace DSP {

template <typename TYPE, int NUM>
class BipBuffer
{
	TYPE buf[2*NUM];
	int pos0, pos1;
public:
	BipBuffer() : pos0(0), pos1(NUM)
	{
		for (int i = 0; i < 2*NUM; ++i)
			buf[i] = 0;
	}
	const TYPE *operator () ()
	{
		return buf + min(pos0, pos1);
	}
	const TYPE *operator () (TYPE input)
	{
		buf[pos0] = buf[pos1] = input;
		if (++pos0 >= 2*NUM)
			pos0 = 0;
		if (++pos1 >= 2*NUM)
			pos1 = 0;
		return operator () ();
	}
};

}

