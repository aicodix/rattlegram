/*
Decibel calculation helpers

Copyright 2018 Ahmet Inan <inan@aicodix.de>
*/

#pragma once

#include "exp.hh"

namespace DSP {

template <typename TYPE>
TYPE decibel(TYPE v)
{
#if 0
	return TYPE(10) * log10(v);
#else
	static constexpr TYPE
		// 2*1024*10/l(10)
		scale = 8894.350989378597430295120259412072085389250678859091274024013479236L,
		inv1 = TYPE(1) / TYPE(1),
		inv3 = TYPE(1) / TYPE(3),
		inv5 = TYPE(1) / TYPE(5),
		inv7 = TYPE(1) / TYPE(7),
		inv9 = TYPE(1) / TYPE(9),
		inv11 = TYPE(1) / TYPE(11),
		inv13 = TYPE(1) / TYPE(13);
	TYPE x = v;
	for (int i = 0; i < 10; ++i)
		x = sqrt(x);
	x = (x - TYPE(1)) / (x + TYPE(1));
	TYPE xx = x * x;
	return scale * x * (xx * (xx * (xx * (xx * (xx * (xx *
		inv13 + inv11) + inv9) + inv7) + inv5) + inv3) + inv1);
#endif
}

template <typename TYPE>
TYPE idecibel(TYPE dB)
{
	return exp10(dB / TYPE(10));
}

}

