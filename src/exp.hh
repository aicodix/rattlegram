/*
Exponentiation approximations

Constants below lifted from the Cephes Mathematical Library:
https://www.netlib.org/cephes/cmath.tgz

Copyright 2018 Ahmet Inan <inan@aicodix.de>
*/

#pragma once

namespace DSP {

template <typename TYPE>
TYPE ldexp(TYPE x, int n)
{
	int a = n < 0 ? -n : n;
	int a8 = a / 8;
	int ar = a - a8 * 8;
	TYPE t = 1 << a8;
	t *= t;
	t *= t;
	t *= t;
	t *= 1 << ar;
	return n < 0 ? x / t : x * t;
}

template <typename TYPE>
TYPE exp10(TYPE x)
{
	static constexpr TYPE
		LOG210 = 3.32192809488736234787e0,
		LG102A = 3.01025390625000000000E-1,
		LG102B = 4.60503898119521373889E-6,
		P0 = 4.09962519798587023075E-2,
		P1 = 1.17452732554344059015E1,
		P2 = 4.06717289936872725516E2,
		P3 = 2.39423741207388267439E3,
		Q0 = 8.50936160849306532625E1,
		Q1 = 1.27209271178345121210E3,
		Q2 = 2.07960819286001865907E3;
	TYPE i = nearbyint(x * LOG210);
	x -= i * LG102A;
	x -= i * LG102B;
	TYPE xx = x * x;
	TYPE py = x * (xx * (xx * (xx * P0 + P1) + P2) + P3);
	TYPE qy = xx * (xx * (xx + Q0) + Q1) + Q2;
	TYPE pq = TYPE(1) + TYPE(2) * py / (qy - py);
	return ldexp(pq, (int)i);
}

}

