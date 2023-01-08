/*
Trigonometric functions on the unit circle

Constants below lifted from the Cephes Mathematical Library:
https://www.netlib.org/cephes/cmath.tgz

Copyright 2019 Ahmet Inan <inan@aicodix.de>
*/

#pragma once

#include "const.hh"

namespace DSP {

template <typename TYPE>
class UnitCircle
{
	static constexpr TYPE
		a1 = 1.0,
		a2 = -0.5,
		a3 = -1.66666666666666307295E-1,
		a4 = 4.16666666666665929218E-2,
		a5 = 8.33333333332211858878E-3,
		a6 = -1.38888888888730564116E-3,
		a7 = -1.98412698295895385996E-4,
		a8 = 2.48015872888517045348E-5,
		a9 = 2.75573136213857245213E-6,
		a10 = -2.75573141792967388112E-7,
		a11 = -2.50507477628578072866E-8,
		a12 = 2.08757008419747316778E-9,
		a13 = 1.58962301576546568060E-10,
		a14 = -1.13585365213876817300E-11;

	static constexpr TYPE cosine_kernel(TYPE xx)
	{
		return a1 + xx * (xx * (xx * (xx * (xx * (xx * (xx * a14 + a12) + a10) + a8) + a6) + a4) + a2);
	}
	static constexpr TYPE sine_kernel(TYPE x, TYPE xx)
	{
		return x * (xx * (xx * (xx * (xx * (xx * (xx * a13 + a11) + a9) + a7) + a5) + a3) + a1);
	}
	static constexpr TYPE cosine_approx(TYPE x)
	{
		return cosine_kernel(x * x);
	}
	static constexpr TYPE sine_approx(TYPE x)
	{
		return sine_kernel(x, x * x);
	}
public:
	static constexpr TYPE rad(int n, int N)
	{
		return Const<TYPE>::TwoPi() * TYPE(n) / TYPE(N);
	}
	static constexpr TYPE cos(int n, int N)
	{
		return
			8*n <-7*N ? cosine_approx(rad(n+N, N)) :
			8*n <-5*N ? -sine_approx(rad(4*n+3*N, 4*N)) :
			8*n <-3*N ? -cosine_approx(rad(2*n+N, 2*N)) :
			8*n <-1*N ? sine_approx(rad(4*n+N, 4*N)) :
			8*n < 1*N ? cosine_approx(rad(n, N)) :
			8*n < 3*N ? -sine_approx(rad(4*n-N, 4*N)) :
			8*n < 5*N ? -cosine_approx(rad(2*n-N, 2*N)) :
			8*n < 7*N ? sine_approx(rad(4*n-3*N, 4*N)) :
			cosine_approx(rad(n-N, N));
	}
	static constexpr TYPE sin(int n, int N)
	{
		return
			8*n <-7*N ? sine_approx(rad(n+N, N)) :
			8*n <-5*N ? cosine_approx(rad(4*n+3*N, 4*N)) :
			8*n <-3*N ? -sine_approx(rad(2*n+N, 2*N)) :
			8*n <-1*N ? -cosine_approx(rad(4*n+N, 4*N)) :
			8*n < 1*N ? sine_approx(rad(n, N)) :
			8*n < 3*N ? cosine_approx(rad(4*n-N, 4*N)) :
			8*n < 5*N ? -sine_approx(rad(2*n-N, 2*N)) :
			8*n < 7*N ? -cosine_approx(rad(4*n-3*N, 4*N)) :
			sine_approx(rad(n-N, N));
	}
};

}

