/*
SIMD wrapper used by polar encoder and decoder

Copyright 2020 Ahmet Inan <inan@aicodix.de>
*/

#pragma once

#include "simd.hh"

namespace CODE {

template <typename TYPE>
struct PolarHelper
{
	typedef TYPE PATH;
	static TYPE one()
	{
		return 1;
	}
	static TYPE zero()
	{
		return 0;
	}
	static TYPE signum(TYPE v)
	{
		return (v > 0) - (v < 0);
	}
	template <typename IN>
	static TYPE quant(IN in)
	{
		return in;
	}
	static TYPE qabs(TYPE a)
	{
		return std::abs(a);
	}
	static TYPE qmin(TYPE a, TYPE b)
	{
		return std::min(a, b);
	}
	static TYPE qadd(TYPE a, TYPE b)
	{
		return a + b;
	}
	static TYPE qmul(TYPE a, TYPE b)
	{
		return a * b;
	}
	static TYPE prod(TYPE a, TYPE b)
	{
		return signum(a) * signum(b) * qmin(qabs(a), qabs(b));
	}
	static TYPE madd(TYPE a, TYPE b, TYPE c)
	{
		return a * b + c;
	}
};

template <typename VALUE, int WIDTH>
struct PolarHelper<SIMD<VALUE, WIDTH>>
{
	typedef SIMD<VALUE, WIDTH> TYPE;
	typedef VALUE PATH;
	typedef SIMD<typename TYPE::uint_type, WIDTH> MAP;
	static TYPE one()
	{
		return vdup<TYPE>(1);
	}
	static TYPE zero()
	{
		return vzero<TYPE>();
	}
	static TYPE signum(TYPE a)
	{
		return vsignum(a);
	}
	static TYPE qabs(TYPE a)
	{
		return vabs(a);
	}
	static TYPE qmin(TYPE a, TYPE b)
	{
		return vmin(a, b);
	}
	static TYPE qadd(TYPE a, TYPE b)
	{
		return vadd(a, b);
	}
	static TYPE qmul(TYPE a, TYPE b)
	{
		return vmul(a, b);
	}
	static TYPE prod(TYPE a, TYPE b)
	{
		return vmul(vmul(vsignum(a), vsignum(b)), vmin(vabs(a), vabs(b)));
	}
	static TYPE madd(TYPE a, TYPE b, TYPE c)
	{
		return vadd(vmul(a, b), c);
	}
};

template <int WIDTH>
struct PolarHelper<SIMD<int8_t, WIDTH>>
{
	typedef SIMD<int8_t, WIDTH> TYPE;
	typedef int PATH;
	typedef SIMD<uint8_t, WIDTH> MAP;
	static TYPE one()
	{
		return vdup<TYPE>(1);
	}
	static TYPE zero()
	{
		return vzero<TYPE>();
	}
	static TYPE signum(TYPE a)
	{
		return vsignum(a);
	}
	static TYPE qabs(TYPE a)
	{
		return vqabs(a);
	}
	static TYPE qadd(TYPE a, TYPE b)
	{
		return vqadd(a, b);
	}
	static TYPE qmul(TYPE a, TYPE b)
	{
#ifdef __ARM_NEON
		return vmul(a, b);
#else
		return vsign(a, b);
#endif
	}
	static TYPE prod(TYPE a, TYPE b)
	{
#ifdef __ARM_NEON
		return vmul(vmul(vsignum(a), vsignum(b)), vmin(vqabs(a), vqabs(b)));
#else
		return vsign(vmin(vqabs(a), vqabs(b)), vsign(vsignum(a), b));
#endif
	}
	static TYPE madd(TYPE a, TYPE b, TYPE c)
	{
#ifdef __ARM_NEON
		return vmax(vqadd(vmul(a, vmax(b, vdup<TYPE>(-127))), c), vdup<TYPE>(-127));
#else
		return vmax(vqadd(vsign(vmax(b, vdup<TYPE>(-127)), a), c), vdup<TYPE>(-127));
#endif
	}
};

template <>
struct PolarHelper<int8_t>
{
	typedef int PATH;
	static int8_t one()
	{
		return 1;
	}
	static int8_t zero()
	{
		return 0;
	}
	static int8_t signum(int8_t v)
	{
		return (v > 0) - (v < 0);
	}
	template <typename IN>
	static int8_t quant(IN in)
	{
		return std::min<IN>(std::max<IN>(std::nearbyint(in), -127), 127);
	}
	static int8_t qabs(int8_t a)
	{
		return std::abs(std::max<int8_t>(a, -127));
	}
	static int8_t qmin(int8_t a, int8_t b)
	{
		return std::min(a, b);
	}
	static int8_t qadd(int8_t a, int8_t b)
	{
		return std::min<int16_t>(std::max<int16_t>(int16_t(a) + int16_t(b), -127), 127);
	}
	static int8_t qmul(int8_t a, int8_t b)
	{
		// return std::min<int16_t>(std::max<int16_t>(int16_t(a) * int16_t(b), -127), 127);
		// only used for hard decision values anyway
		return a * b;
	}
	static int8_t prod(int8_t a, int8_t b)
	{
		return signum(a) * signum(b) * qmin(qabs(a), qabs(b));
	}
	static int8_t madd(int8_t a, int8_t b, int8_t c)
	{
		return std::min<int16_t>(std::max<int16_t>(int16_t(a) * int16_t(b) + int16_t(c), -127), 127);
	}
};

}

