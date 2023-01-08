/*
Class of pseudorandom number generators, discovered by George Marsaglia

Copyright 2018 Ahmet Inan <inan@aicodix.de>
*/

#pragma once

namespace CODE {

class Xorshift32
{
	static const uint32_t Y = 2463534242;
	uint32_t y_;
public:
	typedef uint32_t result_type;
	static constexpr result_type min()
	{
		return 0;
	}
	static constexpr result_type max()
	{
		return UINT32_MAX;
	}
	Xorshift32(uint32_t y = Y) : y_(y) {}
	void reset(uint32_t y = Y)
	{
		y_ = y;
	}
	uint32_t operator()()
	{
		y_ ^= y_ << 13;
		y_ ^= y_ >> 17;
		y_ ^= y_ << 5;
		return y_;
	}
};

class Xorshift64
{
	static const uint64_t X = 88172645463325252;
	uint64_t x_;
public:
	typedef uint64_t result_type;
	static constexpr result_type min()
	{
		return 0;
	}
	static constexpr result_type max()
	{
		return UINT64_MAX;
	}
	Xorshift64(uint64_t x = X) : x_(x) {}
	void reset(uint64_t x = X)
	{
		x_ = x;
	}
	uint64_t operator()()
	{
		x_ ^= x_ << 13;
		x_ ^= x_ >> 7;
		x_ ^= x_ << 17;
		return x_;
	}
};

class Xorwow
{
	static const uint32_t X = 123456789;
	static const uint32_t Y = 362436069;
	static const uint32_t Z = 521288629;
	static const uint32_t W = 88675123;
	static const uint32_t V = 5783321;
	static const uint32_t D = 6615241;
	uint32_t x_, y_, z_, w_, v_, d_;
public:
	typedef uint32_t result_type;
	static constexpr result_type min()
	{
		return 0;
	}
	static constexpr result_type max()
	{
		return UINT32_MAX;
	}
	Xorwow(uint32_t x = X, uint32_t y = Y,
		uint32_t z = Z, uint32_t w = W,
		uint32_t v = V, uint32_t d = D) :
		x_(x), y_(y), z_(z), w_(w), v_(v), d_(d) {}
	void reset(uint32_t x = X, uint32_t y = Y,
		uint32_t z = Z, uint32_t w = W,
		uint32_t v = V, uint32_t d = D)
	{
		x_ = x;
		y_ = y;
		z_ = z;
		w_ = w;
		v_ = v;
		d_ = d;
	}
	uint32_t operator()()
	{
		uint32_t t = x_ ^ (x_ >> 2);
		x_ = y_; y_ = z_; z_ = w_; w_ = v_;
		v_ = (v_ ^ (v_ << 4)) ^ (t ^ (t << 1));
		d_ += 362437;
		return d_ + v_;
	}
};

class Xorshift128
{
	static const uint32_t X = 123456789;
	static const uint32_t Y = 362436069;
	static const uint32_t Z = 521288629;
	static const uint32_t W = 88675123;
	uint32_t x_, y_, z_, w_;
public:
	typedef uint32_t result_type;
	static constexpr result_type min()
	{
		return 0;
	}
	static constexpr result_type max()
	{
		return UINT32_MAX;
	}
	Xorshift128(uint32_t x = X, uint32_t y = Y,
		uint32_t z = Z, uint32_t w = W) :
		x_(x), y_(y), z_(z), w_(w) {}
	void reset(uint32_t x = X, uint32_t y = Y,
		uint32_t z = Z, uint32_t w = W)
	{
		x_ = x;
		y_ = y;
		z_ = z;
		w_ = w;
	}
	uint32_t operator()()
	{
		uint32_t t = (x_ ^ (x_ << 11));
		x_ = y_; y_ = z_; z_ = w_;
		w_ = (w_ ^ (w_ >> 19)) ^ (t ^ (t >> 8));
		return w_;
	}
};

}

