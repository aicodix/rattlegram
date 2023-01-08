/*
Simple moving average

Copyright 2019 Ahmet Inan <inan@aicodix.de>
*/

#pragma once

#include "kahan.hh"
#include "delay.hh"
#include "swa.hh"

namespace DSP {

template <typename TYPE, typename VALUE, int NUM>
class SMA1
{
	TYPE hist_inp[NUM];
	TYPE hist_avg;
	int hist_pos;
public:
	SMA1() : hist_avg(0), hist_pos(0)
	{
		for (int i = 0; i < NUM; ++i)
			hist_inp[i] = 0;
	}
	VALUE abs_dev()
	{
		VALUE sum(abs(hist_inp[0] - hist_avg));
		for (int i = 1; i < NUM; ++i)
			sum += abs(hist_inp[i] - hist_avg);
		return sum / VALUE(NUM);
	}
	TYPE operator () (TYPE input)
	{
		hist_inp[hist_pos] = input;
		if (++hist_pos >= NUM)
			hist_pos = 0;
		TYPE hist_sum(hist_inp[0]);
		for (int i = 1; i < NUM; ++i)
			hist_sum += hist_inp[i];
		return hist_avg = hist_sum / VALUE(NUM);
	}
};

template <typename TYPE, typename VALUE, int NUM, bool NORM = true>
class SMA2
{
	Delay<TYPE, NUM> delay;
	TYPE sum;
public:
	SMA2() : sum(0)
	{
	}
	TYPE operator () (TYPE input)
	{
		sum += input - delay(input);
		if (NORM)
			return sum / VALUE(NUM);
		return sum;
	}
};

template <typename TYPE, typename VALUE, int NUM, bool NORM = true>
class SMA3
{
	Delay<TYPE, NUM> delay;
	Kahan<TYPE> sum;
public:
	SMA3() : sum(0)
	{
	}
	TYPE operator () (TYPE input)
	{
		sum(-delay(input));
		if (NORM)
			return sum(input) / VALUE(NUM);
		return sum(input);
	}
};

template <typename TYPE, typename VALUE, int NUM, bool NORM = true>
class SMA4
{
	SWA<TYPE, std::plus<TYPE>, NUM> swa;
public:
	SMA4() : swa(0)
	{
	}
	TYPE operator () (TYPE input)
	{
		if (NORM)
			return swa(input) / VALUE(NUM);
		return swa(input);
	}
};

}

