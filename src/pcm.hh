/*
Interface for reading and writing PCM data

Copyright 2018 Ahmet Inan <inan@aicodix.de>
*/

#pragma once

namespace DSP {

template <typename TYPE>
struct WritePCM
{
	virtual void write(const TYPE *, int, int = -1) = 0;
	virtual bool good() = 0;
	virtual void silence(int) = 0;
	virtual int rate() = 0;
	virtual int channels() = 0;
	virtual ~WritePCM() = default;
};

template <typename TYPE>
struct ReadPCM
{
	virtual void read(TYPE *, int, int = -1) = 0;
	virtual bool good() = 0;
	virtual void skip(int) = 0;
	virtual int rate() = 0;
	virtual int channels() = 0;
	virtual int frames() { return -1; }
	virtual ~ReadPCM() = default;
};

}

