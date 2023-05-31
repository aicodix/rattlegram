/*
Cyclic redundancy check

Copyright 2018 Ahmet Inan <inan@aicodix.de>
*/

#pragma once

namespace CODE {

template <typename TYPE>
class CRC
{
	TYPE lut[256];
	TYPE poly;
	TYPE crc;
	TYPE update(TYPE prev, bool data)
	{
		TYPE tmp = prev ^ data;
		return (prev >> 1) ^ ((tmp & 1) * poly);
	}
public:
	CRC(TYPE poly, TYPE crc = 0) : poly(poly), crc(crc)
	{
		for (int j = 0; j < 256; ++j) {
			TYPE tmp = j;
			for (int i = 8; i; --i)
				tmp = update(tmp, 0);
			lut[j] = tmp;
		}
	}
	void reset(TYPE v = 0)
	{
		crc = v;
	}
	TYPE operator()()
	{
		return crc;
	}
	TYPE operator()(bool data)
	{
		return crc = update(crc, data);
	}
	TYPE operator()(uint8_t data)
	{
		TYPE tmp = crc ^ data;
		return crc = (crc >> 8) ^ lut[tmp & 255];
	}

	TYPE operator()(uint16_t data)
	{
		(*this)(uint8_t(data & 255));
		(*this)(uint8_t((data >> 8) & 255));
		return crc;
	}
	TYPE operator()(uint32_t data)
	{
		(*this)(uint8_t(data & 255));
		(*this)(uint8_t((data >> 8) & 255));
		(*this)(uint8_t((data >> 16) & 255));
		(*this)(uint8_t((data >> 24) & 255));
		return crc;
	}
	TYPE operator()(uint64_t data)
	{
		(*this)(uint8_t(data & 255));
		(*this)(uint8_t((data >> 8) & 255));
		(*this)(uint8_t((data >> 16) & 255));
		(*this)(uint8_t((data >> 24) & 255));
		(*this)(uint8_t((data >> 32) & 255));
		(*this)(uint8_t((data >> 40) & 255));
		(*this)(uint8_t((data >> 48) & 255));
		(*this)(uint8_t((data >> 56) & 255));
		return crc;
	}
};

template<>
uint8_t CRC<uint8_t>::operator()(uint8_t data)
{
	return crc = lut[crc ^ data];
}

}

