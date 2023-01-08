/*
Bose Chaudhuri Hocquenghem Encoder

Copyright 2018 Ahmet Inan <inan@aicodix.de>
*/

#pragma once

#include <initializer_list>
#include "bitman.hh"

namespace CODE {

template <int LEN, int MSG>
class BoseChaudhuriHocquenghemEncoder
{
public:
	static const int N = LEN, K = MSG, NP = N - K;
	static const int G = ((NP+1)+7)/8;
private:
	uint8_t generator[G];
	static constexpr uint8_t slb1(uint8_t *buf, int pos)
	{
		return (buf[pos]<<1) | (buf[pos+1]>>7);
	}
public:
	BoseChaudhuriHocquenghemEncoder(std::initializer_list<int> minimal_polynomials)
	{
		// $generator(x) = \prod_i(minpoly_i(x))$
		int generator_degree = 1;
		for (int i = 0; i < G; ++i)
			generator[i] = 0;
		set_be_bit(generator, NP, 1);
		for (auto m: minimal_polynomials) {
			assert(0 < m);
			int m_degree = 0;
			while (m>>m_degree)
				++m_degree;
			--m_degree;
			assert(generator_degree + m_degree <= NP + 1);
			for (int i = generator_degree; i >= 0; --i) {
				if (!get_be_bit(generator, NP-i))
					continue;
				set_be_bit(generator, NP-i, m&1);
				for (int j = 1; j <= m_degree; ++j)
					xor_be_bit(generator, NP-(i+j), (m>>j)&1);
			}
			generator_degree += m_degree;
		}
		assert(generator_degree == NP + 1);
		if (0) {
			std::cerr << "generator =";
			for (int i = 0; i <= NP; ++i)
				std::cerr << " " << get_be_bit(generator, NP-i);
			std::cerr << std::endl;
		}
		for (int i = 0; i < NP; ++i)
			set_be_bit(generator, i, get_be_bit(generator, i+1));
		set_be_bit(generator, NP, 0);
	}
	void operator()(const uint8_t *data, uint8_t *parity, int data_len = K)
	{
		assert(0 < data_len && data_len <= K);
		// $code = data * x^{NP} + (data * x^{NP}) \mod{generator}$
		for (int l = 0; l <= (NP-1)/8; ++l)
			parity[l] = 0;
		for (int i = 0; i < data_len; ++i) {
			if (get_be_bit(data, i) != get_be_bit(parity, 0)) {
				for (int l = 0; l < (NP-1)/8; ++l)
					parity[l] = generator[l] ^ slb1(parity, l);
				parity[(NP-1)/8] = generator[(NP-1)/8] ^ (parity[(NP-1)/8]<<1);
			} else {
				for (int l = 0; l < (NP-1)/8; ++l)
					parity[l] = slb1(parity, l);
				parity[(NP-1)/8] <<= 1;
			}
		}
	}
};

template <int ROOTS, int FCR, int MSG, typename GF>
class BoseChaudhuriHocquenghemEncoderReference
{
public:
	typedef typename GF::value_type value_type;
	typedef typename GF::ValueType ValueType;
	typedef typename GF::IndexType IndexType;
	static const int NR = ROOTS;
	static const int N = GF::N, K = MSG, NP = N - K;
private:
	ValueType generator[NP+1];
public:
	BoseChaudhuriHocquenghemEncoderReference(std::initializer_list<int> minimal_polynomials)
	{
		// $generator(x) = \prod_i(minpoly_i(x))$
		int generator_degree = 1;
		generator[0] = ValueType(1);
		for (int i = 1; i <= NP; ++i)
			generator[i] = ValueType(0);
		for (auto m: minimal_polynomials) {
			assert(0 < m && m < 1<<(GF::M+1));
			int m_degree = GF::M;
			while (!(m>>m_degree))
				--m_degree;
			assert(generator_degree + m_degree <= NP + 1);
			for (int i = generator_degree; i >= 0; --i) {
				if (!generator[i])
					continue;
				generator[i] = ValueType(m&1);
				for (int j = 1; j <= m_degree; ++j)
					generator[i+j] += ValueType((m>>j)&1);
			}
			generator_degree += m_degree;
		}
		assert(generator_degree == NP + 1);
		if (0) {
			IndexType root(FCR), pe(1);
			for (int i = 0; i < NR; ++i) {
				ValueType tmp(generator[NP]);
				for (int j = 1; j <= NP; ++j)
					tmp = fma(root, tmp, generator[NP-j]);
				assert(!tmp);
				root *= pe;
			}
			std::cerr << "generator =";
			for (int i = 0; i <= NP; ++i)
				std::cerr << " " << (int)generator[i];
			std::cerr << std::endl;
		}
	}
	void operator()(const ValueType *data, ValueType *parity, int data_len = K)
	{
		assert(0 < data_len && data_len <= K);
		// $code = data * x^{NP} + (data * x^{NP}) \mod{generator}$
		for (int i = 0; i < NP; ++i)
			parity[i] = ValueType(0);
		for (int i = 0; i < data_len; ++i) {
			if (data[i] != parity[0]) {
				for (int j = 1; j < NP; ++j)
					parity[j-1] = generator[NP-j] + parity[j];
				parity[NP-1] = generator[0];
			} else {
				for (int j = 1; j < NP; ++j)
					parity[j-1] = parity[j];
				parity[NP-1] = ValueType(0);
			}
		}
	}
};

}

