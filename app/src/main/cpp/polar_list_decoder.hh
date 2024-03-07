/*
Successive cancellation list decoding of polar codes

Copyright 2020 Ahmet Inan <inan@aicodix.de>
*/

#pragma once

#include "sort.hh"
#include "polar_helper.hh"

namespace CODE {

template <typename TYPE, int M>
struct PolarListNode
{
	typedef PolarHelper<TYPE> PH;
	typedef typename PH::PATH PATH;
	typedef typename PH::MAP MAP;
	static const int N = 1 << M;
	static MAP rate0(PATH *metric, TYPE *hard, TYPE *soft)
	{
		for (int i = 0; i < N; ++i)
			hard[i] = PH::one();
		for (int i = 0; i < N; ++i)
			for (int k = 0; k < TYPE::SIZE; ++k)
				if (soft[i+N].v[k] < 0)
					metric[k] -= soft[i+N].v[k];
		MAP map;
		for (int k = 0; k < TYPE::SIZE; ++k)
			map.v[k] = k;
		return map;
	}
};

template <typename TYPE>
struct PolarListNode<TYPE, 0>
{
	typedef PolarHelper<TYPE> PH;
	typedef typename PH::PATH PATH;
	typedef typename PH::MAP MAP;
	static MAP rate0(PATH *metric, TYPE *hard, TYPE *soft)
	{
		*hard = PH::one();
		for (int k = 0; k < TYPE::SIZE; ++k)
			if (soft[1].v[k] < 0)
				metric[k] -= soft[1].v[k];
		MAP map;
		for (int k = 0; k < TYPE::SIZE; ++k)
			map.v[k] = k;
		return map;
	}
	static MAP rate1(PATH *metric, TYPE *message, MAP *maps, int *count, TYPE *hard, TYPE *soft)
	{
		TYPE sft = soft[1];
		PATH fork[2*TYPE::SIZE];
		for (int k = 0; k < TYPE::SIZE; ++k)
			fork[2*k] = fork[2*k+1] = metric[k];
		for (int k = 0; k < TYPE::SIZE; ++k)
			if (sft.v[k] < 0)
				fork[2*k] -= sft.v[k];
			else
				fork[2*k+1] += sft.v[k];
		int perm[2*TYPE::SIZE];
		CODE::insertion_sort(perm, fork, 2*TYPE::SIZE);
		for (int k = 0; k < TYPE::SIZE; ++k)
			metric[k] = fork[k];
		MAP map;
		for (int k = 0; k < TYPE::SIZE; ++k)
			map.v[k] = perm[k] >> 1;
		TYPE hrd;
		for (int k = 0; k < TYPE::SIZE; ++k)
			hrd.v[k] = 1 - 2 * (perm[k] & 1);
		message[*count] = hrd;
		maps[*count] = map;
		++*count;
		*hard = hrd;
		return map;
	}
};

template <typename TYPE, int M>
struct PolarListTree
{
	typedef PolarHelper<TYPE> PH;
	typedef typename PH::PATH PATH;
	typedef typename PH::MAP MAP;
	static const int N = 1 << M;
	static MAP decode(PATH *metric, TYPE *message, MAP *maps, int *count, TYPE *hard, TYPE *soft, const uint32_t *frozen)
	{
		for (int i = 0; i < N/2; ++i)
			soft[i+N/2] = PH::prod(soft[i+N], soft[i+N/2+N]);
		MAP lmap = PolarListTree<TYPE, M-1>::decode(metric, message, maps, count, hard, soft, frozen);
		for (int i = 0; i < N/2; ++i)
			soft[i+N/2] = PH::madd(hard[i], vshuf(soft[i+N], lmap), vshuf(soft[i+N/2+N], lmap));
		MAP rmap = PolarListTree<TYPE, M-1>::decode(metric, message, maps, count, hard+N/2, soft, frozen+N/2/32);
		for (int i = 0; i < N/2; ++i)
			hard[i] = PH::qmul(vshuf(hard[i], rmap), hard[i+N/2]);
		return vshuf(lmap, rmap);
	}
};

template <typename TYPE>
struct PolarListTree<TYPE, 6>
{
	typedef PolarHelper<TYPE> PH;
	typedef typename PH::PATH PATH;
	typedef typename PH::MAP MAP;
	static const int M = 6;
	static const int N = 1 << M;
	static MAP decode(PATH *metric, TYPE *message, MAP *maps, int *count, TYPE *hard, TYPE *soft, const uint32_t *frozen)
	{
		for (int i = 0; i < N/2; ++i)
			soft[i+N/2] = PH::prod(soft[i+N], soft[i+N/2+N]);
		MAP lmap, rmap;
		if (frozen[0] == 0xffffffff)
			lmap = PolarListNode<TYPE, M-1>::rate0(metric, hard, soft);
		else
			lmap = PolarListTree<TYPE, M-1>::decode(metric, message, maps, count, hard, soft, frozen[0]);
		for (int i = 0; i < N/2; ++i)
			soft[i+N/2] = PH::madd(hard[i], vshuf(soft[i+N], lmap), vshuf(soft[i+N/2+N], lmap));
		if (frozen[1] == 0xffffffff)
			rmap = PolarListNode<TYPE, M-1>::rate0(metric, hard+N/2, soft);
		else
			rmap = PolarListTree<TYPE, M-1>::decode(metric, message, maps, count, hard+N/2, soft, frozen[1]);
		for (int i = 0; i < N/2; ++i)
			hard[i] = PH::qmul(vshuf(hard[i], rmap), hard[i+N/2]);
		return vshuf(lmap, rmap);
	}
};

template <typename TYPE>
struct PolarListTree<TYPE, 5>
{
	typedef PolarHelper<TYPE> PH;
	typedef typename PH::PATH PATH;
	typedef typename PH::MAP MAP;
	static const int M = 5;
	static const int N = 1 << M;
	static MAP decode(PATH *metric, TYPE *message, MAP *maps, int *count, TYPE *hard, TYPE *soft, uint32_t frozen)
	{
		for (int i = 0; i < N/2; ++i)
			soft[i+N/2] = PH::prod(soft[i+N], soft[i+N/2+N]);
		MAP lmap, rmap;
		if ((frozen & ((1<<(1<<(M-1)))-1)) == ((1<<(1<<(M-1)))-1))
			lmap = PolarListNode<TYPE, M-1>::rate0(metric, hard, soft);
		else
			lmap = PolarListTree<TYPE, M-1>::decode(metric, message, maps, count, hard, soft, frozen & ((1<<(1<<(M-1)))-1));
		for (int i = 0; i < N/2; ++i)
			soft[i+N/2] = PH::madd(hard[i], vshuf(soft[i+N], lmap), vshuf(soft[i+N/2+N], lmap));
		if (frozen >> (N/2) == ((1<<(1<<(M-1)))-1))
			rmap = PolarListNode<TYPE, M-1>::rate0(metric, hard+N/2, soft);
		else
			rmap = PolarListTree<TYPE, M-1>::decode(metric, message, maps, count, hard+N/2, soft, frozen >> (N/2));
		for (int i = 0; i < N/2; ++i)
			hard[i] = PH::qmul(vshuf(hard[i], rmap), hard[i+N/2]);
		return vshuf(lmap, rmap);
	}
};

template <typename TYPE>
struct PolarListTree<TYPE, 4>
{
	typedef PolarHelper<TYPE> PH;
	typedef typename PH::PATH PATH;
	typedef typename PH::MAP MAP;
	static const int M = 4;
	static const int N = 1 << M;
	static MAP decode(PATH *metric, TYPE *message, MAP *maps, int *count, TYPE *hard, TYPE *soft, uint32_t frozen)
	{
		for (int i = 0; i < N/2; ++i)
			soft[i+N/2] = PH::prod(soft[i+N], soft[i+N/2+N]);
		MAP lmap, rmap;
		if ((frozen & ((1<<(1<<(M-1)))-1)) == ((1<<(1<<(M-1)))-1))
			lmap = PolarListNode<TYPE, M-1>::rate0(metric, hard, soft);
		else
			lmap = PolarListTree<TYPE, M-1>::decode(metric, message, maps, count, hard, soft, frozen & ((1<<(1<<(M-1)))-1));
		for (int i = 0; i < N/2; ++i)
			soft[i+N/2] = PH::madd(hard[i], vshuf(soft[i+N], lmap), vshuf(soft[i+N/2+N], lmap));
		if (frozen >> (N/2) == ((1<<(1<<(M-1)))-1))
			rmap = PolarListNode<TYPE, M-1>::rate0(metric, hard+N/2, soft);
		else
			rmap = PolarListTree<TYPE, M-1>::decode(metric, message, maps, count, hard+N/2, soft, frozen >> (N/2));
		for (int i = 0; i < N/2; ++i)
			hard[i] = PH::qmul(vshuf(hard[i], rmap), hard[i+N/2]);
		return vshuf(lmap, rmap);
	}
};

template <typename TYPE>
struct PolarListTree<TYPE, 3>
{
	typedef PolarHelper<TYPE> PH;
	typedef typename PH::PATH PATH;
	typedef typename PH::MAP MAP;
	static const int M = 3;
	static const int N = 1 << M;
	static MAP decode(PATH *metric, TYPE *message, MAP *maps, int *count, TYPE *hard, TYPE *soft, uint32_t frozen)
	{
		for (int i = 0; i < N/2; ++i)
			soft[i+N/2] = PH::prod(soft[i+N], soft[i+N/2+N]);
		MAP lmap, rmap;
		if ((frozen & ((1<<(1<<(M-1)))-1)) == ((1<<(1<<(M-1)))-1))
			lmap = PolarListNode<TYPE, M-1>::rate0(metric, hard, soft);
		else
			lmap = PolarListTree<TYPE, M-1>::decode(metric, message, maps, count, hard, soft, frozen & ((1<<(1<<(M-1)))-1));
		for (int i = 0; i < N/2; ++i)
			soft[i+N/2] = PH::madd(hard[i], vshuf(soft[i+N], lmap), vshuf(soft[i+N/2+N], lmap));
		if (frozen >> (N/2) == ((1<<(1<<(M-1)))-1))
			rmap = PolarListNode<TYPE, M-1>::rate0(metric, hard+N/2, soft);
		else
			rmap = PolarListTree<TYPE, M-1>::decode(metric, message, maps, count, hard+N/2, soft, frozen >> (N/2));
		for (int i = 0; i < N/2; ++i)
			hard[i] = PH::qmul(vshuf(hard[i], rmap), hard[i+N/2]);
		return vshuf(lmap, rmap);
	}
};

template <typename TYPE>
struct PolarListTree<TYPE, 2>
{
	typedef PolarHelper<TYPE> PH;
	typedef typename PH::PATH PATH;
	typedef typename PH::MAP MAP;
	static const int M = 2;
	static const int N = 1 << M;
	static MAP decode(PATH *metric, TYPE *message, MAP *maps, int *count, TYPE *hard, TYPE *soft, uint32_t frozen)
	{
		for (int i = 0; i < N/2; ++i)
			soft[i+N/2] = PH::prod(soft[i+N], soft[i+N/2+N]);
		MAP lmap, rmap;
		if ((frozen & ((1<<(1<<(M-1)))-1)) == ((1<<(1<<(M-1)))-1))
			lmap = PolarListNode<TYPE, M-1>::rate0(metric, hard, soft);
		else
			lmap = PolarListTree<TYPE, M-1>::decode(metric, message, maps, count, hard, soft, frozen & ((1<<(1<<(M-1)))-1));
		for (int i = 0; i < N/2; ++i)
			soft[i+N/2] = PH::madd(hard[i], vshuf(soft[i+N], lmap), vshuf(soft[i+N/2+N], lmap));
		if (frozen >> (N/2) == ((1<<(1<<(M-1)))-1))
			rmap = PolarListNode<TYPE, M-1>::rate0(metric, hard+N/2, soft);
		else
			rmap = PolarListTree<TYPE, M-1>::decode(metric, message, maps, count, hard+N/2, soft, frozen >> (N/2));
		for (int i = 0; i < N/2; ++i)
			hard[i] = PH::qmul(vshuf(hard[i], rmap), hard[i+N/2]);
		return vshuf(lmap, rmap);
	}
};

template <typename TYPE>
struct PolarListTree<TYPE, 1>
{
	typedef PolarHelper<TYPE> PH;
	typedef typename PH::PATH PATH;
	typedef typename PH::MAP MAP;
	static MAP decode(PATH *metric, TYPE *message, MAP *maps, int *count, TYPE *hard, TYPE *soft, uint32_t frozen)
	{
		soft[1] = PH::prod(soft[2], soft[3]);
		MAP lmap, rmap;
		if (frozen & 1)
			lmap = PolarListNode<TYPE, 0>::rate0(metric, hard, soft);
		else
			lmap = PolarListNode<TYPE, 0>::rate1(metric, message, maps, count, hard, soft);
		soft[1] = PH::madd(hard[0], vshuf(soft[2], lmap), vshuf(soft[3], lmap));
		if (frozen >> 1)
			rmap = PolarListNode<TYPE, 0>::rate0(metric, hard+1, soft);
		else
			rmap = PolarListNode<TYPE, 0>::rate1(metric, message, maps, count, hard+1, soft);
		hard[0] = PH::qmul(vshuf(hard[0], rmap), hard[1]);
		return vshuf(lmap, rmap);
	}
};

template <typename TYPE, int MAX_M>
class PolarListDecoder
{
	static_assert(MAX_M >= 5 && MAX_M <= 16);
	typedef PolarHelper<TYPE> PH;
	typedef typename TYPE::value_type VALUE;
	typedef typename PH::PATH PATH;
	typedef typename PH::MAP MAP;
	static const int MAX_N = 1 << MAX_M;
	TYPE soft[2*MAX_N];
	TYPE hard[MAX_N];
	MAP maps[MAX_N];
public:
	void operator()(int *rank, TYPE *message, const VALUE *codeword, const uint32_t *frozen, int level)
	{
		assert(level <= MAX_M);
		PATH metric[TYPE::SIZE];
		int count = 0;
		metric[0] = 0;
		for (int k = 1; k < TYPE::SIZE; ++k)
			metric[k] = 1000000;
		int length = 1 << level;
		for (int i = 0; i < length; ++i)
			soft[length+i] = vdup<TYPE>(codeword[i]);

		switch (level) {
		case 5: PolarListTree<TYPE, 5>::decode(metric, message, maps, &count, hard, soft, *frozen); break;
		case 6: PolarListTree<TYPE, 6>::decode(metric, message, maps, &count, hard, soft, frozen); break;
		case 7: PolarListTree<TYPE, 7>::decode(metric, message, maps, &count, hard, soft, frozen); break;
		case 8: PolarListTree<TYPE, 8>::decode(metric, message, maps, &count, hard, soft, frozen); break;
		case 9: PolarListTree<TYPE, 9>::decode(metric, message, maps, &count, hard, soft, frozen); break;
		case 10: PolarListTree<TYPE, 10>::decode(metric, message, maps, &count, hard, soft, frozen); break;
		case 11: PolarListTree<TYPE, 11>::decode(metric, message, maps, &count, hard, soft, frozen); break;
		case 12: PolarListTree<TYPE, 12>::decode(metric, message, maps, &count, hard, soft, frozen); break;
		case 13: PolarListTree<TYPE, 13>::decode(metric, message, maps, &count, hard, soft, frozen); break;
		case 14: PolarListTree<TYPE, 14>::decode(metric, message, maps, &count, hard, soft, frozen); break;
		case 15: PolarListTree<TYPE, 15>::decode(metric, message, maps, &count, hard, soft, frozen); break;
		case 16: PolarListTree<TYPE, 16>::decode(metric, message, maps, &count, hard, soft, frozen); break;
		default: assert(false);
		}

		for (int i = 0, r = 0; rank != nullptr && i < TYPE::SIZE; ++i) {
			if (i > 0 && metric[i-1] != metric[i])
				++r;
			rank[i] = r;
		}
		MAP acc = maps[count-1];
		for (int i = count-2; i >= 0; --i) {
			message[i] = vshuf(message[i], acc);
			acc = vshuf(maps[i], acc);
		}
	}
};

}

