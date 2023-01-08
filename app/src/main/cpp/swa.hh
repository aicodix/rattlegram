/*
Sliding window accelerator

Copyright 2020 Ahmet Inan <inan@aicodix.de>
*/

#pragma once

namespace DSP {

template <typename TYPE, typename OP, int NUM>
class SWA
{
	TYPE tree[2 * NUM];
	int leaf;
	OP op;
public:
	SWA(TYPE ident) : leaf(NUM)
	{
		for (int i = 0; i < 2 * NUM; ++i)
			tree[i] = ident;
	}
	TYPE operator () (TYPE input)
	{
		tree[leaf] = input;
		for (int child = leaf, parent = leaf / 2; parent; child = parent, parent /= 2)
			tree[parent] = op(tree[child], tree[child^1]);
		if (++leaf >= 2 * NUM)
			leaf = NUM;
		return tree[1];
	}
};

}

