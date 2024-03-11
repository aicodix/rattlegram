/*
Quick algorithms for sorting and selecting

Copyright 2024 Ahmet Inan <inan@aicodix.de>
*/

#pragma once

namespace DSP {

namespace QUICK {

template <typename TYPE>
static inline void swap(TYPE *a, int i, int j)
{
	TYPE t = a[i];
	a[i] = a[j];
	a[j] = t;
}

template <typename TYPE>
static inline void median(TYPE *A, int a, int b, int c, int d, int e)
{
	if (A[c] < A[a])
		swap(A, a, c);
	if (A[d] < A[b])
		swap(A, b, d);
	if (A[d] < A[c]) {
		swap(A, c, d);
		swap(A, a, b);
	}
	if (A[e] < A[b])
		swap(A, b, e);
	if (A[e] < A[c]) {
		swap(A, c, e);
		if (A[c] < A[a])
			swap(A, a, c);
	} else if (A[c] < A[b]) {
		swap(A, b, c);
	}
}

template <typename TYPE>
static void partition(TYPE *a, int &l, int &h)
{
	int half = (h - l) / 2;
	int quarter = (h - l) / 4;
	int middle = l + half;
	median(a, l, l + quarter, middle, middle + quarter, h);
	TYPE pivot = a[middle];
	for (int i = l; i <= h;)
		if (a[i] < pivot)
			swap(a, i++, l++);
		else if (a[i] > pivot)
			swap(a, i, h--);
		else
			++i;
}

template <typename TYPE>
static void insertion(TYPE *a, int l, int h)
{
	for (int i = l + 1, j; i <= h; ++i) {
		TYPE t = a[i];
		for (j = i; j > l && t < a[j-1]; --j)
			a[j] = a[j-1];
		a[j] = t;
	}
}

template <typename TYPE>
static void sort(TYPE *a, int l, int h)
{
	if (l < h) {
		if (h - l < 32) {
			insertion(a, l, h);
		} else {
			int lt = l, gt = h;
			partition(a, lt, gt);
			sort(a, l, lt - 1);
			sort(a, gt + 1, h);
		}
	}
}

template <typename TYPE>
static void select(TYPE *a, int l, int h, int k)
{
	while (l < h) {
		if (h - l < 32) {
			insertion(a, l, h);
			break;
		} else {
			int lt = l, gt = h;
			partition(a, lt, gt);
			if (k < lt)
				h = lt - 1;
			else if (k > gt)
				l = gt + 1;
			else
				break;
		}
	}
}

}

template <typename TYPE>
void quick_sort(TYPE *a, int n)
{
	QUICK::sort(a, 0, n-1);
}

template <typename TYPE>
TYPE quick_select(TYPE *a, int k, int n)
{
	assert(n && k < n);
	QUICK::select(a, 0, n-1, k);
	return a[k];
}

}

