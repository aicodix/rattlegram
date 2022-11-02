/*
Generic pixel buffer container

Copyright 2021 Ahmet Inan <inan@aicodix.de>
*/

#pragma once

template<typename TYPE, int WIDTH, int HEIGHT>
struct Image {
	static const int width = WIDTH;
	static const int height = HEIGHT;
	static const int length = WIDTH * HEIGHT;
	TYPE *pixels;

	explicit Image(TYPE *pixels) : pixels(pixels) {}

	static TYPE overwrite(int, int, TYPE, TYPE c) {
		return c;
	}

	template<typename OP>
	void set(int x, int y, TYPE c, OP op) {
		if (0 <= x && x < width && 0 <= y && y < height)
			pixels[width * y + x] = op(x, y, pixels[width * y + x], c);
	}

	void set(int x, int y, TYPE c) {
		if (0 <= x && x < width && 0 <= y && y < height)
			pixels[width * y + x] = c;
	}

	template<typename OP>
	void fill(TYPE c, OP op) {
		for (int y = 0; y < height; ++y)
			for (int x = 0; x < width; ++x)
				pixels[width * y + x] = op(x, y, pixels[width * y + x], c);
	}

	void fill(TYPE c) {
		for (int i = 0; i < length; ++i)
			pixels[i] = c;
	}

	template<typename OP>
	void vline(int x, TYPE c, OP op) {
		if (0 <= x && x < width)
			for (int i = 0; i < height; ++i)
				set(x, i, c, op);
	}

	void vline(int x, TYPE c) {
		if (0 <= x && x < width)
			for (int i = 0; i < height; ++i)
				set(x, i, c);
	}

	template<typename OP>
	void hline(int y, TYPE c, OP op) {
		if (0 <= y && y < height)
			for (int i = 0; i < width; ++i)
				set(i, y, c, op);
	}

	void hline(int y, TYPE c) {
		if (0 <= y && y < height)
			for (int i = 0; i < width; ++i)
				set(i, y, c);
	}

	// only made for abs(x0-x1) <= 1
	template<typename OP>
	void vline(int x0, int y0, int x1, int y1, TYPE c, OP op) {
		int a0 = std::min(y0, (y0 + y1) / 2);
		int a1 = std::max(y0, (y0 + y1) / 2);
		for (int y = a0; y <= a1; ++y)
			set(x0, y, c, op);
		int b0 = std::min((y0 + y1) / 2, y1);
		int b1 = std::max((y0 + y1) / 2, y1);
		for (int y = b0; y <= b1; ++y)
			set(x1, y, c, op);
	}

	void vline(int x0, int y0, int x1, int y1, TYPE c) {
		vline(x0, y0, x1, y1, c, overwrite);
	}

	// only made for abs(y0-y1) <= 1
	template<typename OP>
	void hline(int x0, int y0, int x1, int y1, TYPE c, OP op) {
		int a0 = std::min(x0, (x0 + x1) / 2);
		int a1 = std::max(x0, (x0 + x1) / 2);
		for (int x = a0; x <= a1; ++x)
			set(x, y0, c, op);
		int b0 = std::min((x0 + x1) / 2, x1);
		int b1 = std::max((x0 + x1) / 2, x1);
		for (int x = b0; x <= b1; ++x)
			set(x, y1, c, op);
	}

	void hline(int x0, int y0, int x1, int y1, TYPE c) {
		hline(x0, y0, x1, y1, c, overwrite);
	}

	template<typename OP>
	void low_angle(int x0, int y0, int x1, int y1, TYPE c, OP op) {
		int dx = x1 - x0;
		int dy = y1 - y0;
		int yi = 1;
		if (dy < 0) {
			yi = -1;
			dy = -dy;
		}
		int D = 2 * dy - dx;
		int y = y0;
		for (int x = x0; x <= x1; ++x) {
			set(x, y, c, op);
			if (D > 0) {
				y += yi;
				D -= 2 * dx;
			}
			D += 2 * dy;
		}
	}

	template<typename OP>
	void high_angle(int x0, int y0, int x1, int y1, TYPE c, OP op) {
		int dx = x1 - x0;
		int dy = y1 - y0;
		int xi = 1;
		if (dx < 0) {
			xi = -1;
			dx = -dx;
		}
		int D = 2 * dx - dy;
		int x = x0;
		for (int y = y0; y <= y1; ++y) {
			set(x, y, c, op);
			if (D > 0) {
				x += xi;
				D -= 2 * dy;
			}
			D += 2 * dx;
		}
	}

	template<typename OP>
	void line(int x0, int y0, int x1, int y1, TYPE c, OP op) {
		if (abs(y0 - y1) < abs(x0 - x1)) {
			if (x0 < x1)
				low_angle(x0, y0, x1, y1, c, op);
			else
				low_angle(x1, y1, x0, y0, c, op);
		} else {
			if (y0 < y1)
				high_angle(x0, y0, x1, y1, c, op);
			else
				high_angle(x1, y1, x0, y0, c, op);
		}
	}

	void line(int x0, int y0, int x1, int y1, TYPE c) {
		line(x0, y0, x1, y1, c, overwrite);
	}
};
