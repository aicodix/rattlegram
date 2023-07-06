#pragma once

template<int NUM, typename TYPE, typename CODE>
struct QAM;

template<typename TYPE, typename CODE>
struct QAM<16, TYPE, CODE> {
    static const int NUM = 16;
	static const int BITS = 4;
	typedef TYPE complex_type;
	typedef typename TYPE::value_type value_type;
	typedef CODE code_type;

    static constexpr value_type LEVEL =  0.31622776601683794;
    static constexpr value_type DIST = 2 * LEVEL;

    static constexpr complex_type qam16[16] = {
        complex_type(-3.0 * LEVEL, -3.0 * LEVEL),
        complex_type(3.0 * LEVEL, -3.0 * LEVEL),
        complex_type(-1.0 * LEVEL, -3.0 * LEVEL),
        complex_type(1.0 * LEVEL, -3.0 * LEVEL),
        complex_type(-3.0 * LEVEL, 3.0 * LEVEL),
        complex_type(3.0 * LEVEL, 3.0 * LEVEL),
        complex_type(-1.0 * LEVEL, 3.0 * LEVEL),
        complex_type(1.0 * LEVEL, 3.0 * LEVEL),
        complex_type(-3.0 * LEVEL, -1.0 * LEVEL),
        complex_type(3.0 * LEVEL, -1.0 * LEVEL),
        complex_type(-1.0 * LEVEL, -1.0 * LEVEL),
        complex_type(1.0 * LEVEL, -1.0 * LEVEL),
        complex_type(-3.0 * LEVEL, 1.0 * LEVEL),
        complex_type(3.0 * LEVEL, 1.0 * LEVEL),
        complex_type(-1.0 * LEVEL, 1.0 * LEVEL),
        complex_type(1.0 * LEVEL, 1.0 * LEVEL),
    };

	static code_type quantize(value_type precision, value_type value) {
		value *= DIST * precision;
		if (std::is_integral<code_type>::value)
			value = std::nearbyint(value);
		if (std::is_same<code_type, int8_t>::value)
			value = std::min<value_type>(std::max<value_type>(value, -128), 127);
		return value;
	}

    static void hard(code_type *b, complex_type c) {
		b[0] = (c.real() > value_type(0)) ? code_type(1) : code_type(-1);
        b[1] = (abs(c.real()) < DIST) ? code_type(1) : code_type(-1); 
        b[2] = (c.imag() > value_type(0)) ? code_type(1) : code_type(-1); 
        b[3] = (abs(c.imag()) < DIST) ? code_type(1) : code_type(-1);
	}

	static void soft(code_type *b, complex_type c, value_type precision) {
		b[0] = quantize(precision, c.real());
        b[1] = quantize(precision, (value_type(DIST) - abs(c.real())));
        b[2] = quantize(precision, c.imag());
        b[3] = quantize(precision, (value_type(DIST) - abs(c.imag())));
	}

	static complex_type map(code_type *b) {
        uint8_t index = 0;
        for (int i = 0; i < BITS; i++) {
            if (b[i] > code_type(0)) {
                index |= 1 << i;
            }
        }
        complex_type ret = qam16[index];
        return ret;
	}
};