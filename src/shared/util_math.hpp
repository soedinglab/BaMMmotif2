/** \file util_math.hpp
 * \brief header file for math functions
 *
 * \author Mark Heron
 * \date 2012
 */

#ifndef UTIL_MATH_H_
#define UTIL_MATH_H_

#include <cmath>
#include <cassert>

#include <vector>

#include <emmintrin.h>	// SSE2

// util_math.c
/** \defgroup u_math Utility Math Functions
 * Common math functions and wrappers that don't exist in C
 */

template<typename numeric>
struct LOGSPACE {
    static const numeric ZERO;
    static const numeric MIN;
    static const numeric ONE;
    static const numeric BOOL[2];
};

/** \defgroup u_math_simple Simple math functions
 * Simple functions that (for some unknown reason) are not implemented in math.h
 * \ingroup u_math
 * @{
 */
template<typename numeric>
int sign(const numeric x);
int powInt(const int base, const int exp);
/** @} */

/** \defgroup u_math_array Functions for arrays and vectors
 * \ingroup u_math
 * @{
 */
template<typename numeric>
numeric sumArrayList(numeric const*const*float_arraylist, const int rows, const int *cols);

template<typename numeric>
numeric sumArray(numeric const* const array, const int array_length);

template<typename numeric>
numeric sumVector(const std::vector<numeric> vec);

template<typename numeric>
numeric meanArray(numeric const* const array, const int array_length);


template<typename numeric>
void normalizeSumToOne(numeric *array, const int length);
/** @} */


/** \defgroup u_math_fast Fast function wrappers
 * wrappers to easily switch between implementations of  fast log2 and exp2 functions
 * \ingroup u_math
 * @{
 */


/**
 * \brief Fast logarithm wrapper (uses base 2 at the moment, but that shouldn't be important for your use case!)
 *
 * Wrapper to simply switch between precise and approximate log computations in code,
 * and has additional asserts.
 *
 * @param x
 * @return basically log2(x) (at the moment)
 */
template<typename numeric>
inline numeric log_fast(const numeric x) {
	assert(!(x != x) && "log_fast: x == NaN!");
	assert(x >= 0 && "log_fast: x < 0");
//	assert(x > LOGSPACE<numeric>::ZERO && "log_fast: x < log2(0)!");
	if(x == 0) {
		return( LOGSPACE<numeric>::ZERO);
	} else {
		return( (numeric) log( (double) x));
	}
}
template float log_fast<float>(const float x);
template double log_fast<double>(const double x);

template<>
inline long double log_fast<long double>(const long double x) {
	assert(!(x != x) && "log_fast: x == NaN!");
	assert(x >= 0 && "log_fast: x < 0");
	//	assert(x > LOGSPACE<long double>::ZERO && "log_fast: x < log2(0)!");
	if(x == 0) {
		return( LOGSPACE<double>::ZERO);
	} else {
		return( logl(x) );
	}
}


/**
 * \brief Fast exponentiation wrapper (uses base 2 at the moment, but that shouldn't be important for your use case!)
 *
 * Wrapper to simplify switching between precise and approximate exponentiation computation in code.
 *
 * @param x exponent
 * @return basically 2^x
 */
template<typename numeric>
inline numeric exp_fast(const numeric x) {
	return( (numeric) exp( (double) x));
}
template float exp_fast<float>(const float x);
template double exp_fast<double>(const double x);

template<>
inline long double exp_fast<long double>(const long double x) {
	return( expl(x) );
}
template long double exp_fast<long double>(const long double x);
/** @} */

/** \defgroup u_math_logspace Log-space functions
 * Addition and subtraction in log-space
 * \ingroup u_math
 * @{
 */
template<typename numeric>
numeric logAdd(const numeric x, const numeric y);

template<typename numeric>
numeric logSub(const numeric x, const numeric y);
/** @} */

/** \defgroup u_math_sse SSE functions
 * exp2, log2 and logAdd implemented using SSE2
 * \ingroup u_math
 * @{
 */
__m128 _mm_exp2_ps(__m128 X);
__m128 _mm_log2_ps(__m128 X);
__m128 _mm_log2Add_ps( __m128 X, __m128 Y);
/** @}/ */

#endif /* UTIL_MATH_H_ */
