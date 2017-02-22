/** \file util_math.cpp
 * \brief math utility functions
 *
 * Selection of math utility functions, mostly fast log/exp and logAdd (including some SSE2 versions). \n
 * Written by Mark Heron, based partially on util.C of Johannes Soeding.
 * \author Mark Heron
 * \author Johannes Soeding
 * \date 2012
 */

#include <cmath>
#include <cfloat>
#include <cstdio>
#include <cassert>

#include <vector>
#include <numeric>

#include <emmintrin.h>   // SSE2

#include "util_math.hpp"


template<> const float LOGSPACE<float>::ZERO = -252.0;
template<> const float LOGSPACE<float>::MIN  = -126.0;
template<> const float LOGSPACE<float>::ONE  =    0.0;
template<> const float LOGSPACE<float>::BOOL[2] = {-252.0, 0.0};

template<> const double LOGSPACE<double>::ZERO = -2044.0;
template<> const double LOGSPACE<double>::MIN  = -1022.0;
template<> const double LOGSPACE<double>::ONE  =     0.0;
template<> const double LOGSPACE<double>::BOOL[2] = {-2044.0, 0.0};

template<> const long double LOGSPACE<long double>::ZERO = -32764.0;
template<> const long double LOGSPACE<long double>::MIN  = -16382.0;
template<> const long double LOGSPACE<long double>::ONE  =      0.0;
template<> const long double LOGSPACE<long double>::BOOL[2] = {-32764.0, 0.0};



/**
 * \brief Sign function
 *
 * Return -1 for negative, 0 for 0 and 1 for positive numbers.
 * @param x
 * @return sign(x) (i.e. -1,0 or 1)
 */
template<typename numeric>
inline int sign(const numeric x) {
	return((x > 0) - (x < 0));
}
template int sign<float>(float x);
template int sign<double>(double x);
template int sign<long double>(long double x);



/**
 * \brief Power function for integers
 *
 * @param base
 * @param exp exponent
 * @return base^exponent
 */
int powInt(const int base, const int exp) {
	int result = 1;
	for(int i = 0; i < exp; i++) {
		result *= base;
	}
	return(result);
}



/**
 * \brief returns the sum of all values in the array
 *
 * Adds the first *length* values of the array.
 * Does not check if *length* is smaller than the actual array length!
 * @param array
 * @param length
 * @return sum of array elements
 */
template<typename numeric>
inline numeric sumArray(numeric const* const array, const int length) {
	return(std::accumulate(&array[0], &array[length], (numeric) 0.0));
}
template int sumArray<int>(int const* const x, const int length);
template float sumArray<float>(float const* const x, const int length);
template double sumArray<double>(double const* const x, int length);
template long double sumArray<long double>(long double const* const x, const int length);



/**
 * \brief returns the sum of all values in the vector
 *
 * @param vec the vector whose elements are summed
 * @return sum of vector elements
 */
template<typename numeric>
inline numeric sumVector(const std::vector<numeric> vec) {
	return(std::accumulate(vec.begin(), vec.end(), (numeric) 0.0));
}
template int sumVector<int>(const std::vector<int> vec);
template float sumVector<float>(const std::vector<float> vec);
template double sumVector<double>(const std::vector<double> vec);
template long double sumVector<long double>(const std::vector<long double> vec);


/**
 * \brief returns the sum of all values in all arrays
 *
 * Adds the first *cols* values of each of the first *rows* arrays in arraylist.
 * Does not check if *cols* or *rows* is smaller than the actual lengths!
 * @param arraylist
 * @param rows
 * @param cols
 * @return sum over the list of all array elements
 */
template<typename numeric>
numeric sumArrayList(numeric const*const*arraylist, const int rows, const int *cols) {
	numeric sum = 0.0;
	for(int i = 0; i < rows; i++) {
		sum += sumArray(arraylist[i], cols[i]);
	}
	return(sum);
}
template float sumArrayList<float>(float const*const*arraylist, const int rows, const int *cols);
template double sumArrayList<double>(double const*const*arraylist, const int rows, const int *cols);



/**
 * \brief returns the mean of all values in the array
 *
 * Adds the first *length* values of the array and divides by length.
 * Does not check if *length* is smaller than the actual array length!
 * @param array
 * @param length
 * @return mean of the array elements
 */
template<typename numeric>
inline numeric meanArray(numeric const* const array, const int length) {
	return(sumArray(array,length) / (numeric) length);
}
template float meanArray<float>(float const* const x, const int length);
template double meanArray<double>(double const* const x, const int length);
template long double meanArray<long double>(long double const* const x, const int length);


/**
 * \brief Divides an array by its sum so it sums to 1
 *
 * Normalizes the first *length* values of the array, doesn't check if *length* makes sense.
 *
 * @param array that is normalized
 * @param length of the array
 */
template<typename numeric>
void normalizeSumToOne(numeric *array, const int length) {
	numeric sum_corrector = (numeric) 1.0 / sumArray<numeric>(array, length);
	for(int i = 0; i < length; i++) {
		array[i] *= sum_corrector;
	}
}
template void normalizeSumToOne<float>(float *array, const int length);
template void normalizeSumToOne<double>(double *array, const int length);
template void normalizeSumToOne<long double>(long double *array, const int length);




/**
 * \brief log2(2^x+2^y)
 *
 * Adds two values in binary log space, i.e. log2(2^x+2^y) \n
 * where both x and y are already in log space, i.e. x=log2(a) and y=log2(b) and we really want log2(a+b) \n
 * It does not compute: log(x) + log(y) \n
 * IMPORTANT The parameters x and y have to be in log space. \n
 * Formula used: \n
 * log2(2^x+2^y) \n = log2(2^x * (1+2^(y-x)) \n = log2(2^x) + log2(1+2^(y-x)) \n = x + log2(1+2^(y-x))
 *
 * @param x in log space
 * @param y in log space
 * @return log2(2^x+2^y)
 */
template<typename numeric>
inline numeric logAdd(const numeric x, const numeric y) {
	assert(!(x != x) && "logAdd: x == NaN!");
	assert(!(y != y) && "logAdd: y == NaN!");

	if (x <= LOGSPACE<numeric>::MIN) {
		return(y);
	} else if (y <= LOGSPACE<numeric>::MIN) {
		return(x);
	} else if (x < y) {	// extract the larger
		return(y + log_fast<numeric>( (numeric) 1.0 + exp_fast<numeric>(x-y)));
	} else {
		return(x + log_fast<numeric>( (numeric) 1.0 + exp_fast<numeric>(y-x)));
	}
}
template float logAdd<float>(const float x, const float y);
template double logAdd<double>(const double x, const double y);
template long double logAdd<long double>(const long double x, const long double y);



/**
 * \brief log2(2^x-2^y)
 *
 * Subtracts two values in binary log space, i.e. log2(2^x-2^y) \n
 * where both x and y are already in log space, i.e. x=log2(a) and y=log2(b) and we really want log2(a-b) \n
 * It does not compute: log(x) - log(y) \n
 * IMPORTANT The parameters x and y have to be in log space. \n
 * Formula used: \n
 * log2(2^x-2^y) \n = log2(2^x * (1-2^(y-x)) \n = log2(2^x) + log2(1-2^(y-x)) \n = x + log2(1-2^(y-x))
 *
 * @param x in log space
 * @param y in log space
 * @return log2(2^x-2^y)
 */
template<typename numeric>
numeric logSub(const numeric x, const numeric y) {
	assert(!(x != x) && "logSub: x == NaN!");
	assert(!(y != y) && "logSub: y == NaN!");
	assert(x > y && "logSub: y > x -> log(<0)");

	if (y <= LOGSPACE<numeric>::MIN) {
		return(x);
	} else {
		return(x + log_fast<numeric>( (numeric) 1.0 - exp_fast<numeric>(y-x)));
	}
}
template float logSub<float>(const float x, const float y);
template double logSub<double>(const double x, const double y);
template long double logSub<long double>(const long double x, const long double y);



/* Code copied from Johannes Soeding as blueprint/inspiration for fast (SSE) exp2 and log2 code

// This function returns log2 with a max absolute deviation of +/- 1.5E-5 (typically 0.8E-5).
// It takes 0.80E-8 s  whereas log2(x) takes 5.4E-7 s. It is hence 9.4 times faster.
// It makes use of the representation of 4-byte floating point numbers:
// seee eeee emmm mmmm mmmm mmmm mmmm mmmm
// s is the sign,
// the following 8 bits, eee eee e, give the exponent + 127 (in hex: 0x7f).
// The following 23 bits give the mantisse, the binary digits after the decimal
// point:  x = (-1)^s * 1.mmmmmmmmmmmmmmmmmmmmmmm * 2^(eeeeeeee-127)
// In the code, *(int *)&x is an integer which contains the bytes as the
// floating point variable x is represented in memory. The expression
//     (((*(int *)&x) & 0x7f800000 ) >>23 )-0x7f is the exponent eeeeeeee,
// i.e., the largest integer that is smaller than log2(x) (e.g. -1 for 0.9).
float fast_log2(float x) {
	int i;
	static float lg2[1025]; // lg2[i] = log2[1+x/1024]
	static float diff[1025]; // diff[i]= (lg2[i+1]-lg2[i])/8096 (for interpolation)
	static char initialized = 0;
	if (x <= 0) {
		perror("took log of negative or 0 value");
		return -100000;
	}
	if (!initialized) //First fill in the arrays lg2[i] and diff[i]
	{
		float prev = 0.0f;
		lg2[0] = 0.0f;
		for (i = 1; i <= 1024; ++i) {
			lg2[i] = log( (float) (1024+i))*1.442695041-10.0f;
			diff[i - 1] = (lg2[i] - prev) * 1.2352E-4;
			prev = lg2[i];
		}
		initialized = 1;
	}
	int a = (((*((int *) &x)) & 0x7F800000) >> 23) - 0x7f; // exponent
	int b = ((*((int *) &x)) & 0x007FE000) >> 13; // first 10 bits of mantisse
	int c = ((*((int *) &x)) & 0x00001FFF); // further 13 bits of mantisse
	return a + lg2[b] + diff[b] * (float) (c);
}

/////////////////////////////////////////////////////////////////////////////////////
// fast 2^x
// ATTENTION: need to compile with g++ -fno-strict-aliasing when using -O2 or -O3!!!
// Relative deviation < 4.6E-6  (< 2.3E-7 with 5'th order polynomial)
// Speed: 2.1E-8s (2.3E-8s) per call! (exp(): 8.5E-8, pow(): 1.7E-7)
// Internal representation of float number according to IEEE 754:
//   1bit sign, 8 bits exponent, 23 bits mantissa: seee eeee emmm mmmm mmmm mmmm mmmm mmmm
//                                    0x4b400000 = 0100 1011 0100 0000 0000 0000 0000 0000
//   In summary: x = (-1)^s * 1.mmmmmmmmmmmmmmmmmmmmmm * 2^(eeeeeee-127)
/////////////////////////////////////////////////////////////////////////////////////
float fpow2(float x) {
	if (x > FLT_MAX_EXP)
		return FLT_MAX;
	if (x < FLT_MIN_EXP)
		return 0.0f;
	int *px = (int*) (&x); // store address of float as pointer to long int
	float tx = (x - 0.5f) + (3 << 22); // temporary value for truncation: x-0.5 is added to a large integer (3<<22),
									   // 3<<22 = (1.1bin)*2^23 = (1.1bin)*2^(150-127),
									   // which, in internal bits, is written 0x4b400000 (since 10010110bin = 150)
	int lx = *((int*) &tx) - 0x4b400000; // integer value of x
	float dx = x - (float) (lx); // float remainder of x
//   x = 1.0f + dx*(0.69606564f           // cubic approximation of 2^x for x in the range [0, 1]
//            + dx*(0.22449433f           // Gives relative deviation < 1.5E-4
//            + dx*(0.07944023f)));       // Speed: 1.9E-8s
	x = 1.0f + dx * (0.693019f // polynomial approximation of 2^x for x in the range [0, 1]
	+ dx * (0.241404f // Gives relative deviation < 4.6E-6
	+ dx * (0.0520749f // Speed: 2.1E-8s
	+ dx * 0.0134929f)));
//   x = 1.0f + dx*(0.693153f             // polynomial apporoximation of 2^x for x in the range [0, 1]
//            + dx*(0.240153f             // Gives relative deviation < 2.3E-7
//            + dx*(0.0558282f            // Speed: 2.3E-8s
//            + dx*(0.00898898f
//            + dx* 0.00187682f ))));
	*px += (lx << 23); // add integer power of 2 to exponent
	return x;
}

*/



/**
 * \brief SSE2 exp2(float)
 *
 * Computes an exp2 approximation for 4 floats at once using SSE2. \n
 * Intrinsic adaptation of assembler code published at:	\n
 * http://devmaster.net/forums/topic/6679-approximate-math-library/page__p__39242#entry39242
 * @param X exponent vector
 * @return vector of 2^exponents approximations
 */
__m128 _mm_exp2_ps(const __m128 X) {

	const __m128 LIMIT_MAX 	= _mm_set1_ps(129.0f); // 129.00000e+0f
	const __m128 LIMIT_MIN 	= _mm_set1_ps(-126.99999f);// -126.99999e+0f
	const __m128 HALF 		= _mm_set1_ps(0.5f);
	const __m128i EXPONENTIAL_BASE = _mm_set1_epi32(0x0000007F);// 127

	const __m128 C5 = (__m128) _mm_set1_epi32(0x3AF61905);// 1.8775767e-3f
	const __m128 C4 = (__m128) _mm_set1_epi32(0x3C134806);// 8.9893397e-3f
	const __m128 C3 = (__m128) _mm_set1_epi32(0x3D64AA23);// 5.5826318e-2f
	const __m128 C2 = (__m128) _mm_set1_epi32(0x3E75EAD4);// 2.4015361e-1f
	const __m128 C1 = (__m128) _mm_set1_epi32(0x3F31727B);// 6.9315308e-1f
	const __m128 C0 = (__m128) _mm_set1_epi32(0x3F7FFFFF);// 9.9999994e-1f

	__m128 X0;
	__m128 X1;
	__m128i X2;

	X0 = X;
	X0 = _mm_min_ps(X0, LIMIT_MAX);
	X0 = _mm_max_ps(X0, LIMIT_MIN);

	X1 = X0;
	X1 = _mm_sub_ps(X1, HALF);
	X2 = _mm_cvtps_epi32(X1);
	X1 = _mm_cvtepi32_ps(X2);
	X2 = _mm_add_epi32(X2, EXPONENTIAL_BASE);
	X2 = _mm_slli_epi32(X2, 23);


	// compute polynomial approximation of exp2 mantissa
	X0 = _mm_sub_ps(X0, X1);

	X1 = C5;
	X1 = _mm_mul_ps(X1, X0);
	X1 = _mm_add_ps(X1, C4);

	X1 = _mm_mul_ps(X1, X0);
	X1 = _mm_add_ps(X1, C3);

	X1 = _mm_mul_ps(X1, X0);
	X1 = _mm_add_ps(X1, C2);

	X1 = _mm_mul_ps(X1, X0);
	X1 = _mm_add_ps(X1, C1);

	X1 = _mm_mul_ps(X1, X0);
	X1 = _mm_add_ps(X1, C0);

	// multiply mantissa and exponent together
	X1 = _mm_mul_ps(X1, (__m128) X2);

	return(X1);
}



// Fast SSE2 log2 for four floats
// copied from Johannes Soeding
// Calculate integer of log2 for four floats in parallel with SSE2
// Maximum deviation: +/- 2.1E-5
// Run time: ~5.6ns on Intel core2 2.13GHz.
// For a negative argument, nonsense is returned. Otherwise, when <1E-38, a value
// close to -126 is returned and when >1.7E38, +128 is returned.
// The function makes use of the representation of 4-byte floating point numbers:
// seee eeee emmm mmmm mmmm mmmm mmmm mmmm
// s is the sign, eee eee e gives the exponent + 127 (in hex: 0x7f).
// The following 23 bits give the mantisse, the binary digits after the decimal
// point:  x = (-1)^s * 1.mmmmmmmmmmmmmmmmmmmmmmm * 2^(eeeeeeee-127)
// Therefore,  log2(x) = eeeeeeee-127 + log2(1.mmmmmm...)
//                     = eeeeeeee-127 + log2(1+y),  where y = 0.mmmmmm...
//                     ~ eeeeeeee-127 + ((a*y+b)*y+c)*y
// The coefficients a, b  were determined by a least squares fit, and c=1-a-b to get 1 at y=1.
// Lower/higher order polynomials may be used for faster or more precise calculation:
// Order 1: log2(1+y) ~ y
// Order 2: log2(1+y) = (a*y + 1-a)*y, a=-0.3427
//  => max dev = +/- 8E-3, run time ~ 3.8ns
// Order 3: log2(1+y) = ((a*y+b)*y + 1-a-b)*y, a=0.1564, b=-0.5773
//  => max dev = +/- 1E-3, run time ~ 4.4ns
// Order 4: log2(1+y) = (((a*y+b)*y+c)*y + 1-a-b-c)*y, a=-0.0803 b=0.3170 c=-0.6748
//  => max dev = +/- 1.4E-4, run time ~ 5.0ns?
// Order 5: log2(1+y) = ((((a*y+b)*y+c)*y+d)*y + 1-a-b-c-d)*y, a=0.0440047 b=-0.1903190 c=0.4123442 d=-0.7077702
//  => max dev = +/- 2.1E-5, run time ~ 5.6ns?
/**
 * \brief SSE2 log2(float)
 *
 * Computes an log2 approximation for 4 floats at once using SSE2. \n
 * Copied (with minor changes) from Johannes Soeding, see comment in code.	\n
 *
 * @param X vector
 * @return vector of log2(X) approximations
 */
__m128 _mm_log2_ps(const __m128 X) {

	const __m128i EXPONENTIAL_BASE 	 = _mm_set1_epi32(0x0000007F);// 127
	const __m128i CONST32_0x7fffff	 = _mm_set1_epi32(0x7fffff);
	const __m128i CONST32_0x3f800000 = _mm_set1_epi32(0x3f800000);
	const __m128  ONE = _mm_set1_ps(1.0);
	// const float a=0.1564, b=-0.5773, c=1.0-a-b;  // third order
	const float a=0.0440047f, b=-0.1903190f, c=0.4123442f, d=-0.7077702f, e=1.0f-a-b-c-d; // fifth order

	const __m128  C_A = _mm_set1_ps(a);
	const __m128  C_B = _mm_set1_ps(b);
	const __m128  C_C = _mm_set1_ps(c);
	const __m128  C_D = _mm_set1_ps(d);
	const __m128  C_E = _mm_set1_ps(e);

	__m128i E; // exponents of X
	__m128 Q; // intermediary result
	__m128 R; // result

	E = _mm_srli_epi32((__m128i) X, 23);    // shift right by 23 bits to obtain exponent+127
	E = _mm_sub_epi32(E, EXPONENTIAL_BASE);     // subtract 127 = 0x7f
	Q = (__m128) _mm_and_si128((__m128i) X, CONST32_0x7fffff);  // mask out exponent => mantisse
	Q = (__m128) _mm_or_si128((__m128i) Q, CONST32_0x3f800000); // set exponent to 127 (i.e., 0)
	Q = _mm_sub_ps(Q, ONE);          // subtract one from mantisse
	R = _mm_mul_ps(Q, C_A);           // R = a*X
	R = _mm_add_ps(R, C_B);           // R = a*X+b
	R = _mm_mul_ps(R, Q);                   // R = (a*X+b)*X
	R = _mm_add_ps(R, C_C);           // R = (a*X+b)*X+c
	R = _mm_mul_ps(R, Q);                   // R = ((a*X+b)*X+c)*X
	R = _mm_add_ps(R, C_D);           // R = ((a*X+b)*X+c)*X+d
	R = _mm_mul_ps(R, Q);                   // R = (((a*X+b)*X+c)*X+d)*X
	R = _mm_add_ps(R, C_E);           // R = (((a*X+b)*X+c)*X+d)*X+e
	R = _mm_mul_ps(R, Q);                   // R = ((((a*X+b)*X+c)*X+d)*X+e)*X ~ log2(1+X) !!

	R = _mm_add_ps(R, _mm_cvtepi32_ps(E));  // convert integer exponent to float and add to mantisse
	return(R);
}



/**
 * \brief SSE2 log2(2^X+2^Y)
 *
 * Approximation of addition of two numbers in log space (see formula below). \n
 * X and Y have to already be in log space.\n
 * \n
 * log2(2^X+2^Y) \n = log2(2^X * (1+2^(Y-X)) \n = log2(2^X) + log2(1+2^(Y-X)) \n = X + log2(1+2^(Y-X))
 *
 * @param X in log space
 * @param Y in log space
 * @return approximation of log2(2^X+2^Y)
 */
__m128 _mm_log2Add_ps(const __m128 X, const __m128 Y ) {

	const __m128  CONST32_1f = _mm_set1_ps(1.0);
	__m128 R;

	__m128 MAX = _mm_max_ps(X, Y);		// get the max which is extracted from the log
	__m128 MIN = _mm_min_ps(X, Y);		// get the min

	// computes:	MAX + log2( 1 + exp( MIN - MAX ))
	R = _mm_add_ps( MAX, _mm_log2_ps( _mm_add_ps( CONST32_1f, _mm_exp2_ps( _mm_sub_ps(MIN, MAX) ))));

	return(R);
}



/**
 * \brief SSE2 log2(2^X+2^Y)
 *
 * Approximation of addition of two numbers in log space. \n
 * Compared to _mm_logAdd_ps this function does not call _mm_exp2_ps and _mm_log2_ps, but has them implicitly. \n
 * However this doesn't make it faster so just use _mm_logAdd_ps. \n
 * X and Y have to already be in log space.
 *
 * @param X in log space
 * @param Y in log space
 * @return approximation of log2(2^X+2^Y)
 */
__m128 _mm_logAdd2_ps( const __m128 X, const __m128 Y ) {	// TODO ideally X wouldn't be modified but that

	const __m128  CONST32_1f = _mm_set1_ps(1.0);
	const __m128 limit_min = (__m128) _mm_set1_epi32(0xC2FDFFFF);// -126.99999e+0f
	const __m128 half = (__m128) _mm_set1_epi32(0x3F000000);// 0.5e+0f
	const __m128i exponent_base = _mm_set1_epi32(0x0000007F);// 127

	const __m128 C5 = (__m128) _mm_set1_epi32(0x3AF61905);// 1.8775767e-3f
	const __m128 C4 = (__m128) _mm_set1_epi32(0x3C134806);// 8.9893397e-3f
	const __m128 C3 = (__m128) _mm_set1_epi32(0x3D64AA23);// 5.5826318e-2f
	const __m128 C2 = (__m128) _mm_set1_epi32(0x3E75EAD4);// 2.4015361e-1f
	const __m128 C1 = (__m128) _mm_set1_epi32(0x3F31727B);// 6.9315308e-1f
	const __m128 C0 = (__m128) _mm_set1_epi32(0x3F7FFFFF);// 9.9999994e-1f

	const __m128i CONST32_0x7fffff = _mm_set_epi32(0x7fffff,0x7fffff,0x7fffff,0x7fffff);
	const __m128i CONST32_0x3f800000 = _mm_set_epi32(0x3f800000,0x3f800000,0x3f800000,0x3f800000);
	// const float a=0.1564, b=-0.5773, c=1.0-a-b;  // third order
	const float a=0.0440047f, b=-0.1903190f, c=0.4123442f, d=-0.7077702f, e=1.0f-a-b-c-d; // fifth order

	const __m128  CONST32_A = _mm_set1_ps(a);	// _mm_set_ps(a,a,a,a);
	const __m128  CONST32_B = _mm_set1_ps(b);
	const __m128  CONST32_C = _mm_set1_ps(c);
	const __m128  CONST32_D = _mm_set1_ps(d);
	const __m128  CONST32_E = _mm_set1_ps(e);
//	__m128i E; // exponents of X
	__m128 Q; //  intermediary result
	__m128 R; //  result

	__m128 MAX = _mm_max_ps(X, Y);		// get the max which is extracted from the log
	__m128 X0 = _mm_min_ps(X, Y);		// get the min
	__m128 X1;
	__m128i X2;

	//TODO add limits/check X and Y if they are basically log2(0)

	X0 = _mm_sub_ps(X0, MAX);
//	// X0 is negativ
	X0 = _mm_max_ps(X0, limit_min);

	X1 = X0;
	X1 = _mm_sub_ps(X1, half);
	X2 = _mm_cvtps_epi32(X1);
	X1 = _mm_cvtepi32_ps(X2);
	X2 = _mm_add_epi32(X2, exponent_base);
	X2 = _mm_slli_epi32(X2, 23);

	// compute polynomial approximation of exp2 mantissa
	X0 = _mm_sub_ps(X0, X1);

	X1 = C5;
	X1 = _mm_mul_ps(X1, X0);
	X1 = _mm_add_ps(X1, C4);

	X1 = _mm_mul_ps(X1, X0);
	X1 = _mm_add_ps(X1, C3);

	X1 = _mm_mul_ps(X1, X0);
	X1 = _mm_add_ps(X1, C2);

	X1 = _mm_mul_ps(X1, X0);
	X1 = _mm_add_ps(X1, C1);

	X1 = _mm_mul_ps(X1, X0);
	X1 = _mm_add_ps(X1, C0);

	// exp(MIN - MAX)
	Q = _mm_mul_ps(X1, (__m128) X2);

	Q = _mm_add_ps( Q, CONST32_1f);
	// X between 2 and 1

//	E = _mm_srli_epi32((__m128i) X, 23);    // shift right by 23 bits to obtain exponent+127
//	E = _mm_sub_epi32(E, CONST32_0x7f);     // subtract 127 = 0x7f
	// E = 0 because X > 0 and X < 2

	Q = (__m128) _mm_and_si128((__m128i) Q, CONST32_0x7fffff);  // mask out exponent => mantisse
	Q = (__m128) _mm_or_si128((__m128i) Q, CONST32_0x3f800000); // set exponent to 127 (i.e., 0)

	Q = _mm_sub_ps(Q, CONST32_1f);          // subtract one from mantisse
	R = _mm_mul_ps(Q, CONST32_A);           // R = a*X
	R = _mm_add_ps(R, CONST32_B);           // R = a*X+b
	R = _mm_mul_ps(R, Q);                   // R = (a*X+b)*X
	R = _mm_add_ps(R, CONST32_C);           // R = (a*X+b)*X+c
	R = _mm_mul_ps(R, Q);                   // R = ((a*X+b)*X+c)*X
	R = _mm_add_ps(R, CONST32_D);           // R = ((a*X+b)*X+c)*X+d
	R = _mm_mul_ps(R, Q);                   // R = (((a*X+b)*X+c)*X+d)*X
	R = _mm_add_ps(R, CONST32_E);           // R = (((a*X+b)*X+c)*X+d)*X+e
	R = _mm_mul_ps(R, Q);                   // R = ((((a*X+b)*X+c)*X+d)*X+e)*X ~ log2(1+X) !!

	// E is still 0
//	R = _mm_add_ps(R, _mm_cvtepi32_ps(E));  // convert integer exponent to float and add to mantisse

	R = _mm_add_ps( R, MAX);

	return(R);
}
