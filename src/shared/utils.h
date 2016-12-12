#ifndef UTILS_H_
#define UTILS_H_

#include <algorithm>	// e.g. std::sort
#include <limits>		// e.g. std::numeric_limits
#include <numeric>		// e.g. std::numeric
#include <vector>
#include <memory>

#include <sys/stat.h>	// e.g. stat

#ifndef M_PIl
/** The constant Pi in high precision */
#define M_PIl 3.1415926535897932384626433832795029L
#endif

#ifndef M_GAMMAl
/** Euler's constant in high precision */
#define M_GAMMAl 0.5772156649015328606065120900824024L
#endif

#ifndef M_LN2l
/** the natural logarithm of 2 in high precision */
#define M_LN2l 0.6931471805599453094172321214581766L
#endif

#ifndef M_PIf
/** The constant Pi for float */
#define M_PIf 3.14159265f
#endif

#ifndef M_GAMMAf
/** Euler's constant for float */
#define M_GAMMAf 0.57721566f
#endif

#ifndef M_LN2f
/** the natural logarithm of 2 for float */
#define M_LN2f 0.69314718f
#endif

#ifndef CALL_EM_FN
#define CALL_EM_FN( object, ptrToFunc )( ( object ).*( ptrToFunc ) )
#endif

#ifndef make_unique
// note: this implementation does not disable this overload for array types
template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args)
{
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}
#endif

static std::string						baseName( const char* filePath );
// calculate posterior probabilities from log likelihoods
std::vector<double>						calculatePosteriorProbabilities( std::vector<double> lLikelihoods );
static void								createDirectory( char* dir );
static std::vector< std::vector<int> >	generateFoldIndices( unsigned int N, unsigned int folds );
// calculate the power for integer base
static int								ipow( unsigned int base, int exp );

template <typename T>
std::vector<size_t> sortIndices( const std::vector<T> &v ); // returns a permutation which rearranges v into ascending order

inline std::string baseName( const char* filePath ){

	int i = 0, start = 0, end = 0;

	while( filePath[++i] != '\0' ){
		if( filePath[i] == '.' ){
			end = i - 1;
		}
	}
	while( --i != 0 && filePath[i] != '/' ){
		;
	}
	if( i != 0 ){
		start = i + 1;
	}

	std::string basename( filePath, start, end-start+1 );

	return basename;
}

inline std::vector<double> calculatePosteriorProbabilities( std::vector<double> lLikelihoods ){

	// see http://stats.stackexchange.com/questions/66616/converting-normalizing-very-small-likelihood-values-to-probability/66621#66621

	int d = std::numeric_limits<double>::digits10; // digits of precision

	double epsilon = pow( 10, -d );
	long unsigned int N = lLikelihoods.size();

	double limit = log( epsilon ) - log( static_cast<double>( N ) );

	// sort indices into ascending order
	std::vector<size_t> order = sortIndices( lLikelihoods );
	// sort likelihoods into ascending order
	std::sort( lLikelihoods.begin(), lLikelihoods.end() );

	for( unsigned int i = 0; i < N; i++ ){
		lLikelihoods[i] = lLikelihoods[i] - lLikelihoods[N-1];
	}

	double partition = 0.0;

	for( unsigned int i = 0; i < N; i++ ){
		if( lLikelihoods[i] >= limit ){
			lLikelihoods[i] = exp( lLikelihoods[i] );
			partition += lLikelihoods[i];
		}
	}

	std::vector<double> posteriors( N );

	for( unsigned int i = 0; i < N; i++ ){
		if( lLikelihoods[i] >= limit ){
			posteriors[order[i]] = lLikelihoods[i] / partition;
		} else{
			posteriors[order[i]] = 0.0;
		}
	}

	return( posteriors );
}

inline void createDirectory( char* dir ){

	struct stat fileStatus;

	if( stat( dir, &fileStatus ) != 0 ){
		std::cout << "Status: Output directory does not exist. "
				"New directory is created automatically.\n";
		if( system( ( "mkdir " + std::string( dir ) ).c_str() ) != 0 ){
			std::cerr << "Error: Directory " << dir << " could not be created" << std::endl;
			exit( -1 );
		}
	}
}

inline std::vector< std::vector<int> > generateFoldIndices( unsigned int N, unsigned int folds ){

	std::vector< std::vector<int> > indices( folds );

	for( unsigned int i = 0; i < N; i += folds ){
		for( unsigned int j = 0; j < folds; j++ ){
		    if( i+j < N ){
				indices[j].push_back( i+j );
		    }
		}
	}
	return indices;
}

inline int ipow( unsigned int base, int exp ){

    int res = 1;

    while( exp ){
        if( exp & 1 )
        	res *= base;
        exp >>= 1;
        base *= base;
    }

    return res;
}

template <typename T>
inline std::vector<size_t> sortIndices( const std::vector<T>& v ){

  // initialize with original indices
  std::vector<size_t> idx( v.size() );
  iota( idx.begin(), idx.end(), 0 );

  // sort indices based on comparing values in v
  sort( idx.begin(), idx.end(),
		[&v]( size_t i1, size_t i2 ){
	        return v[i1] < v[i2];
        }
      );

  return idx;
}

inline float sign( float a, float b ){
    // return value of a with sign of b
    if( ( a < 0.0 && b < 0.0 ) || (a > 0.0 && b > 0.0 ) ){
        return a;
    } else {
        return -a;
    }
}


template <class T> inline float zbrent(T &obj, float ( T::*func )( float , int ), const float x1, const float x2, const float tol, int k){
        // Using Brents method, return the root of a function of function func known to lie between x1 and x2. The root will be refined until its accuracy it tol.

        // Maximum allowed number of iterations.
        const int ITMAX=std::numeric_limits<int>::max();
        //Machine floating-point precision.
        const float EPS = std::numeric_limits<float>::epsilon();
        float a=x1 , b=x2, c=x2, d=b-a, e=b-a, fa=CALL_EM_FN( obj, func )( a ,k ), fb=CALL_EM_FN( obj, func )( b ,k ), fc, p, q, r, s, tol1, xm;

        if( fa < fb ){
            printf("\n Root defines a minimum -> do not optimize.\n ");
            // root defines a minimum -> do not optimize
            return -1;
        }
        if(( fa > 0.0 && fb > 0.0 ) || ( fa < 0.0 && fb < 0.0 )){
            printf("\n Root must be bracketed for order %d ", k);
            printf("\n fa = %f ; fb = %f ", fa,fb);
            // find out which border ( min = x1 or max = x2 ) results in a higher likelihood.
            // do this in EM an here just return a flag
//            if( fabsf(fa) < fabsf(fb)){
//                printf(" --> choose %f ", a);
//                return a;
//            }else{
//                printf(" --> choose %f ", b);
//                return b;
//            }
            return -1;
        }

        fc = fb;
        for( int iter = 0; iter < ITMAX; iter++){
            if(( fb > 0.0 && fc > 0.0 ) || ( fb < 0.0 && fc < 0.0 )){
                // Rename a, b, c and adjust bounding interval
                c = a;
                fc = fa;
                e = d = b-a;
            }
            if( fabsf(fc) < fabsf(fb) ){
                a = b;
                b = c;
                c = a;
                fa = fb;
                fb = fc;
                fc = fa;
            }
            // Convergence check.
            tol1 = 2.0f * EPS * fabsf(b) + 0.5f * tol;
            xm = 0.5f * ( c-b );
            if( fabsf(xm) <= tol1 || fb == 0.0 ){
                return b;
            }
            if( fabsf(e) >= tol1 && fabsf(fa) > fabsf(fb) ){
                // Attempt inverse quadratic interpolation
                s = fb / fa;
                if( a == c ){
                    p = 2.0f * xm * s;
                    q = 1.0f -s;
                } else {
                    q = fa / fc;
                    r = fb / fc;
                    p = s * (2.0f * xm * q * ( q - r ) - ( b - a ) * ( r - 1.0f ));
                    q = ( q - 1.0f ) * ( r - 1.0f ) * ( s * 1.0f );
                }
                if( p > 0.0 ){
                    q = -q;
                }
                p = fabsf(p);
                float min1 = 3.0f * xm * q - fabsf( tol1 * q );
                float min2 = fabsf(e * q);
                if(2.0 * p < ( min1 < min2 ? min1 : min2 )){
                    // Accept interpolation.
                    e = d;
                    d = p/q;
                } else{
                    // Interpolation failed, use bisection
                    d = xm;
                    e = d;
                }
            } else {
                // Bounds decreasing to slowly, use bisection
                d = xm;
                e = d;
            }
            // move last best guess to a
            a = b;
            fa = fb;
            if( fabsf(d) > tol1 ){
                // Evaluate new trial root
                b += d;
                fb = CALL_EM_FN( obj, func )( b , k ); // THIS WAS ADDED BY ME!
            } else{
                b += sign( tol1, xm );
                fb = CALL_EM_FN( obj, func )( b , k);
            }
        }
        printf( "Maximum number of iterations exceeded in zbrent ");
        exit(0);
}

/** The digamma function in long double precision.
* @param x the real value of the argument
* @return the value of the digamma (psi) function at that point
* @author Richard J. Mathar
* @since 2005-11-24
*
* EDIT: precision down scaled to floats only
* -> this function is needed within the gradient calculation of the Q-function
*/

inline float digammaf( float x ){
	/* force into the interval 1..3 */
	if( x < 0.0f )
		return digammaf(1.0f-x)+M_PIf/tanf(M_PIf*(1.0f-x)) ;	/* reflection formula */
	else if( x < 1.0f )
		return digammaf(1.0f+x)-1.0f/x ;
	else if ( x == 1.0f)
		return -M_GAMMAf ;
	else if ( x == 2.0f)
		return 1.0f-M_GAMMAf ;
	else if ( x == 3.0f)
		return 1.5f-M_GAMMAf ;
	else if ( x > 3.0f)
		/* duplication formula */
		return 0.5f*(digammaf(x/2.0f)+digammaf((x+1.0f)/2.0f))+M_LN2f ;
	else
	{
		/* Just for your information, the following lines contain
		* the Maple source code to re-generate the table that is
		* eventually becoming the Kncoe[] array below
		* interface(prettyprint=0) :
		* Digits := 63 :
		* r := 0 :
		*
		* for l from 1 to 60 do
		* 	d := binomial(-1/2,l) :
		* 	r := r+d*(-1)^l*(Zeta(2*l+1) -1) ;
		* 	evalf(r) ;
		* 	print(%,evalf(1+Psi(1)-r)) ;
		*o d :
		*
		* for N from 1 to 28 do
		* 	r := 0 :
		* 	n := N-1 :
		*
 		*	for l from iquo(n+3,2) to 70 do
		*		d := 0 :
 		*		for s from 0 to n+1 do
 		*		 d := d+(-1)^s*binomial(n+1,s)*binomial((s-1)/2,l) :
 		*		od :
 		*		if 2*l-n > 1 then
 		*		r := r+d*(-1)^l*(Zeta(2*l-n) -1) :
 		*		fi :
 		*	od :
 		*	print(evalf((-1)^n*2*r)) ;
 		*od :
 		*quit :
		*/
		static float Kncoe[] = {
		.30459198f,
		.72037977f,		-.12454959f,
		.27769457e-1f, 	-.67762371e-2f,
		.17238755e-2f, 	-.44817699e-3f,
		.11793660e-3f, 	-.31253894e-4f,
		.83173997e-5f, 	-.22191427e-5f,
		.59302266e-6f, 	-.15863051e-6f,
		.42459203e-7f, 	-.11369129e-7f,
		.30450221e-8f, 	-.81568455e-9f,
		.21852324e-9f, 	-.58546491e-10f,
		.15686348e-10f, -.42029496e-11f,
		.11261435e-11f, -.30174353e-12f,
		.80850955e-13f, -.21663779e-13f,
		.58047634e-14f, -.15553767e-14f,
		.41676108e-15f,	-.11167065e-15f } ;

		register float Tn_1 = 1.0f ;	/* T_{n-1}(x), started at n=1 */
		register float Tn = x-2.0f ;	/* T_{n}(x) , started at n=1 */
		register float result = Kncoe[0] + Kncoe[1]*Tn ;

		x -= 2.0f;

		for( size_t n = 2; n < sizeof( Kncoe ) / sizeof( float ); n++ ){
			const float Tn1 = 2.0f * x * Tn - Tn_1;	/* Chebyshev recursion, Eq. 22.7.4 Abramowitz-Stegun */
			result += Kncoe[n] * Tn1;
			Tn_1 = Tn;
			Tn = Tn1;
		}
		return result;
	}
}
#endif /* UTILS_H_ */
