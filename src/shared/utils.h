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
#define M_PIf 3.141592653589f
#endif

#ifndef M_PId
/** The constant Pi for double */
#define M_PId 3.14159265358979363
#endif

#ifndef M_GAMMAf
/** Euler's constant for float */
#define M_GAMMAf 0.577215664901f
#endif

#ifndef M_GAMMAd
/** Euler's constant for double */
#define M_GAMMAd 0.5772156649015328
#endif

#ifndef M_LN2f
/** the natural logarithm of 2 for float */
#define M_LN2f 0.69314718f
#endif

#ifndef M_LN2d
/** the natural logarithm of 2 for double */
#define M_LN2d 0.693147180559945309
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
				.72037977f, -.12454959f,
				.27769457e-1f, -.67762371e-2f,
				.172387552e-2f, -.4481770e-3f,
				.11793660e-3f, -.312538943e-4f,
				.83173997e-5f, -.221914276e-5f,
				.59302267e-6f, -.158630512e-6f,
				.42459204e-7f, -.113691296e-7f,
				.30450222e-8f, -.815684551e-9f,
				.21852325e-9f, -.585464914e-10f,
				.15686348e-10f, -.420294963e-11f,
				.11261436e-11f, -.301743536e-12f,
				.80850955e-13f, -.216637798e-13f,
				.58047634e-14f, -.155537672e-14f,
				.41676109e-15f, -.111670651e-15f } ;


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

inline double digamma(double x)
{
	/* force into the interval 1..3 */
	if( x < 0.0 )
		return digamma(1.0-x)+M_PId/tan(M_PId*(1.0-x)) ;	/* reflection formula */
	else if( x < 1.0 )
		return digamma(1.0+x)-1.0/x ;
	else if ( x == 1.0)
		return -M_GAMMAd ;
	else if ( x == 2.0)
		return 1.0-M_GAMMAd ;
	else if ( x == 3.0)
		return 1.5-M_GAMMAd ;
	else if ( x > 3.0)
		/* duplication formula */
		return 0.5*(digamma(x/2.0)+digamma((x+1.0)/2.0))+M_LN2d ;
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
		static double Kncoe[] = {
				.30459198558715155634315638246624251,
				.72037977439182833573548891941219706, -.12454959243861367729528855995001087,
				.27769457331927827002810119567456810e-1, -.67762371439822456447373550186163070e-2,
				.17238755142247705209823876688592170e-2, -.44817699064252933515310345718960928e-3,
				.11793660000155572716272710617753373e-3, -.31253894280980134452125172274246963e-4,
				.83173997012173283398932708991137488e-5, -.22191427643780045431149221890172210e-5,
				.59302266729329346291029599913617915e-6, -.15863051191470655433559920279603632e-6,
				.42459203983193603241777510648681429e-7, -.11369129616951114238848106591780146e-7,
				.304502217295931698401459168423403510e-8, -.81568455080753152802915013641723686e-9,
				.21852324749975455125936715817306383e-9, -.58546491441689515680751900276454407e-10,
				.15686348450871204869813586459513648e-10, -.42029496273143231373796179302482033e-11,
				.11261435719264907097227520956710754e-11, -.30174353636860279765375177200637590e-12,
				.80850955256389526647406571868193768e-13, -.21663779809421233144009565199997351e-13,
				.58047634271339391495076374966835526e-14, -.15553767189204733561108869588173845e-14,
				.41676108598040807753707828039353330e-15, -.11167065064221317094734023242188463e-15 } ;

		register double Tn_1 = 1.0 ;	/* T_{n-1}(x), started at n=1 */
		register double Tn = x-2.0 ;	/* T_{n}(x) , started at n=1 */
		register double resul = Kncoe[0] + Kncoe[1]*Tn ;

		x -= 2.0 ;

		for(size_t n = 2 ; n < sizeof(Kncoe)/sizeof(double) ;n++)
		{
			const double Tn1 = 2.0 * x * Tn - Tn_1 ;	/* Chebyshev recursion, Eq. 22.7.4 Abramowitz-Stegun */
			resul += Kncoe[n]*Tn1 ;
			Tn_1 = Tn ;
			Tn = Tn1 ;
		}
		return resul ;
	}
}

#endif /* UTILS_H_ */
