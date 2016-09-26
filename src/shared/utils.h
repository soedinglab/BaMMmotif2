#ifndef UTILS_H_
#define UTILS_H_

#include <algorithm>	// e.g. std::sort
#include <limits>		// e.g. std::numeric_limits
#include <numeric>		// e.g. std::numeric
#include <vector>

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
#endif‚

static char*							baseName( const char* filePath );
// calculate posterior probabilities from log likelihoods
std::vector<double>						calculatePosteriorProbabilities( std::vector<double> lLikelihoods );
static void								createDirectory( char* dir );
static std::vector< std::vector<int> >	generateFoldIndices( unsigned int N, unsigned int folds );
// calculate the power for integer base
static int								ipow( unsigned int base, int exp );
// sort in descending order using
static void								quickSort( std::vector<float> arr, int left, int right );

template <typename T>
std::vector<size_t> sortIndices( const std::vector<T> &v ); // returns a permutation which rearranges v into ascending order

inline char* baseName( const char* filePath ){

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

	char* baseName = ( char* )malloc( ( end-start+2 ) * sizeof( char ) );
	for( i = start; i <= end; i++ ){
		baseName[i-start] = filePath[i];
	}
	baseName[i-start] = '\0';

	return baseName;
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

		char* cmd = ( char* )calloc( 1024, sizeof( char ) );
		if( system( cmd ) != 0 ){
			fprintf( stderr, "Error: Directory %s could not be created\n", dir );
			exit( -1 );
		}
		free( cmd );
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

inline void quickSort( std::vector<float> arr, int left, int right ){

	int i = left, j = right;
	float tmp;
	float pivot = arr[( left + right ) / 2];

	/* partition */
	while( i <= j ){
		while( arr[i] - pivot > 0 )	i++;
		while( arr[j] - pivot < 0 )	j--;
		if( i <= j ){
			tmp = arr[i];
			arr[i] = arr[j];
			arr[j] = tmp;
			i++;
			j--;
		}
	}

	/* recursion */
	if( left < j )	quickSort( arr, left, j );
	if( i < right )	quickSort( arr, i, right );
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

inline double sign( double a, double b){
    // return value of a with sign of b
    if( ( a < 0.0 && b < 0.0 ) || (a > 0.0 && b > 0.0 ) ){
        return a;
    }else{
        return -a;
    }
}

template <class T> inline double zbrent(T &func, const double x1, const double x2, const double tol){
        // Using Brents method, return the root of a function of functior func known to lie between x1 and x2. The root will be refined until its accuracy it tol.

        // Maximum allowed number of iterations.
        const int ITMAX=100;
        //Machine floating-point precision.
        const double EPS = std::numeric_limits<double>::epsilon();
        double a=x1 , b=x2, c=x2, d, e, fa=func( a ), fb=func( b ), fc, p, q, r, s, tol1, xm;
        if(( fa > 0.0 && fb > 0.0 ) || ( fa < 0.0 && fb < 0.0 )){
            printf("\n Root must be bracketed in zbrent ");
            printf("\n fa = %f ; fb = %f \n ", fa,fb);
//            // find out which border ( min = x1 or max = x2 ) results in a higher likelihood.
            // do this in EM an here just return a flag
            return -1;
        }
        if( fa < fb ){
            printf("\n Root defines a minimum -> do not optimize.\n ");
            // root defines a minimum -> do not optimize
            return 0;
        }
        fc = fb;
        for( int iter = 0; iter < ITMAX; iter++){
            if(( fb > 0.0 && fc > 0.0 ) || ( fb < 0.0 && fc < 0.0 )){
                // Rename a, b, c and adjust bounding interval
                c = a;
                fc = fa;
                e = d = b-a;
            }
            if( fabs(fc) < fabs(fb) ){
                a = b;
                b = c;
                c = a;
                fa = fb;
                fb = fc;
                fc = fa;
            }
            // Convergence check.
            tol1 = 2.0* EPS * fabs(b) + 0.5 * tol;
            xm = 0.5 * ( c-b );
            if( fabs(xm) <= tol1 || fb == 0.0 ){
                return b;
            }
            if( fabs(e) >= tol1 && fabs(fa) > fabs(fb) ){
                // Attempt inverse quadratic interpolation
                s = fb / fa;
                if( a == c ){
                    p = 2.0 * xm * s;
                    q = 1.0 -s;
                } else {
                    q = fa / fc;
                    r = fb / fc;
                    p = s * (2.0 * xm * q * ( q - r ) - ( b - a ) * ( r - 1.0 ));
                    q = ( q - 1.0 ) * ( r - 1.0 ) * ( s * 1.0 );
                }
                if( p > 0.0 ){
                    q = -q;
                }
                p = fabs(p);
                double min1 = 3.0 * xm * q - fabs(tol1 * q);
                double min2 = fabs(e * q);
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
            if( fabs(d) > tol1 ){
                // Evaluate new trial root
                b += d;
                fb = func( b ); // THIS WAS ADDED BY ME!
            } else{
                b += sign( tol1, xm );
                fb = func( b );
            }
        }
        printf( "Maximum number of iterations exceeded in zbrent ");
        exit(0);
}


//
/** The digamma function in long double precision.
* @param x the real value of the argument
* @return the value of the digamma (psi) function at that point
* @author Richard J. Mathar
* @since 2005-11-24
*
* this piece of code is needed for the gradient calculation
* within the conjugate gradient for the alpha learning
*/
inline long double digamma(long double x)
{
	/* force into the interval 1..3 */
	if( x < 0.0L )
		return digamma(1.0L-x)+M_PIl/tanl(M_PIl*(1.0L-x)) ;	/* reflection formula */
	else if( x < 1.0L )
		return digamma(1.0L+x)-1.0L/x ;
	else if ( x == 1.0L)
		return -M_GAMMAl ;
	else if ( x == 2.0L)
		return 1.0L-M_GAMMAl ;
	else if ( x == 3.0L)
		return 1.5L-M_GAMMAl ;
	else if ( x > 3.0L)
		/* duplication formula */
		return 0.5L*(digamma(x/2.0L)+digamma((x+1.0L)/2.0L))+M_LN2l ;
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
		static long double Kncoe[] = { .30459198558715155634315638246624251L,
		.72037977439182833573548891941219706L, -.12454959243861367729528855995001087L,
		.27769457331927827002810119567456810e-1L, -.67762371439822456447373550186163070e-2L,
		.17238755142247705209823876688592170e-2L, -.44817699064252933515310345718960928e-3L,
		.11793660000155572716272710617753373e-3L, -.31253894280980134452125172274246963e-4L,
		.83173997012173283398932708991137488e-5L, -.22191427643780045431149221890172210e-5L,
		.59302266729329346291029599913617915e-6L, -.15863051191470655433559920279603632e-6L,
		.42459203983193603241777510648681429e-7L, -.11369129616951114238848106591780146e-7L,
		.304502217295931698401459168423403510e-8L, -.81568455080753152802915013641723686e-9L,
		.21852324749975455125936715817306383e-9L, -.58546491441689515680751900276454407e-10L,
		.15686348450871204869813586459513648e-10L, -.42029496273143231373796179302482033e-11L,
		.11261435719264907097227520956710754e-11L, -.30174353636860279765375177200637590e-12L,
		.80850955256389526647406571868193768e-13L, -.21663779809421233144009565199997351e-13L,
		.58047634271339391495076374966835526e-14L, -.15553767189204733561108869588173845e-14L,
		.41676108598040807753707828039353330e-15L, -.11167065064221317094734023242188463e-15L } ;

		register long double Tn_1 = 1.0L ;	/* T_{n-1}(x), started at n=1 */
		register long double Tn = x-2.0L ;	/* T_{n}(x) , started at n=1 */
		register long double resul = Kncoe[0] + Kncoe[1]*Tn ;

		x -= 2.0L ;

		for(unsigned int n = 2 ; n < sizeof(Kncoe)/sizeof(long double) ;n++)
		{
			const long double Tn1 = 2.0L * x * Tn - Tn_1 ;	/* Chebyshev recursion, Eq. 22.7.4 Abramowitz-Stegun */
			resul += Kncoe[n]*Tn1 ;
			Tn_1 = Tn ;
			Tn = Tn1 ;
		}
		return resul ;
	}
}

#endif /* UTILS_H_ */
