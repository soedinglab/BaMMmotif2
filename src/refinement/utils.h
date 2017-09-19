#ifndef UTILS_H_
#define UTILS_H_

#include <algorithm>	// e.g. std::sort
#include <limits>		// e.g. std::numeric_limits
#include <numeric>		// e.g. std::numeric
#include <vector>
#include <memory>

#include <sys/stat.h>	// e.g. stat

#ifndef M_GAMMAl
/** Euler's constant in high precision */
#define M_GAMMAl 0.5772156649015328606065120900824024L
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

static std::string						baseName( const char* filePath );
// calculate posterior probabilities from log likelihoods
std::vector<double>						calculatePosteriorProbabilities( std::vector<double> lLikelihoods );
static void								createDirectory( char* dir );

// calculate the power for integer base
static size_t							ipow( size_t base, size_t exp );

// returns a permutation which rearranges v into ascending order
template <typename T> std::vector<size_t> sortIndices( const std::vector<T> &v );

inline std::string baseName( const char* filePath ){

	size_t i = 0, start = 0, end = 0;

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

	for( size_t i = 0; i < N; i++ ){
		if( lLikelihoods[i] >= limit ){
			lLikelihoods[i] = exp( lLikelihoods[i] );
			partition += lLikelihoods[i];
		}
	}

	std::vector<double> posteriors( N );

	for( size_t i = 0; i < N; i++ ){
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
		std::cout << "New output directory is created automatically.\n";
		if( system( ( "mkdir " + std::string( dir ) ).c_str() ) != 0 ){
			std::cerr << "Error: Directory " << dir << " could not be created." << std::endl;
			exit( -1 );
		}
	}
}

inline size_t ipow( size_t base, size_t exp ){

	size_t res = 1;

    while( exp ){
        if( exp & 1 )
        	res *= base;
        exp >>= 1;
        base *= base;
    }

    return res;
}

template <typename T> inline std::vector<size_t> sortIndices( const std::vector<T>& v ){

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


		float Tn_1 = 1.0f ;	/* T_{n-1}(x), started at n=1 */
		float Tn = x-2.0f ;	/* T_{n}(x) , started at n=1 */
		float result = Kncoe[0] + Kncoe[1]*Tn ;

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

		double Tn_1 = 1.0 ;	/* T_{n-1}(x), started at n=1 */
		double Tn = x-2.0 ;	/* T_{n}(x) , started at n=1 */
		double resul = Kncoe[0] + Kncoe[1]*Tn ;

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
