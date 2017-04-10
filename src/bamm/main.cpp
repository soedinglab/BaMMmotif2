#include <time.h>		// time(), clock_t, clock, CLOCKS_PER_SEC
#include <stdio.h>

#include "Global.h"
#include "../shared/BackgroundModel.h"
#include "../shared/utils.h"
#include "MotifSet.h"
#include "ModelLearning.h"
#include "ScoreSeqSet.h"
#include "SeqGenerator.h"

#include "FDR.h"

int main( int nargs, char* args[] ){

	// seed random number
	srand( 42 );

	clock_t t0 = clock();

	fprintf( stderr, "\n" );
	fprintf( stderr, "======================================\n" );
	fprintf( stderr, "=      Welcome to use BaMMmotif      =\n" );
	fprintf( stderr, "=                   Version 1.0      =\n" );
	fprintf( stderr, "=              by Soeding Group      =\n" );
	fprintf( stderr, "=  http://www.mpibpc.mpg.de/soeding  =\n" );
	fprintf( stderr, "======================================\n" );

	// initialization
	Global::init( nargs, args );

	Global::rngx.seed( 42 );

	fprintf( stderr, "\n" );
	fprintf( stderr, "************************\n" );
	fprintf( stderr, "*   Background Model   *\n" );
	fprintf( stderr, "************************\n" );
	BackgroundModel* bgModel = new BackgroundModel( *Global::negSequenceSet,
													Global::bgModelOrder,
													Global::bgModelAlpha,
													Global::interpolateBG );
	if( Global::saveBgModel ){
		bgModel->write( Global::outputDirectory );
	}

	fprintf( stderr, "\n" );
	fprintf( stderr, "*********************\n" );
	fprintf( stderr, "*   Initial Motif   *\n" );
	fprintf( stderr, "*********************\n" );
	MotifSet motifs;

	int motifNum = ( Global::num ) ? Global::num : motifs.getN();

	if( motifNum > motifs.getN() ){
		std::cout << "--num is larger than the number of the initial motifs." << std::endl;
		exit( -1 );
	}

	for( int n = 0; n < motifNum; n++ ){
		// initialize the model
		Motif* motif = new Motif( *motifs.getMotifs()[n] );

		// train the model with either EM or Gibbs sampling
		ModelLearning model( motif, bgModel );

		if( Global::EM ){

			// learn motifs by EM
			model.EM();

		} else if ( Global::CGS ){

			// learn motifs by collapsed Gibbs sampling
			model.GibbsSampling();

		} else {

			std::cout << "Model is not optimized!\n";

			//exit( -1 );	// allow to do FDR on Initial Model i.e. for PWMs or precalculated BaMMs

		}

		if( Global::generatePseudoSet ){

			// optimize motifs by EM
			model.EM();

			// generate pesudo positive sequence set based on the learned motif
			SeqGenerator seqset( Global::posSequenceSet->getSequences(), model.getMotif() );

			seqset.sample_pseudo_seqset();

			seqset.write( seqset.sample_pseudo_seqset() );

		}

		// write model parameters on the disc
		if( Global::saveBaMMs ){
			model.write( n );
		}

		// write out the learned model
		motif->write( n );

		if ( Global::bammSearch ){
			// 1. generate MN negative sequences of same size and length as posSet
			//    base on 2nd order homogeneous IMM background model
			//    M ~ min{ 10^6/N , 1 }
			SeqGenerator neg_seqset( Global::posSequenceSet->getSequences(), model.getMotif() );
			neg_seqset.sample_negative_seqset();

			// 2. Score each sequence and obtain N- and N+ scores
			//	  sort the scores in descending order
			neg_seqset.score();

			// 3. Sl_lower  = max{ Sn- : Sn- <= Sl }
			// 4. Sl_higher = min{ Sn+ : Sn- >= Sl }
			//    -> as the nearest ranked negative scores.

			// 5. P-value is interpolated between those scores with FPl and FPl+1 false positive counts
			//	  P-value( Sl ) = ( FPl / N- ) + 1 / N-  * ( Sl_higher - Sl + e ) / ( Sl_higher - Sl_lower + e )
			//
			//	  e : very small i.e. 1e-5
			//    --> only works if Sl in between Sl_lower and Sl_higher
			// 6. --> if Sl is lower than any negatve score -> P-value is simply 1
			//
			// 7. --> if Sl is higher than any negative score:
			//	  P-value( Sl ) = 1 / N- exp( - ( Sl - Smax- ) / lambda )
			//	  lambda = 1 / n_top  sum_from n=1 to ntop( Sn- - sntop- )
			//    ntop = min { 100, 0.1 X N- }
			//
			// 8. calculate E-value:
			//	  -> Evalue = N+ * P-Value
		}

		if ( Global::scoreSeqset){
			// score the model on sequence set
			ScoreSeqSet seqset( motif, bgModel, Global::posSequenceSet->getSequences() );
			seqset.score();
			seqset.write( n, Global::scoreCutoff );
		}

		delete motif;
	}

	// evaluate motifs
	if( Global::FDR ){
		fprintf( stderr, "\n" );
		fprintf( stderr, "***********\n" );
		fprintf( stderr, "*   FDR   *\n" );
		fprintf( stderr, "***********\n" );
		for( int n = 0; n < motifNum; n++ ){
			Motif* motif = new Motif( *motifs.getMotifs()[n] );
			FDR fdr( motif );
			fdr.evaluateMotif( n );
			fdr.write( n );
			delete motif;
		}
	}

	fprintf( stderr, "\n" );
	fprintf( stderr, "******************\n" );
	fprintf( stderr, "*   Statistics   *\n" );
	fprintf( stderr, "******************\n" );
	std::cout << "Given alphabet type is " << Alphabet::getAlphabet();
	// for positive sequence set
	std::cout << "\nGiven positive sequence set is " << Global::posSequenceBasename
			<< "\n	"<< Global::posSequenceSet->getN() << " sequences. max.length: " <<
			Global::posSequenceSet->getMaxL() << ", min.length: " <<
			Global::posSequenceSet->getMinL() << "\n	base frequencies:";
	for( int i = 0; i < Alphabet::getSize(); i++ ){
		std::cout << ' ' << Global::posSequenceSet->getBaseFrequencies()[i]
		          << "(" << Alphabet::getAlphabet()[i] << ")";
	}
	// for negative sequence set
	if( Global::negSeqGiven ){
		std::cout << "\nGiven negative sequence set is " << Global::negSequenceBasename
				<< "\n	"<< Global::negSequenceSet->getN() << " sequences. max.length: "
				<< Global::negSequenceSet->getMaxL() << ", min.length: " <<
				Global::negSequenceSet->getMinL() << "\n	base frequencies:";
		for( int i = 0; i < Alphabet::getSize(); i++ )
			std::cout << ' ' << Global::negSequenceSet->getBaseFrequencies()[i]
					  << "(" << Alphabet::getAlphabet()[i] << ")";
	}
	std::cout << "\nGiven initial model is " << Global::initialModelBasename;
	if( Global::FDR ){
		std::cout << "\nGiven folds for FDR estimation: " << Global::cvFold;
	}

	fprintf( stdout, "\n-------------- Runtime: %.2f seconds (%0.2f minutes) --------------\n",
			( ( float )( clock() - t0 ) ) / CLOCKS_PER_SEC, ( ( float )( clock() - t0 ) ) / ( CLOCKS_PER_SEC * 60.0f ) );

	// free memory
	if( bgModel ) delete bgModel;
	Global::destruct();

	return 0;
}
