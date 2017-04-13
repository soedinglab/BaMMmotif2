#include <time.h>		// time(), clock_t, clock, CLOCKS_PER_SEC
#include <stdio.h>
#include <algorithm>    // std::min

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
			BackgroundModel* bg = bgModel;
			// 0. Depending on motif input, and learning, define bgModel
			// use bgModel generated from input sequences wehn prediction is turned on
			if( ! Global::EM and ! Global::CGS ){
				// use provided bgModelFile if initialized with bamm format
				if( Global::BaMMFilename != NULL ){
					if( Global::bgModelFile == NULL ){
						std::cout << "No background Model File provided for initial search motif!\n";
						exit(-1);
					}
					bg = new BackgroundModel( Global::bgModelFile );
				}
				// use bgModel generated when reading in PWM File
				if( Global::PWMFilename != NULL ){
					bg = new BackgroundModel( Global::PWMFilename , 0, 1);
				}
				else{
					std::cout << "No background Model found for provided initial search motif!\n";
					exit(-1);
				}
			}


			int posN = Global::posSequenceSet->getN();
			int M = std::min(int(pow(10,6)/posN),1);
			int negN = posN * M;

			// 1. generate MN negative sequences of same size and length as posSet
			//    base on 2nd order homogeneous IMM background model
			//    M ~ min{ 10^6/N , 1 }
			SeqGenerator neg_seqs( Global::posSequenceSet->getSequences(), model.getMotif() );
			neg_seqs.sample_negative_seqset( negN );

			// 2. Score each sequence and obtain N- and N+ scores
			//	  sort the scores in descending order
			ScoreSeqSet neg_seqset( model.getMotif(), bg, neg_seqs.getSeqs());
			neg_seqset.score();
			ScoreSeqSet pos_seqset( model.getMotif(), bg, Global::posSequenceSet->getSequences());
			pos_seqset.score();


			// obtain p_value and e_value for each positive_mopsscore
			std::vector<float> pos_scores = pos_seqset.getScoreAll();
			std::vector<float> neg_scores = neg_seqset.getScoreAll();

			std::sort( pos_scores.begin(), pos_scores.end(), std::greater<float>() );
			std::sort( neg_scores.begin(), neg_scores.end(), std::greater<float>() );



			int FPl = 0;
			float Sl, SlLower, SlHigher, p_value;
			float lambda = 0;
			float eps = (float) 1e-5;
			int nTop = std::min( 100, int( 0.1*negN ));
			std::vector<float> p_values;
			std::vector<float> e_values;

			for( int n = 0; n < nTop; n++ ){
				lambda =+ (neg_scores[n] - neg_scores[nTop]);
			}
			lambda = lambda / (float)nTop;

			for( int l = 0; l < posN; l++ ){
				Sl = pos_scores[l];
				while( Sl <= neg_scores[FPl] && FPl < negN){
					SlHigher = neg_scores[FPl];
					FPl++;
				}
				// Sl is lower than worst negScore
				if( FPl == negN ){
					p_value = (float) 1;
					p_values.push_back( p_value );
					e_values.push_back( p_value * (float)posN );
					continue;
				}
				// Sl is higher than best negScore
				if ( FPl == 0 ){
					p_value =  float(1 / negN) * (float) exp( - ( Sl - neg_scores[0] ) / lambda );
					p_values.push_back( p_value );
					e_values.push_back( p_value * (float)posN );
					continue;
				}else{
					// SlHigher and SlLower can be defined
					SlLower = neg_scores[FPl];
					p_value = float(FPl / negN) + float(1 / negN) * ( SlHigher - Sl + eps ) / ( SlHigher - SlLower + eps );
					p_values.push_back( p_value );
					e_values.push_back( p_value * (float)posN );
					continue;
				}
			}
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
