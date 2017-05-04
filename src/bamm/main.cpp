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
			fprintf( stderr, "\n" );
    		fprintf( stderr, "**********************\n" );
			fprintf( stderr, "*   Score Sequences  *\n" );
			fprintf( stderr, "**********************\n" );

			BackgroundModel* bg = bgModel;
			// 0. Depending on motif input, and learning, define bgModel
			// use bgModel generated from input sequences wehn prediction is turned on
			std::cout << "bgmodel is copied\n";

			if( ! Global::EM and ! Global::CGS ){

				//std::cout << "no Optimization takes place EM! CGS!\n";

				// use provided bgModelFile if initialized with bamm format
				if( Global::BaMMFilename != NULL ){
					//std::cout << "init with BaMMfile \n";
					if( Global::bgModelFile == NULL ){
						//std::cout << "bgModel not provided! \n";
						std::cout << "No background Model File provided for initial search motif!\n";
						exit(-1);
					}
					bg = new BackgroundModel( Global::bgModelFile );
				}
				else{
					//std::cout << "not optimized with BammFile\n";
					// use bgModel generated when reading in PWM File
					if( Global::PWMFilename != NULL ){
						//std::cout << "Init with PWM instead \n";
						bg = new BackgroundModel( Global::PWMFilename , 0, 1);
						// this means that also the global motif order needs to be adjusted;
						Global::modelOrder = 0;
					}
					else{
						//std::cout << "Other init, no implementation yet\n";
						std::cout << "No background Model found for provided initial search motif!\n";
						exit(-1);
					}
				}

			}

			// generate M * Npos negative sequences
			int M = std::min(int(pow(10,6)/Global::posSequenceSet->getN()),5);

			SeqGenerator neg_seqs( Global::posSequenceSet->getSequences(), model.getMotif() );
			neg_seqs.sample_negative_seqset( M );

			// generate and score negative sequence set
			ScoreSeqSet neg_seqset( model.getMotif(), bg, neg_seqs.getSeqs());
			neg_seqset.score();

			// score positive sequence set
			ScoreSeqSet pos_seqset( model.getMotif(), bg, Global::posSequenceSet->getSequences());
			pos_seqset.score();

			// collapse and rank negative scores
			std::vector<float> neg_scores = neg_seqset.getScoreAll();
			std::sort( neg_scores.begin(), neg_scores.end(), std::greater<float>() );

			fprintf( stderr, "\n" );
    		fprintf( stderr, "***************************\n" );
    		fprintf( stderr, "*   Evaluate Occurrences  *\n" );
			fprintf( stderr, "***************************\n" );

			// calculate p- and e-values for positive scores based on negative scores
			pos_seqset.calcPvalues(neg_scores);

			// print out scores and p-/e-values larger than p-values based cutoff
			pos_seqset.writePvalues( n, Global:: pvalCutoff );
		}

		if ( Global::scoreSeqset){
			// score the model on sequence set
			ScoreSeqSet seqset( motif, bgModel, Global::posSequenceSet->getSequences() );
			seqset.score();
			seqset.writePvalues( n, Global::scoreCutoff );
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
