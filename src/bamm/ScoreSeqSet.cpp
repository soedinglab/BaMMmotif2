/*
 * ScoreSeqSet.cpp
 *
 *  Created on: Dec 13, 2016
 *      Author: wanwan
 */

#include "ScoreSeqSet.h"

#include <float.h>		// -FLT_MAX

ScoreSeqSet::ScoreSeqSet( Motif* motif, BackgroundModel* bg, std::vector<Sequence*> seqSet ){

	for( int k = 0; k < Global::Yk; k++ ){
		Y_.push_back( ipow( Alphabet::getSize(), k ) );
	}

	motif_ = motif;
	bg_ = bg;
	seqSet_ = seqSet;

}


ScoreSeqSet::~ScoreSeqSet(){

}

// store the log odds scores at all positions of each sequence
void ScoreSeqSet::score(){

	int K = Global::modelOrder;
	int W = motif_->getW();

	// pre-calculate log odds scores given motif and bg model
	motif_->calculateS( bg_ );
	float** s = motif_->getS();

	mops_scores_.resize( seqSet_.size() );

	for( size_t n = 0; n < seqSet_.size(); n++ ){

		int 	LW1 = seqSet_[n]->getL() - W + 1;
		int* 	kmer = seqSet_[n]->getKmer();
		float 	maxScore = -FLT_MAX;

		for( int i = 0; i < LW1; i++ ){

			float logOdds = 0.0f;

			for( int j = 0; j < W; j++ ){

				int y = ( kmer[i+j] >= 0 ) ? kmer[i+j] % Y_[K+1] : -1;

				logOdds += ( y >= 0 ) ? s[y][j] : 0;

			}

			// take all the log odds scores for MOPS model:
			mops_scores_[n].push_back( logOdds );

			// take the largest log odds score for ZOOPS model:
			maxScore = ( logOdds > maxScore ) ? logOdds : maxScore;
		}
		zoops_scores_.push_back( maxScore );
	}
	// fill up ALL and Max
	for( size_t n = 0; n < seqSet_.size(); n++ ){
		ScoreAll_.insert( std::end( ScoreAll_ ), std::begin( mops_scores_[n] ), std::end( mops_scores_[n] ) );
	}
}


// compute p_values for motif scores based on negative sequeunce scores
void ScoreSeqSet::calcPvalues( std::vector<float> neg_scores ){

	std::vector<std::vector<float>> pos_scores = mops_scores_;
	mops_p_values_.resize( seqSet_.size() );
	mops_e_values_.resize( seqSet_.size() );

	int FPl = 0;
	float SlLower, SlHigher, p_value, Sl;
	float lambda = 0;
	float eps = (float) 1e-5;
	int negN = (int) neg_scores.size();
	int posN = (int) ScoreAll_.size();
	int nTop = std::min( 100, int( 0.1*negN ));

	for( int n = 0; n < nTop; n++ ){
		//fprintf(stderr, "%f\n", neg_scores[n]);
		lambda =+ (neg_scores[n] - neg_scores[nTop]);
	}
	lambda = lambda / (float)nTop;


	//fprintf(stderr, "\n negN = %d  posN = %d  nTop = %d  lambda = %f \n", negN, posN, nTop, lambda);

	for( size_t n = 0; n < seqSet_.size(); n++ ){
		int seqlen = seqSet_[n]->getL();
		if( !Global::ss ){
			seqlen = seqlen / 2;
		}
		int LW1 = seqSet_[n]->getL() - motif_->getW() + 1;

		for( int i = 0; i < LW1; i++ ){
			Sl = pos_scores[n][i];

			//fprintf(stderr, "\n Sl_prev = %f  Sl = %f \n", Sl_prev, Sl);

			FPl = (int) count_if(neg_scores.begin(), neg_scores.end(), [Sl](float n) { return n >= Sl; } );

			//fprintf(stderr, "\n FPl =  %d  \n", FPl);

			// Sl is lower than worst negScore
			if( FPl == negN ){
				p_value = (float) 1;
				mops_p_values_[n].push_back( (float)p_value );
				mops_e_values_[n].push_back( (float)p_value * (float)posN );
				//fprintf(stderr, "p_value = %f \n", (float)p_value );
				continue;
			}
			// Sl is higher than best negScore
			if ( FPl == 0 ){
				p_value =  float(1) / (float)negN * (float) exp( - ( Sl - neg_scores[0] ) / lambda );
				mops_p_values_[n].push_back( (float)p_value );
				mops_e_values_[n].push_back( (float)p_value * (float)posN );
				//fprintf(stderr, "p_value = %f \n", p_value );
				continue;
			}else{
				// SlHigher and SlLower can be defined
				SlHigher = neg_scores[FPl-1];
				SlLower = neg_scores[FPl];
				p_value = float(FPl) / float(negN) + float(1) /float(negN) * ( SlHigher - Sl + eps ) / ( SlHigher - SlLower + eps );
				mops_p_values_[n].push_back( (float)p_value );
				mops_e_values_[n].push_back( (float)p_value * (float)posN );
				//fprintf(stderr, "FPl / negN   = %f \n", (float) FPl / (float)negN  );
				//fprintf(stderr, "  1 / negN   = %f \n", (float) 1 / (float)negN );
				//fprintf(stderr, "SlH - Sl +e  = %f \n",  SlHigher - Sl + eps );
				//fprintf(stderr, "SlH - SlL +e = %f \n",  SlHigher - SlLower + eps );
				//fprintf(stderr, "Bruch        = %f \n", ( SlHigher - Sl + eps ) / ( SlHigher - SlLower + eps ) );

				//fprintf(stderr, "SlHigher = %f \n", SlHigher );
				//fprintf(stderr, "SlLower = %f \n", SlLower );
				//fprintf(stderr, "p_value = %f \n", p_value );
				continue;
			}
		}
	}
}

std::vector<float> ScoreSeqSet::getScoreAll(){
	return ScoreAll_;
}

std::vector<std::vector<float>> ScoreSeqSet::getMopsScores(){
	return mops_scores_;
}

std::vector<float> ScoreSeqSet::getZoopsScores(){
	return zoops_scores_;
}


void ScoreSeqSet::write( int N, float cutoff ){

	/**
	 * save log odds scores in one flat file:
	 * posSequenceBasename.logOdds
	 */

	bool 	first_hit = true;
	int 	end; 				// end of motif match

	std::string opath = std::string( Global::outputDirectory ) + '/'
			+ Global::posSequenceBasename +  "_motif_" + std::to_string( N+1 ) + ".logOdds";

	std::ofstream ofile( opath );

	for( size_t n = 0; n < seqSet_.size(); n++ ){
		first_hit = true;
		int seqlen = seqSet_[n]->getL();
		if( !Global::ss ){
			seqlen = seqlen / 2;
		}
		int LW1 = seqSet_[n]->getL() - motif_->getW() + 1;
		for( int i = 0; i < LW1; i++ ){

			if( mops_scores_ [n][i] > cutoff ){
				if( first_hit ){
					// >header:sequence_length
					ofile << '>' << seqSet_[n]->getHeader() <<  ':' << seqlen << std::endl;
					first_hit = false;
				}
				// start:end:score:strand:sequence_matching
				end = i + motif_->getW()-1;

				ofile << i << ':' << end << ':' << std::setprecision( 3 ) << mops_scores_[n][i] << ':' <<
						( ( i < seqlen ) ? '+' : '-' ) << ':' ;
				for( int m = i; m <= end; m++ ){
					ofile << Alphabet::getBase( seqSet_[n]->getSequence()[m] );
				}
				ofile << std::endl;
			}
		}

		//ofile << std::endl;
	}

}

void ScoreSeqSet::writePvalues( int N, float cutoff ){

	/**
	 * save log odds scores in one flat file:
	 * posSequenceBasename.logOdds
	 */

	bool 	first_hit = true;
	int 	end; 				// end of motif match

	std::string opath = std::string( Global::outputDirectory ) + '/'
			+ Global::posSequenceBasename +  "_motif_" + std::to_string( N+1 ) + ".scores";

	std::ofstream ofile( opath );

	for( size_t n = 0; n < seqSet_.size(); n++ ){
		first_hit = true;
		int seqlen = seqSet_[n]->getL();
		if( !Global::ss ){
			seqlen = seqlen / 2;
		}
		int LW1 = seqSet_[n]->getL() - motif_->getW() + 1;
		for( int i = 0; i < LW1; i++ ){

			if( mops_p_values_ [n][i] < cutoff ){
				if( first_hit ){
					// >header:sequence_length
					ofile << '>' << seqSet_[n]->getHeader() <<  ':' << seqlen << std::endl;
					first_hit = false;
				}
				// start:end:score:strand:sequence_matching
				end = i + motif_->getW()-1;

				ofile << i << ':' << end << ':' << std::setprecision( 3 ) << mops_scores_[n][i] << ':' << std::setprecision( 3 ) << mops_p_values_[n][i] << ':'  << std::setprecision( 3 ) << mops_e_values_[n][i] << ':'<<
						( ( i < seqlen ) ? '+' : '-' ) << ':' ;
				for( int m = i; m <= end; m++ ){
					ofile << Alphabet::getBase( seqSet_[n]->getSequence()[m] );
				}
				ofile << std::endl;
			}
		}
	}
}
