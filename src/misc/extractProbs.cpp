//
// Created by wanwan on 25.06.18.
// This script is aimed to extract model probabilities from
// conditional probabilities in BaMM format
// Input: file with suffix .ihbcp & .hbcp
// Output: file with suffix .hbcp & .hbp
//

#include "../init/Alphabet.h"
#include "../init/BackgroundModel.h"
#include "../init/Motif.h"

void constructBaMM( char* indir, char* odir, float** v_bg, size_t k_bg) {
    // each BaMM file contains one optimized motif model
    // read file to calculate motif length
    std::ifstream file;
    file.open(indir, std::ifstream::in);

    if (!file.good()) {
        std::cerr << "Error: Cannot open BaMM file: " << indir << std::endl;
        exit(1);
    } else {

        size_t model_length = 0;
        size_t model_order = 0;
        size_t check_lines = 0;
        std::string line;

        // read file to calculate motif length, order and alphabet size
        while (getline(file, line)) {
            if (line.empty()) {
                // count the number of empty lines as the motif length
                model_length++;
                // check if the input BaMM file has the right format
                if (model_length > 1 and check_lines != model_order) {
                    std::cerr << "This is not a BaMM-format file: " << indir << std::endl;
                    exit(1);
                }
                check_lines = 0;
            } else if (model_length == 0) {
                // count the lines of the first position
                model_order++;
            } else {
                check_lines++;
            }
        }

        // adjust model order, extra 1
        model_order -= 1;
        if (model_order > 8) {
            std::cerr << "The input BaMM model order is too high: " << indir << std::endl;
            exit(1);
        }

        std::vector<float> 	alphas( model_order+1, 1.f );

        // construct an initial motif
        Motif *motif = new Motif(model_length, model_order, alphas, v_bg, k_bg, 1);

        // initialize motif from file
        motif->initFromBaMM(indir, 0, 0);

        // write out the foreground model
        motif->write( odir, baseName( indir ) );

        if( motif )		delete motif;

    }
}

int mainExtractPos( int nargs, char* args[] ){

    /**
     * read in input files
     */
    if( nargs < 3 ) {
        std::cerr << "Error: Arguments are missing!" << std::endl
                  << "Usage: extractProbs <outdir> <bamm foreground model(.ihbcp)> <bamm background model(.hbcp)>"
                  << std::endl;
        exit( 1 );
    }

    // read in the output directory and create it
    char* outputDirectory = args[1];
    createDirectory( outputDirectory );

    // read in foreground model
    char* BaMMFilename = args[2];

    // read in background model
    char* BgFilename = args[3];

    // define Alphabet type, just for getting the size of alphabet table
    char* alphabetType = new char[9];
    strcpy( alphabetType, "STANDARD" );
    Alphabet::init( alphabetType );

    // build background model
    std::string BgBasename = baseName( BgFilename );
    BackgroundModel* bgModel;
    bgModel = new BackgroundModel( BgFilename );
    // save background model
    bgModel->write( outputDirectory, BgBasename );

    // construct foreground model and save it
    constructBaMM(BaMMFilename, outputDirectory, bgModel->getV(), bgModel->getOrder());

    return 0;
}

