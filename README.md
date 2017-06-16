# BaMM!motif

**Ba**yesian **M**arkov **M**odel **motif** discovery software.

## Requirements
To compile from source you need

  * [GCC](https://gcc.gnu.org/) compiler 4.7 or later
  * [CMake](http://cmake.org/) 2.8.11 or later

To plot BaMM logos you need

  * [R](https://cran.r-project.org/) 2.14.1 or later

## How to compile BaMM!motif?

      mkdir build
      cd build
      cmake ..
      make

## How to use BaMM!motif from the command line?
SYNOPSIS

      BaMMmotif DIRPATH FILEPATH [OPTIONS]

DESCRIPTION

      Bayesian Markov Model motif discovery software.

      DIRPATH
          Output directory for the results.

      FILEPATH
          FASTA file with positive sequences of equal length.

OPTIONS

Sequence options

      --alphabet <STRING>
          STANDARD.  	For alphabet type ACGT, default setting;
          METHYLC.   	For alphabet type ACGTM;
          HYDROXYMETHYLC.  For alphabet type ACGTH;
          EXTENDED.  	For alphabet type ACGTMH.
      
      --ss
          Search motif only on single strand strands (positive sequences).
          This option is not recommended for analyzing ChIP-seq data.
          By default, BaMM searches motifs on both strands.
          
      --negSeqSet <FILEPATH>
          FASTA file with negative/background sequences used to learn the
          (homogeneous) background BaMM. If not specified, the background BaMM
          is learned from the positive sequences.

  Options to initialize BaMM(s) from file

      --bindingSiteFile <FILEPATH>
          File with binding sites of equal length (one per line).
      
      --PWMFile <STRING>
          File that contains position weight matrices (PWMs).
      
      --BaMMFile <STRING>
          File that contains a model in bamm file format.

      --num <INTEGER>
          Number of models to be learned by BaMM!motif, specific for PWMs.

  Options for the (inhomogeneous) motif BaMMs

      -k|--order <INTEGER>
          Model order. The default is 2.

      -a|--alpha <FLOAT> [<FLOAT>...]
          Order-specific prior strength. The default is 1.0 (for k = 0) and
          beta x gamma^k (for k > 0). The options -b and -g are ignored.

      -b|--beta <FLOAT>
          Calculate order-specific alphas according to beta x gamma^k (for
          k > 0). The default is 7.0.

      -g|--gamma <FLOAT>
          Calculate order-specific alphas according to beta x gamma^k (for
          k > 0). The default is 3.0.

      --extend <INTEGER>{1,2}
          Extend BaMMs by adding uniformly initialized positions to the left
          and/or right of initial BaMMs. Invoking e.g. with --extend 0 2 adds
          two positions to the right of initial BaMMs. Invoking with --extend 2
          adds two positions to both sides of initial BaMMs. By default, BaMMs
          are not being extended.

  Options for the (homogeneous) background BaMM

      -K <INTEGER>
          Order. The default is 2.

      -A|--Alpha <FLOAT>
          Prior strength. The default is 10.0.
      
      --bgModelFile <STRING>
          Read in background model from a bamm-formatted file. 

  EM options

      --EM
          Triggers Expectation Maximization (EM) algorithm.
      
      -q <FLOAT>
          Prior probability for a positive sequence to contain a motif. The
          default is 0.9.

      -e|--epsilon <FLOAT>
          The EM algorithm is deemed to be converged when the sum over the
          absolute differences in BaMM probabilities from successive EM rounds
          is smaller than epsilon. The default is 0.001.

  Gibbs sampling options:

      --CGS
          Triggers Collapsed Gibbs Sampling (CGS) algorithm.
      
      --maxCGSIterations <INTEGER> 
          Limit the number of CGS iterations.
          It should be larger than 5 and defaults to 100.

  Options for model evaluation:
      
      --FDR
          Triggers False-Discovery-Rate (FDR) estimation.
        
      -m|--mFold <INTEGER>
          Number of negative sequences as multiple of positive sequences.
          The default is 10.
      
      -n, --cvFold <INTEGER>
          Fold number for cross-validation. 
          The default is 5, which means the training set is 4-fold of the test set.
          
      -s, --sOrder <INTERGER>
          The order of k-mer for sampling pseudo/negative set. The default is 2.

  Output options

      --saveBaMMs
          Write optimized BaMM(s) to disk.

      --saveInitBaMMs
          Write initialized BaMM(s) to disk.
      
      --saveBgModel
          Write background model to disk.
          
      --verbose
          Verbose terminal printouts.

      -h, --help
          Printout this help.

## BaMM flat file format

BaMMs are written to flat file when invoking BaMM!motif with the output option `--saveInitBaMMs` and/or `--saveBaMMs`. In this case, BaMM!motif generates two files for each inhomogeneous BaMM &ndash; one containing the probabilities (filename extension: `ihbp`) and one containing the conditional probabilities (filename extension: `.ihbcp`). The format is the same for these two files. While blank lines separate BaMM positions, lines 1 to *k*+1 of each BaMM position contain the (conditional) probabilities for order 0 to order *k*. For instance, the format for a BaMM of order 2 and length *W* is as follows:

Filename extension: .ihbp

P<sub>1</sub>(A) P<sub>1</sub>(C) P<sub>1</sub>(G) P<sub>1</sub>(T)<br>
P<sub>1</sub>(AA) P<sub>1</sub>(AC) P<sub>1</sub>(AG) P<sub>1</sub>(AT) P<sub>1</sub>(CA) P<sub>1</sub>(CC) P<sub>1</sub>(CG) ... P<sub>1</sub>(TT)<br>
P<sub>1</sub>(AAA) P<sub>1</sub>(AAC) P<sub>1</sub>(AAG) P<sub>1</sub>(AAT) P<sub>1</sub>(ACA) P<sub>1</sub>(ACC) P<sub>1</sub>(ACG) ... P<sub>1</sub>(TTT)<br>

P<sub>2</sub>(A) P<sub>2</sub>(C) P<sub>2</sub>(G) P<sub>2</sub>(T)<br>
P<sub>2</sub>(AA) P<sub>2</sub>(AC) P<sub>2</sub>(AG) P<sub>2</sub>(AT) P<sub>2</sub>(CA) P<sub>2</sub>(CC) P<sub>2</sub>CG) ... P<sub>2</sub>(TT)<br>
P<sub>2</sub>(AAA) P<sub>2</sub>(AAC) P<sub>2</sub>(AAG) P<sub>2</sub>(AAT) P<sub>2</sub>(ACA) P<sub>2</sub>(ACC) P<sub>2</sub>(ACG) ... P<sub>2</sub>(TTT)<br>
...

P<sub>W</sub>(A) P<sub>W</sub>(C) P<sub>W</sub>(G) P<sub>W</sub>(T)<br>
P<sub>W</sub>(AA) P<sub>W</sub>(AC) P<sub>W</sub>(AG) P<sub>W</sub>(AT) P<sub>W</sub>(CA) P<sub>W</sub>(CC) P<sub>W</sub>CG) ... P<sub>W</sub>(TT)<br>
P<sub>W</sub>(AAA) P<sub>W</sub>(AAC) P<sub>W</sub>(AAG) P<sub>W</sub>(AAT) P<sub>W</sub>(ACA) P<sub>W</sub>(ACC) P<sub>W</sub>(ACG) ... P<sub>W</sub>(TTT)<br>

Filename extension: .ihbcp

P<sub>1</sub>(A) P<sub>1</sub>(C) P<sub>1</sub>(G) P<sub>1</sub>(T)<br>
P<sub>1</sub>(A|A) P<sub>1</sub>(C|A) P<sub>1</sub>(G|A) P<sub>1</sub>(T|A) P<sub>1</sub>(A|C) P<sub>1</sub>(C|C) P<sub>1</sub>(G|C) ... P<sub>1</sub>(T|T)<br>
P<sub>1</sub>(A|AA) P<sub>1</sub>(C|AA) P<sub>1</sub>(G|AA) P<sub>1</sub>(T|AA) P<sub>1</sub>(A|AC) P<sub>1</sub>(C|AC) P<sub>1</sub>(G|AC) ... P<sub>1</sub>(T|TT)<br>

P<sub>2</sub>(A) P<sub>2</sub>(C) P<sub>2</sub>(G) P<sub>2</sub>(T)<br>
P<sub>2</sub>(A|A) P<sub>2</sub>(C|A) P<sub>2</sub>(G|A) P<sub>2</sub>(T|A) P<sub>2</sub>(A|C) P<sub>2</sub>(C|C) P<sub>2</sub>(G|C) ... P<sub>2</sub>(T|T)<br>
P<sub>2</sub>(A|AA) P<sub>2</sub>(C|AA) P<sub>2</sub>(G|AA) P<sub>2</sub>(T|AA) P<sub>2</sub>(A|AC) P<sub>2</sub>(C|AC) P<sub>2</sub>(G|AC) ... P<sub>2</sub>(T|TT)<br>
...

P<sub>W</sub>(A) P<sub>W</sub>(C) P<sub>W</sub>(G) P<sub>W</sub>(T)<br>
P<sub>W</sub>(A|A) P<sub>W</sub>(C|A) P<sub>W</sub>(G|A) P<sub>W</sub>(T|A) P<sub>W</sub>(A|C) P<sub>W</sub>(C|C) P<sub>W</sub>(G|C) ... P<sub>W</sub>(T|T)<br>
P<sub>W</sub>(A|AA) P<sub>W</sub>(C|AA) P<sub>W</sub>(G|AA) P<sub>W</sub>(T|AA) P<sub>W</sub>(A|AC) P<sub>W</sub>(C|AC) P<sub>W</sub>(G|AC) ... P<sub>W</sub>(T|TT)<br>


In addition, BaMM!motif generates two files for the homogeneous background BaMM &ndash; one containing the probabilities (filename extension: `.hbp`) and the other containing the conditional probabilities (filename extension: `.hbcp`). For instance, the format for a background BaMM of order 2 is as follows:

Filename extension: .hbp

P(A) P(C) P(G) P(T)<br>
P(AA) P(AC) P(AG) P(AT) P(CA) P(CC) P(CG) ... P(TT)<br>
P(AAA) P(AAC) P(AAG) P(AAT) P(ACA) P(ACC) P(ACG) ... P(TTT)<br>

Filename extension: .hbcp

P(A) P(C) P(G) P(T)<br>
P(A|A) P(C|A) P(G|A) P(T|A) P(A|C) P(C|C) P(G|C) ... P(T|T)<br>
P(A|AA) P(C|AA) P(G|AA) P(T|AA) P(A|AC) P(C|AC) P(G|AC) ... P(T|TT)<br>


## How to plot BaMM logos?

R scripts are provided in directory R to plot the BaMM logo from a BaMM flat file. To create a BaMM logo, edit the parameter setting in `plotBaMM.wrapper.R` and source the code in the R session using

    source( "plotBaMM.wrapper.R" )

Please find comments on available plotting options in the wrapper.

## License
BaMM!motif is released under the GNU General Public License v3 or later. See LICENSE for more details.
