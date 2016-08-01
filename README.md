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

      --negSequenceSet <FILEPATH>
          FASTA file with negative/background sequences used to learn the
          (homogeneous) background BaMM. If not specified, the background BaMM
          is learned from the positive sequences.

      --reverseComp
          Search motifs on both strands (positive sequences and reverse
          complements). This option is e.g. recommended when using sequences
          derived from ChIP-seq experiments.

  Options to initialize a single BaMM from file

      --bindingSiteFile <FILEPATH>
          File with binding sites of equal length (one per line).

      --markovModelFile <FILEPATH>
          File with BaMM probabilities as obtained from BaMM!motif (omit
          filename extension).

  Options to initialize one or more BaMMs from XXmotif PWMs

      --minPWMs <INTEGER>
          Minimum number of PWMs. The options --maxPValue and --minOccurrence
          are ignored. The default is 1.

      --maxPWMs <INTEGER>
          Maximum number of PWMs.

      --maxPValue <FLOAT>
          Maximum p-value of PWMs. This filter is not applied to the top
          minimum number of PWMs (see --minPWMs). The default is 1.0.

      --minOccurrence <FLOAT>
          Minimum fraction of sequences that contain the motif. This filter is
          not applied to the top minimum number of PWMs (see --minPWMs). The
          default is 0.05.

      --rankPWMs <INTEGER> [<INTEGER>...]
          PWM ranks in XXmotif results. The former options to initialize BaMMs
          from PWMs are ignored.

  Options for (inhomogeneous) motif BaMMs

      -k <INTEGER>
          Order. The default is 2.

      -a|--alpha <FLOAT> [<FLOAT>...]
          Order-specific prior strength. The default is 1.0 (for k = 0) and
          20 x 3^(k-1) (for k > 0). The options -b and -g are ignored.

      -b|--beta <FLOAT>
          Calculate order-specific alphas according to beta x gamma^(k-1) (for
          k > 0). The default is 20.0.

      -g|--gamma <FLOAT>
          Calculate order-specific alphas according to beta x gamma^(k-1) (for
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

  EM options

      -q <FLOAT>
          Prior probability for a positive sequence to contain a motif. The
          default is 0.9.

      -e|--epsilon <FLOAT>
          The EM algorithm is deemed to be converged when the sum over the
          absolute differences in BaMM probabilities from successive EM rounds
          is smaller than epsilon. The default is 0.001.

  XXmotif options

      --XX-ZOOPS
          Use the zero-or-one-occurrence-per-sequence model (default).

      --XX-MOPS
          Use the multiple-occurrence-per-sequence model.

      --XX-OOPS
          Use the one-occurrence-per-sequence model.

      --XX-seeds ALL|FIVEMERS|PALINDROME|TANDEM|NOPALINDROME|NOTANDEM
          Define the nature of seed patterns. The default is to start using ALL
          seed pattern variants.

      --XX-gaps 0|1|2|3
          Maximum number of gaps used for seed patterns. The default is 0.

      --XX-pseudoCounts <FLOAT>
          Percentage of pseudocounts. The default is 10.0.

      --XX-mergeMotifsThreshold LOW|MEDIUM|HIGH
          Define the similarity threshold used to merge PWMs. The default is to
          merge PWMs with LOW similarity in order to reduce runtime.

      --XX-maxPositions <INTEGER>
          Limit the number of motif positions to reduce runtime. The default is
          17.

      --XX-noLengthOptimPWMs
          Omit the length optimization of PWMs.

      --XX-K <INTEGER>
          Order of the (homogeneous) background BaMM. The default is either 2
          (when learned on positive sequences) or 8 (when learned on background
          sequences).

      --XX-A <FLOAT>
          Prior strength of the (homogeneous) background BaMM. The default is
          10.0.

      --XX-jumpStartPatternStage <STRING>
          Jump-start pattern stage using an IUPAC pattern string.

      --XX-jumpStartPWMStage <FILEPATH>
          Jump-start PWM stage reading in a PWM from file.

      --XX-localization
          Calculate p-values for positional clustering of motif occurrences in
          positive sequences of equal length. Improves the sensitivity to find
          weak but positioned motifs.

      --XX-localizationRanking
          Rank motifs according to localization statistics.

      --XX-downstreamPositions <INTEGER>
          Distance between the anchor position (e.g. the transcription start
          site) and the last positive sequence nucleotide. Corrects motif
          positions in result plots. The default is 0.

      --XX-batch
          Suppress progress bars.

  Options to score sequences

      --scorePosSequenceSet
          Score positive (training) sequences with optimized BaMMs.

      --scoreNegSequenceSet
          Score background (training) sequences with optimized BaMMs.

      --scoreTestSequenceSet <FILEPATH> [<FILEPATH>...]
          Score test sequences with optimized BaMMs. Test sequences can be
          provided in a single or multiple FASTA files.

  Output options

      --saveInitBaMMs
          Write initialized BaMM(s) to disk.

      --saveBaMMs
          Write optimized BaMM(s) to disk.

      --verbose
          Verbose terminal printouts.

      -h, --help
          Printout this help.

## BaMM flat file format

BaMMs are written to flat file when invoking BaMM!motif with the output option `--saveInitBaMMs` and/or `--saveBaMMs`. In this case, BaMM!motif generates three files for each (inhomogeneous) BaMM &ndash; one containing the probabilities (filename extension: probs), one containing the conditional probabilities (filename extension: conds), and one containing the background frequencies of mononucleotides in the positive sequences (file extension: freqs). The format is the same for the first two. While blank lines separate BaMM positions, lines 1 to *k*+1 of each BaMM position contain the (conditional) probabilities for order 0 to order *k*. For instance, the format for a BaMM of order 2 and length *W* is as follows:

Filename extension: probs

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

Filename extension: conds

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

Filename extension: freqs

P(A) P(C) P(G) P(T)<br>

Note that contexts are restricted to the binding site. For instance, P<sub>1</sub>(G|AC) and P<sub>2</sub>(G|AC) are defined as P<sub>1</sub>(G) and P<sub>2</sub>(G|C), respectively.

In addition, BaMM!motif generates three files for the (homogeneous) background BaMM &ndash; one containing the probabilities (filename extension: probsBg), one containing the conditional probabilities (filename extension: condsBg), and one containing the background frequencies of mononucleotides (file extension: freqs). For instance, the format for a background BaMM of order 2 is as follows:

Filename extension: probsBg

P(A) P(C) P(G) P(T)<br>
P(AA) P(AC) P(AG) P(AT) P(CA) P(CC) P(CG) ... P(TT)<br>
P(AAA) P(AAC) P(AAG) P(AAT) P(ACA) P(ACC) P(ACG) ... P(TTT)<br>

Filename extension: condsBg

P(A) P(C) P(G) P(T)<br>
P(A|A) P(C|A) P(G|A) P(T|A) P(A|C) P(C|C) P(G|C) ... P(T|T)<br>
P(A|AA) P(C|AA) P(G|AA) P(T|AA) P(A|AC) P(C|AC) P(G|AC) ... P(T|TT)<br>

Filename extension: freqsBg

P(A) P(C) P(G) P(T)<br>

Note that the background frequencies of mononucleotides are identical to the probabilities of mononucleotides in the other two files.

## How to plot BaMM logos?

R scripts are provided in directory R to plot the BaMM logo from a BaMM flat file. To create a BaMM logo, edit the parameter setting in `plotBaMM.wrapper.R` and source the code in the R session using

    source( "plotBaMM.wrapper.R" )

Please find comments on available plotting options in the wrapper.

## License
BaMM!motif is released under the GNU General Public License v3 or later. See LICENSE for more details.
