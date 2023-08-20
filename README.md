- [eGADA](#egada)
- [Prerequisites to build eGADA](#prerequisites-to-build-egada)
- [Build the eGADA using g++](#build-the-egada-using-g)
- [To run](#to-run)
  - [The C++ binary](#the-c-binary)
  - [Calling the dynamic library from Python](#calling-the-dynamic-library-from-python)
  - [Docker image](#docker-image)
- [References](#references)

# eGADA

eGADA (enhanced Genomic Alteration Detection Algorithm) is an enhanced version of GADA (Pique-Regi 2007, Pique-Regi 2008).

It is a fast segmentation algorithm utilizing the Sparse Bayesian Learning (or Relevance Vector Machine) technique from Tipping 2001. It can be applied to array intensity data, sequencing coverage data, or any sequential data that displays characteristics of step-wise functions.

eGADA (Huang 2023) vs GADA (Pique-Regi 2007):

- Use a customized Red-Black tree to expedite the final backward elimination step.
- Code in C++, not C.
- Use Boost libraries extensively.
- More friendly help and commandline-argument processing.
- More friendly input and output formats.
- A dynamic library eGADA.so (packaged via Boost.Python) that offers API to Python.
- Minor changes/bugfixes.

# Prerequisites to build eGADA

- libboost-program-options-dev: Boost program-options dev library
- libboost-iostreams-dev: Boost iostreams dev library
- libboost-python-dev: Boost python dev library


# Build the eGADA using g++


```bash
cd src;

make
```


# To run

## The C++ binary

```bash
./src/eGADA -i ../data/input.txt -o ../data/output2.txt
```
- The input is a single-column plain text file.
- The output is a 4-column tsv file:
  - Start is the starting index (1-based) of the segment.
  - Stop is the ending index (1-based and inclusive) of the segment.
  - Length is the number of data points/bins/probes included in the segment.
  - Ampl (Amplitude) is the average amplitude/coverage of the segment.

```
# eGADA: enhanced Genome Alteration Detection Algorithm
# Authors: Yu Huang polyactis@gmail.com, Roger Pique-Regi piquereg@usc.edu
# Parameters: a=0.2,T=5,MinSegLen=0,sigma2=0.207949,BaseAmp=0, convergenceDelta=1e-08, maxNoOfIterations=50000, convergenceMaxAlpha=1e+08, convergenceB=1e-20.
# Reading M=80000 probes in input file
# Overall mean 0.00135799
# Sigma^2=0.207949
# Convergence: delta=9.88473e-09 after 924 EM iterations.
# Found 1915 breakpoints after SBL
# Kept 442 breakpoints after BE
Start   Stop    Length  Ampl
1       25      25      0.881578
26      75      50      -1.08218
76      125     50      1.04644
126     175     50      -0.973855
176     226     51      0.912763
227     275     49      -0.954001
276     325     50      0.963878
326     375     50      -1.16553
376     427     52      0.948964
...
```

- Type ```eGADA``` or ```eGADA -h ``` for more help.

```bash
yh@fusilier:~/src/eGADA/src$ ./eGADA
program name is ./eGADA.
Usage:
./eGADA -i INPUTFNAME -o OUTPUTFNAME [OPTIONS]

  -h [ --help ]                         produce help message
  -T [ --TBackElim ] arg (=5)            is the backward elimination critical 
                                        value for a breakpoint. i.e. minimum 
                                        (mean1-mean2)/stddev difference between
                                        two adjacent segments.
  -a [ --aAlpha ] arg (=0.5)            is the SBL hyper prior parameter for a 
                                        breakpoint. It is the  shape parameter 
                                        of the Gamma distribution. Higher 
                                        (lower) value means less (more) 
                                        breakpoints expected a priori.
  -M [ --MinSegLen ] arg (=0)           is the minimum size in number of probes
                                        for a segment to be deemed significant.
  --BaseAmp arg (=0)                    Mean amplitude associated to the 
                                        Neutral state. If not provided, and c 
                                        option is used, then it is estimated as
                                        the median value of all probes 
                                        hybridization values after running the 
                                        algorithm. We recomend to estimate this
                                        on chromosomes that are known to have a
                                        Neutral state on most areas. In some 
                                        cases this value may be known if we 
                                        have been applied some normalization, 
                                        preprocessing or using another sample 
                                        as ref.
  -s [ --sigma2 ] arg (=-1)             Variance observed, if negative value, 
                                        it will be estimated by the mean of the
                                        differences. I would recommend to be 
                                        estimated on all the chromosomes and as
                                        a trimmed mean.
  -c [ --SelectClassifySegments ] arg (=0)
                                        Classify segment into altered state 
                                        (1), otherwise 0
  --SelectEstimateBaseAmp arg (=1)      toggle this to estimate BaseAmp from 
                                        data, rather than user-supplied.
  --convergenceDelta arg (=1e-08)       a delta number controlling convergence
  --maxNoOfIterations arg (=50000)      maximum number of iterations for EM 
                                        convergence algorithm to run before 
                                        being stopped
  --convergenceMaxAlpha arg (=100000000)
                                        one convergence related number, not 
                                        sure what it does.
  --convergenceB arg (=9.9999999999999995e-21)
                                        one convergence related number, not 
                                        sure what it does
  -b [ --debug ]                        toggle debug mode
  -r [ --report ]                       toggle report mode
  --reportIntervalDuringBE arg (=100000)
                                        how often to report the break point to 
                                        be removed during backward elimination
  -i [ --inputFname ] arg               input filename, gzipped or not. could 
                                        be specified as option or positional 
                                        argument.It is a single column text 
                                        file with no header.
  -o [ --outputFname ] arg              output filename

```

## Calling the dynamic library from Python

```bash
./src/testGADA.py -i ./data/input.txt -o ./data/output_a0.5T4M5.tsv
```
## Docker image

- https://hub.docker.com/repository/docker/polyactis/egada
- Binary and testGADA.py are deposited in /opt/eGADA.
- Run ```docker pull polyactis/egada``` or ```singularity pull docker://polyactis/egada```.
# References

1. Huang YS. eGADA: enhanced Genomic Alteration Detection Algorithm. bioRxiv. 2023
2. Tipping ME. [Sparse Bayesian learning and the relevance vector machine](https://www.jmlr.org/papers/volume1/tipping01a/tipping01a.pdf). Journal of machine learning research. 2001;1(Jun):211-44.
3. Pique-Regi R, Tsau ES, Ortega A, Seeger R, Asgharzadeh S. [Wavelet footprints and sparse bayesian learning for DNA copy number change analysis](https://ieeexplore.ieee.org/abstract/document/4217089). In 2007 IEEE International Conference on Acoustics, Speech and Signal Processing-ICASSP'07 2007 Apr 15 (Vol. 1, pp. I-353). IEEE.
4. Pique-Regi R, Monso-Varona J, Ortega A, Seeger RC, Triche TJ, Asgharzadeh S. [Sparse representation and Bayesian detection of genome copy number alterations from microarray data](https://academic.oup.com/bioinformatics/article/24/3/309/253648). https://github.com/rpique/GADA. Bioinformatics. 2008 Feb 1;24(3):309-18.
5. https://github.com/isglobal-brge/R-GADA: R implementation of GADA.
