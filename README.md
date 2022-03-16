[[_TOC_]]

# eGADA

eGADA is an enhanced version of GADA (Pique-Regi 2007).

It is a fast segmentation algorithm utilizing the Sparse Bayesian Learning (or Relevance Vector Machine) technique from Tipping 2001. It can be applied to array intensity data, sequencing coverage data, or any sequential data that displays characteristics of step-wise functions.

eGADA enhancements vs GADA:

- Use a customized Red-Black tree to expedite the final backward elimination step.
- Code in C++, not C.
- Use Boost libraries extensively.
- More friendly help and commandline-argument processing.
- More friendly input and output formats.
- A dynamic library eGADA.so (packaged via Boost.Python) that can be imported in Python.
- Minor changes/bugfixes.

# Prerequisites to build eGADA

- libboost-program-options-dev: Boost program-options dev library
- libboost-iostreams-dev: Boost iostreams dev library
- libboost-python-dev: Boost python dev library


# To build the eGADA using g++


```bash
cd src;

make
```


# To run

## The C++ binary

```bash
./src/eGADA -a 0.8 -T 5 -M 3 -s -0.2 -b 0.0 -c < ../test/input.txt > ../test/output2.txt  
```
- The input is a single column plain text file.
- The output is a several column plain file with several segments.
- Type ```eGADA -h``` for more help.

## Calling the dynamic library from Python

```bash
./src/testGADA.py -i ./data/input.txt -o ./data/output_a0.5T4M5.tsv
```

# References

1. Huang YS. eGADA: an enhanced Genomic Alteration Detection Algorithm. 2022 (in preparation)
1. Tipping ME. [Sparse Bayesian learning and the relevance vector machine](https://www.jmlr.org/papers/volume1/tipping01a/tipping01a.pdf). Journal of machine learning research. 2001;1(Jun):211-44.
1. Pique-Regi R, Tsau ES, Ortega A, Seeger R, Asgharzadeh S. [Wavelet footprints and sparse bayesian learning for DNA copy number change analysis](https://ieeexplore.ieee.org/abstract/document/4217089). In 2007 IEEE International Conference on Acoustics, Speech and Signal Processing-ICASSP'07 2007 Apr 15 (Vol. 1, pp. I-353). IEEE.
1. Pique-Regi R, Monso-Varona J, Ortega A, Seeger RC, Triche TJ, Asgharzadeh S. [Sparse representation and Bayesian detection of genome copy number alterations from microarray data](https://academic.oup.com/bioinformatics/article/24/3/309/253648). Bioinformatics. 2008 Feb 1;24(3):309-18.
