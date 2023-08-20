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
./src/eGADA -a 0.8 -T 5 -M 3 -s -0.2 -b 0.0 -c < ../test/input.txt > ../test/output2.txt  
```
- The input is a single column plain text file.
- The output is a several column plain file with several segments.
- Type ```eGADA -h``` for more help.

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
