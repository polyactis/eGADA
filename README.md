[[_TOC_]]

# eGADA
eGADA is an enhanced version of GADA, https://github.com/rpique/GADA.

- Use a Red-Black tree to expedite the final backward elimination step.
- Code in C++, not C.
- Use Boost libraries extensively.
- More friendly help and commandline-argument processing.
- More friendly input and output formats.
- A dynamic library eGADA.so (packaged via Boost.Python) that can be imported by Python2.

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