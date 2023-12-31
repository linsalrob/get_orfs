[![Edwards Lab](https://img.shields.io/badge/Bioinformatics-EdwardsLab-03A9F4)](https://edwards.flinders.edu.au)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![GitHub language count](https://img.shields.io/github/languages/count/linsalrob/get_orfs)




# get_orfs

C code to translate a DNA sequence using different translation tables. Designed to be fast and lightweight, with few dependencies (only kseq.h, zlib and pthreads).

This is a fast multi-threaded DNA sequence translation application that will convert any DNA sequence to protein sequence.

# Should I use this code?

If you want a light-weight but powerful DNA -&gt; amino-acid translation, then perhaps. You probably want to check out the [Python](https://github.com/linsalrob/genetic_codes) library that includes this C code for an easier interface to access the data, but if you want to process things in snakemake or bash, this may be your easiest bet.


# Building from source

At the moment, the only option I provide is to build this from source, because there is the above-mentioned Python package. If there is interest I might package it using conda.

```
mkdir build
cd build/
cmake ../src/
make
cmake --install . --prefix ..
cd ..
bin/get_orfs 
```

# Using get_orfs

There are only a few options

   - `-f` or `--fasta` (required) the `fasta` format DNA sequence file to translate
   - `-t` or `--translation_table` the [NCBI translation table](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=tgencodes#SG1) to use in the translation
   - `-l` or `--length` the minimum length of ORFs to be printed. The default is to print anything one codon or longer, even if that is a stop codon (so we terminate with \*).

> Note:
> An polypeptide sequence that terminates with a stop codon includes a \* at the end (regardless of the stop codon). If the polypeptide terminates at the end of the sequence, it does not include that termination character (unless the last codon happens to be a stop codon in your translation table!).




 
