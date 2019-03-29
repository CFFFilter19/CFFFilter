
We design and implement a novel cuckoo filter called cuckoo filter with a fingerprint family (CFFF). The CFFF is a generalized variant of the CF that generalizes to more than two hash functions for improving the space efficiency while retaining high performance. The key to the CFFF is to use a family of fingerprints for each item, instead of a single fingerprint used by the CF. We design the addition and subtraction (ADD/SUB) operations for the CFFF, instead of the XOR operation used in the CF.


For comparsion, we implement the counting quotient filter (CQF) based on the paper recently published in SIGMOD 2017. The key to the CQF is to combine both rank-and-select-based metadata and counter embedding for improving the time and space efficiency of counting Bloom filters (CBFs).


One can use a IDE software tool called Code::Blocks to edit, compile, debug, and run the source codes of the CFFF and CQF. There are two categories of experiments for performance evaluation. In the first experiment, each filter is tested for applications that require a maximum number of items stored in the filter. In the second experiment, each filter is tested for applications that require a fixed-size memory space used for the filter.


Compile:

  g++ ./*.cpp -mbmi -mbmi2 -o CFFF

or

  g++ ./*.cpp -mbmi -mbmi2 -o CQF


Run:

  ./CFFF

or

  ./CQF


