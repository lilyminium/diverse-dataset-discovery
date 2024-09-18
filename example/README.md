# Example usage

This directory contains an example run of the script. To write out the top 100 molecules by number of matches to groups of interest to OpenFF:

```
python ../bin/select-interesting-molecules_2.2.0_v1.py -i test.smi -o top-100.smi -n 100
```

This reads in the input `test.smi` document and writes out `top-100.smi`, which contains 100 molecules from `test.smi` sorted by number of matches to groups where OpenFF currently has low data coverage.
