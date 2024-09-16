# Example usage

This directory contains an example run of the script. To write out the top 100 molecules by number of matches to groups with low coverage:

```
python ../bin/sort-by-rare-groups_2.2.0_v1.py -i test.smi -o sorted_by_rarity.csv --no-write-all -n 100
```

This reads in the input `test.smi` document and writes out `sorted_by_rarity.csv`, which contains 100 molecules from `test.smi` sorted by number of matches to groups with low coverage.
