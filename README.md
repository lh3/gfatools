## Getting Started
```sh
git clone https://github.com/lh3/gfatools
cd gfatools && make
# Extract a subgraph
./gfatools view -l MTh4502 -r 1 test/MT.gfa > sub.gfa
# Convert GFA to segment FASTA
./gfatools gfa2fa test/MT.gfa > MT-seg.fa
# Convert rGFA to stable FASTA or BED
./gfatools gfa2fa -s test/MT.gfa > MT.fa
./gfatools gfa2bed -m test/MT.gfa > MT.bed
```

## Introduction

gfatools is a set of tools for manipulating sequence graphs in the GFA or the
[rGFA][rGFA] format. It has implemented parsing, subgraph and conversion to
FASTA/BED. More functionality may be added in future.

[rGFA]: https://github.com/lh3/gfatools/blob/master/doc/rGFA.md
