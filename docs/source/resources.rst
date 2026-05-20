==========================
Runtime and memory usage
==========================

Runtime and memory usage depend on the number of variants genotyped and the number of haplotypes present in the graph.

The largest dataset that we have tested contained around 39 million variants, 462 haplotypes and around 30x read coverage. With 24 cores, ``PanGenie-index`` ran in 2:48 hours using 106 GB of RAM. ``PanGenie`` with option ``-f`` ran in 55 minutes using 24 cores (around 14 CPU hours) and used 53 GB of RAM.
