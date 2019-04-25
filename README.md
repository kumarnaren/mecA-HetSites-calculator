# mecA-HetSites-calculator
calculates mecA coverage, coverage depth and reference genome coverage and number heterozygous sites

Essential requirements
1. python2.7
2. samtools version 1.3.1 or higher
3. bcftools version 1.3.1 or higher
4. a bed file with coordinates of mobile genetic elements
5. Reference fasta file

Please edit the path of the respective files in the mecAndHetsites.py at lines 9, 10, 11, 12 and 82.

Example:
mecAndHetsites.py <folder> <outputName>
