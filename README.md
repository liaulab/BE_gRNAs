# BE_gRNAs
# Created 17 Jan 2023
# Author: Calvin XiaoYang Hu


You can run this program from terminal or through the jupyter notebook. 
This program is used to generate all sgRNAs for CBE and ABE base editor systems for Sp, SpG, and SpRY cas variants. The editing window is also inputed by the user. 
The code has options to generate all sgRNAs, all sgRNAs which would generate a stop codon, or let the use input which codons they want to generate with base editing (last option has not been exhasutively tested). 


The code takes in a .fasta file with exons of a gene, and a 20 base pair intron segment on either side of the exon. 


1. 
To run from jupyter notebook, open the .ipynb file and input all the appropriate fields, following the example in the last cell. 

2. 
To run from terminal, cd into the BE_gRNAs directory and run
"python main.py" 
and follow the prompts. 


The output is by default saved to the results directory. 


Future directions for Calvin: 
Add option to input which codons to target, instead of which codons to generate. 
Look into how to score the sgRNAs with API or just write the code. (CRISPOR, BE Hive)
