# Created 17 Jan 2023
# Last Edit 14 Apr 2023
# Author: Calvin XiaoYang Hu

from base_editing import Gene_BaseEditing, find_gRNAs, process_PAM

print('Welcome')

# exon_filename
print('Input the exon filename (including path from current directory): ')
exon_filename = input()
print("exon filename: ", exon_filename)

# output_dir
output_dir = 'results/'

# run code
gene = Gene_BaseEditing(filepath=exon_filename, 
                      output_dir=output_dir
                     )
gene.find_all_guides(n=23)

print("guides found successfully")

###

# BE_type
print('Select the base editor type (ABE, CBE, etc): ')
be_type = input()
print("Base editor type: ", be_type)

# cas_type
print('Select the cas type (Sp, SpG, SpRY): ')
cas_type = input()
print("Cas type: ", cas_type)

# cas_type
print('Select the PAM (None, NGNG, etc): ')
pam = input()
print("PAM type: ", pam)

# editing_window
print('Enter the editing window: ')
window = list(map(int, input("Enter multiple values: ").split()))
print(f"Editing window {window}")

# target_codons_type

# output csv filename
print('Name of output .csv file: ')
output_filename = input()
print("filename: ", output_filename)

# run code
gene.find_all_guides(n=23)
find_gRNAs(gene_object= gene, 
           mode=        be_type, 
           cas_type=    cas_type, 
           PAM=         pam
           window=      window
          ).to_csv(output_filename)

print(f"Results have been outputed to {output_filename}")
