# Created 17 Jan 2023
# Author: Calvin XaoYang Hu

from simple_term_menu import TerminalMenu
from BE_gRNAs import BE_gRNAs

print('Welcome')

# gene_name
print('Name of gene: ')
gene_name = input()
print("Notes: ", gene_name)

# BE_type
BE_options = ["ABE", "CBE"]
BE_terminal_menu = TerminalMenu(BE_options, title='Select the base editor type: ')
BE_choice_index = BE_terminal_menu.show()
BE_type = BE_options[BE_choice_index]
print(f"You have selected {BE_type}")

# cas_type
cas_options = ["Sp (NGG)", "SpG (NGN)", "SpRY (NNN)"]
cas_terminal_menu = TerminalMenu(cas_options, title='Select the cas type: ')
cas_choice_index = cas_terminal_menu.show()
cas_type = cas_options[cas_choice_index].split()[0]
print(f"You have selected {cas_type}")

# editing_window
print('Enter the editing window: ')
editing_window = list(map(int, input("Enter multiple values: ").split()))
print(f"Editing window {editing_window}")

# target_codons_type
target_codons_options = ['[] no target', '[TAG, TAA, TGA] stop codons', 'Other']
target_codons_terminal_menu = TerminalMenu(target_codons_options, title='Select the target codons to mutate into: ')
target_codons_choice_index = target_codons_terminal_menu.show()
target_codons_type = target_codons_options[target_codons_choice_index]
if target_codons_type == '[] no target':
    target_codons = []
    print(f"You have selected {cas_type}")
elif target_codons_type == '[TAG, TAA, TGA] stop codons':
    print(f"You have selected {cas_type}")
    target_codons = ['TAG', 'TAA', 'TGA']
else:
    target_codons = list(map(string, input("Enter the target codons: ").split()))

# exon_filename
print('Input the exon filename (including path from current directory): ')
exon_filename = input()
print("exon filename: ", exon_filename)

# notes
print('Notes for file naming: ')
notes = input()
print("Notes: ", notes)

# output_dir
output_dir = 'results/'


# run code
gene_BE = BE_gRNAs(BE_type, editing_window, gene_name, target_codons, cas_type, exon_filename)
#


output_path = output_dir + gene_name + '_' + cas_type + BE_type + '_' + notes + '_gRNAs.csv'
gene_BE.save_data(output_path)

print(f"Results have been outputed to {output_path}")
