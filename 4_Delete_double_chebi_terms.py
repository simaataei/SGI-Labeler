import ast
from Methods.Static_methods import read_acc_GO_term, dataset_to_latex
from Constants.Constants import path_acc2go_trans_swissprot, path_acc2go_trans_uniprot



# read CHEBI leaf file
chebi_leaf = []
with open('Ref_files/chebi_leaf.txt') as f:
    data = f.readlines()
    for row in data:
        chebi_leaf.append(row.strip('\n'))

# read C2GO map concatinated from both GO side and CHEBI side
file = open("Mid_files/dict_C2GO_overall_CHEBI_filtered_organic.txt", "r")
contents = file.read()
C2GO = ast.literal_eval(contents)

# reverse the C2GO map to GO2C
GO2C = {}
flat_list = [item for sublist in C2GO.values() for item in sublist]
for GO in set(flat_list):
    GO2C[GO] = []

for chebi in C2GO.keys():
    for item in C2GO[chebi]:
        GO2C[item].append(chebi)

# For GO terms with more than one CHEBI term, delete non-leaf CHEBIs
drop_list = []
for GO in GO2C.keys():
    if len(GO2C[GO]) > 1:
        remove_list = []
        for c in GO2C[GO]:
            if c not in chebi_leaf:
                if c !='CHEBI:30616' and c != 'CHEBI:456216':
                    remove_list.append(c)
        for r in remove_list:
            GO2C[GO].remove(r)
        if len(GO2C[GO]) == 0:
            drop_list.append(GO)
# if a GO term lost all its CHEBIs drop it from the map
[GO2C.pop(key) for key in drop_list]


# reconstruct the C2GO map
C2GO_new = {}
chebis = [item for sublist in GO2C.values() for item in sublist]
for GO in set(chebis):
    C2GO_new[GO] = []

for GO in GO2C.keys():
    for item in GO2C[GO]:
        C2GO_new[item].append(GO)

with open('Mid_files/dict_C2GO_overall_CHEBI_filtered_without_double_chebi_organic.txt', 'w') as d:
    d.write(str(C2GO_new))
#build new dataset based on this map
acc2GO_swiss = read_acc_GO_term(path_acc2go_trans_swissprot)
acc2GO_uni = read_acc_GO_term(path_acc2go_trans_uniprot)

dataset_to_latex(C2GO_new, acc2GO_swiss,acc2GO_uni,'trans_all_CHEBI_filtered_without_double_chebi_organic')

