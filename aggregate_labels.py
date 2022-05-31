import pandas as pd
from Bio import SeqIO
from Methods.Static_methods import label_data,delete_multilabel_data,get_chebi_name, find_dag_chebi, find_children_chebis,dataset_to_latex,read_acc_GO_term
from astropy.table import Table
import ast
import  Constants.Constants as Constants




def delete_parent_chebis(latex_address, dict_address):
    # drop CHEBI terms with children in the map
    table = Table.read(latex_address, format='latex').to_pandas()
    table.columns = ['CHEBI', 'Substrate', 'UniProt', 'SwissProt']
    table = table[table.UniProt > 0]
    chebi_list = list(table.CHEBI)
    chebi_dag = find_dag_chebi(chebi_list[0:-1])
    rm_list = find_children_chebis(chebi_list[0:-1], chebi_dag)
    [chebi_list.pop(chebi_list.index(key)) for key in rm_list]
    file = open(dict_address, "r")
    contents = file.read()
    C2GO = ast.literal_eval(contents)
    delete_chebies= list(set(C2GO.keys()) - set(chebi_list))
    [C2GO.pop(key) for key in delete_chebies]

    return C2GO

def remove_ident(Dataset, acc2GO):
    records = list(SeqIO.parse(Dataset, "fasta"))

    # remove sequences with ident100 from acc2go
    r_acc = []
    for r in records:
        r_acc.append(str(r.id).split('|')[1])
    remove_acc = []

    kh=9
    for acc in acc2GO.keys():
        if acc not in r_acc:
            remove_acc.append(acc)
    [acc2GO.pop(key) for key in remove_acc]
    return acc2GO



#C2GO = delete_parent_chebis(Constants.latex_organic, Constants.dict_C2GO_organic_with_parent)
#with open(Constants.dict_C2GO_organic_final, 'w') as f:
#    f.write(str(C2GO))
#file = open(Constants.dict_C2GO_organic_final, "r")
#contents = file.read()
#C2GO = ast.literal_eval(contents)

file = open(Constants.dict_C2GO_inorganic_final, "r")
contents = file.read()
C2GO = ast.literal_eval(contents)

#C2GO = {**C2GO_inorganic, **C2GO_organic}

#acc2GO_swiss = read_acc_GO_term(Constants.path_acc2go_trans_swissprot)
acc2GO_uni = read_acc_GO_term(Constants.path_acc2go_trans_uniprot)
acc2GO = remove_ident(Constants.UniProt_60, acc2GO_uni)
records = list(SeqIO.parse(Constants.UniProt_60, "fasta"))

labeled = label_data(C2GO, acc2GO)
one_label = delete_multilabel_data(labeled)

df_data = pd.DataFrame(one_label.items())
df_data.columns = ['seq', 'label']
df_data = df_data.groupby(by='label').count().sort_values(by='seq', ascending=False)
names = get_chebi_name(list(df_data.index))
df_data['Substrate'] = names.values()
df_data = df_data[['Substrate','seq']]
df_data['seq'] = df_data['seq'].astype(int)
cut = df_data[df_data.seq >= 10]
df_data.columns =['Substrate','UniProt']
sum = df_data.sum(numeric_only=True)
sum.name = 'Total'
df_data = df_data.append(sum, ignore_index=False)

df_data.to_csv('Dataset/UniProt/ident-60/UniProt-60-inorganic_df.csv')


chebis = list(cut.index)
label_names = list(cut.Substrate)
label_dict = {}
i = 0

with open('Dataset/UniProt/ident-60/Label_name_list_inorganic_uni_ident60_t10','w') as f:
    for c in chebis:
        label_dict[c] = i
        f.write(c + ',' + label_names[i] + ',' + str(i) + '\n')
        i += 1

final_data_label = {}
seq_label = {}
for acc in one_label.keys():
    if one_label[acc] in chebis:
        final_data_label[acc] = one_label[acc]
for r in records:
    if str(r.id).split('|')[1] in final_data_label.keys():
        chebi_seq = final_data_label[str(r.id).split('|')[1]]
        label = label_dict[chebi_seq]
        seq_label[str(r.seq)] = label
with open(Constants.data_label_uni_60_inorganic, 'w') as f:
    for key, val in seq_label.items():
        f.write(key)
        f.write(',')
        f.write(str(val)+'\n')

a=1