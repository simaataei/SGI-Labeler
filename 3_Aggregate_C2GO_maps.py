import ast
from Methods.Static_methods import read_acc_GO_term,get_chebi_name
import pandas as pd
def intersect(a, b):
    sb = set(b)
    return [x for x in a if x in sb]
file = open("Mid_files/dict_C2GO_from_CHEBI_with_DAG_organic.txt", "r")
contents = file.read()
C2GO_from_CHEBI = ast.literal_eval(contents)

file = open("Mid_files/dict_C2GO_from_GO_organic_leaf_CHEBI_filtered.txt", "r")
contents = file.read()
C2GO_from_GO = ast.literal_eval(contents)
keys = list(set(list(C2GO_from_CHEBI.keys())+ list(C2GO_from_GO.keys())))
inter = set(C2GO_from_GO.keys()).intersection(C2GO_from_CHEBI.keys())
C2GO = {}
for k in keys:
    C2GO[k] =[]
    if k in C2GO_from_CHEBI.keys():
        for item in C2GO_from_CHEBI[k]:
            C2GO[k].append(item)
    if k in C2GO_from_GO.keys():
        for item in C2GO_from_GO[k]:
            C2GO[k].append(item)
    C2GO[k] = list(set(C2GO[k]))
with open('Mid_files/dict_C2GO_overall_CHEBI_filtered_organic.txt', 'w') as d:
    d.write(str(C2GO))

'''
acc2GO_swiss = read_acc_GO_term('SwissProt')
acc2GO_uni = read_acc_GO_term('UniProt')

data_label_dict_uni = {}
data_label_dict_swiss = {}
for acc in acc2GO_swiss.keys():
    data_label_dict_swiss[acc] = []
for acc in acc2GO_uni.keys():
    data_label_dict_uni[acc] = []
label_list = {}
i = 0
for c in C2GO:
    label_list[c] = i
    i += 1

C2GO_df = pd.DataFrame([C2GO]).T
C2GO_df.columns = ['GO_term']
with open('C2GO_latex.txt', 'w') as f:
    with pd.option_context("max_colwidth", 1000):
        f.write(C2GO_df.to_latex(index=True))

names = get_chebi_name(C2GO_df.index)
C2GO_df.insert(loc=1, column='Substrate', value=names)
C2GO_df.insert(loc=2, column='SwissProt', value=0)
C2GO_df.insert(loc=3, column='UniProt', value=0)

for item in C2GO.keys():
    for seq in acc2GO_uni.keys():
        if intersect(C2GO[item], acc2GO_uni[seq]):
            #   C2GO_df.at[item,'Accession_list'].append(seq)
            C2GO_df.at[item, 'UniProt'] += 1
            data_label_dict_uni[seq].append(label_list[item])

    for seq in acc2GO_swiss.keys():
        if intersect(C2GO[item], acc2GO_swiss[seq]):
            #   C2GO_df.at[item,'Accession_list'].append(seq)
            C2GO_df.at[item, 'SwissProt'] += 1
            data_label_dict_swiss[seq].append(label_list[item])

C2GO_df = C2GO_df.sort_values(by='UniProt', ascending=False)
C2GO_df.astype({'UniProt': int})
C2GO_df.astype({'SwissProt': int})
df_sum = {'GO_term': '-', 'Substrate': '-', 'SwissProt': sum(C2GO_df.SwissProt), 'UniProt': sum(C2GO_df.UniProt)}
C2GO_df.append(df_sum, ignore_index=True)
df1 = C2GO_df[['GO_term']]
df2 = C2GO_df[['Substrate', 'UniProt', 'SwissProt']]
with open('Mid_files/CHEBI2GO_trans_all_latex.txt', 'w') as f:
    with pd.option_context("max_colwidth", 1000):
        f.write(df1.to_latex(index=True))
with open('Dataset/Latex/CHEBI_from_trans_all.txt', 'w') as f:
    with pd.option_context("max_colwidth", 1000):
        f.write(df2.to_latex(index=True))

with open('./Dataset/UniProt/Data_label_UniProt_trans_all.csv', 'w') as f:
    for key in data_label_dict_uni.keys():
        f.write("%s,%s\n" % (key, data_label_dict_uni[key]))
with open('./Dataset/SwissProt/Data_label_SwissProt_trans_all.csv', 'w') as f:
    for key in data_label_dict_swiss.keys():
        f.write("%s,%s\n" % (key, data_label_dict_swiss[key]))
with open('./Dataset/label_list_trans_all.csv.csv', 'w') as f:
    for key in label_list.keys():
        f.write("%s,%s\n" % (C2GO_df.at[key, 'Substrate'], label_list[key]))

a=1
'''