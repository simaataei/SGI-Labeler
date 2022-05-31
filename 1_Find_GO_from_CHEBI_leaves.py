from Methods.Static_methods import dag_to_list,read_dag, find_dag_chebi, read_acc_GO_term, map_chebi2GO, get_chebi_name,dataset_to_latex
import ast
import csv
from astropy.table import Table
import pandas as pd

# find the DAG of inorganic cation and inorganic anion
#inorganic_cation = ['CHEBI:36915']
#inorganic_anion = ['CHEBI:24834']


organic_cation = ['CHEBI:25697']
organic_anion = ['CHEBI:25696']


#root_list = [inorganic_cation, inorganic_anion]
#cations = list(find_dag_chebi(organic_cation).values())[0]
#anions = list(find_dag_chebi(organic_anion).values())[0]


# read CHEBI leaf list
with open('Ref_files/chebi_leaf.txt') as f:
    chebi_leaf = f.readlines()
    chebi_leaf = [ch.strip('\n') for ch in chebi_leaf]

# find leaf CHEBI terms for cations and anions
#leaf_cations = [c for c in cations if c in chebi_leaf]
#leaf_anions = [a for a in anions if a in chebi_leaf]
#leaf_cation_names = get_chebi_name(leaf_cations)
#leaf_anion_names = get_chebi_name(leaf_anions)



# find C2GO map for leaf cation and anion
#C2GO_cation = map_chebi2GO(leaf_cations)
#C2GO_anion = map_chebi2GO(leaf_anions)
C2GO = map_chebi2GO(chebi_leaf)

# delete chebis that are not in transmembrane transporter
transporters = read_dag('GO:0022857')

'''
with open('Mid_files/dict_C2GO_cation.txt', 'w') as d:
    d.write(str(C2GO_cation))


with open('Mid_files/dict_C2GO_anion.txt', 'w') as d:
    d.write(str(C2GO_anion))



# read C2GO map for cation and anion
file = open("Mid_files/dict_C2GO_cation.txt", "r")
contents = file.read()
C2GO_cation = ast.literal_eval(contents)

file = open("Mid_files/dict_C2GO_anion.txt", "r")
contents = file.read()
C2GO_anion = ast.literal_eval(contents)
'''


# combine two maps to one C2GO
C2GO = {**C2GO_cation, **C2GO_anion}


GO_terms = [item for sublist in list(C2GO.values()) for item in sublist]
go_dag_dict = {}
for go in GO_terms:
    go_dag_dict[go] = []
    go_dag_dict[go] = dag_to_list(read_dag(go),[])
#with open('Mid_files/GO_DAG_dict_organic.txt', 'w') as d:
  #  d.write(str(go_dag_dict))

#file = open("Mid_files/GO_DAG_dict.txt", "r")
#contents = file.read()
#go_dag_dict = ast.literal_eval(contents)
for chebi in C2GO.keys():
    gos = C2GO[chebi]
    dag_list = []
    for go in gos:
        for g in go_dag_dict[go]:
            dag_list.append(g)
    C2GO[chebi] = C2GO[chebi] + dag_list
    C2GO[chebi] = list(set(C2GO[chebi]))
a = 1
with open('Mid_files/dict_C2GO_from_CHEBI_with_DAG_organic.txt', 'w') as d:
    d.write(str(C2GO))

'''

#cation_names = get_chebi_name(C2GO_cation.keys())
#anion_names = get_chebi_name(C2GO_anion.keys())


acc2GO_swiss = read_acc_GO_term('SwissProt')
acc2GO_uni = read_acc_GO_term('UniProt')

data_label_dict_uni = {}
data_label_dict_swiss = {}
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
        if set(C2GO[item]).intersection(acc2GO_uni[seq]):
            #   C2GO_df.at[item,'Accession_list'].append(seq)
            C2GO_df.at[item, 'UniProt'] += 1
            data_label_dict_uni[seq].append(label_list[item])

    for seq in acc2GO_swiss.keys():
        if set(C2GO[item]).intersection(acc2GO_swiss[seq]):
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
with open('Mid_files/Transmembrane_transporter/CHEBI2GO_trans_Inorganic_anion_cation.txt', 'w') as f:
    with pd.option_context("max_colwidth", 1000):
        f.write(df1.to_latex(index=True))
with open('Datasets/Latex/CHEBI_from_trans_Inorganic_anion_cation.txt', 'w') as f:
    with pd.option_context("max_colwidth", 1000):
        f.write(df2.to_latex(index=True))

with open('./Datasets/Transmembrane_transport/UniProt/Data_label_UniProt_trans_Inorganic_cation_anion.csv', 'w') as f:
    for key in data_label_dict_uni.keys():
        f.write("%s,%s\n" % (key, data_label_dict_uni[key]))
with open('./Datasets/Transmembrane_transport/SwissProt/Data_label_SwissProt_trans_Inorganic_cation_anion.csv', 'w') as f:
    for key in data_label_dict_swiss.keys():
        f.write("%s,%s\n" % (key, data_label_dict_swiss[key]))
with open('./Datasets/label_list_trans_Inorganic_cation_anion.csv.csv', 'w') as f:
    for key in label_list.keys():
        f.write("%s,%s\n" % (C2GO_df.at[key, 'Substrate'], label_list[key]))
'''