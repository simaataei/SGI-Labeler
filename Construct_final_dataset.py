import pandas as pd
from astropy.table import Table
import ast
from Methods.Static_methods import dataset_to_latex,get_chebi_name,read_acc_GO_term
from Bio import SeqIO
import  Constants.Constants as Constant
def C2G_build(Dataset):
    table = Table.read('Dataset/Latex/Modified_datasets_inorganic.txt', format='latex').to_pandas()
    table.columns = ['Substrate', 'UniProt', 'SwissProt']
    if 'SwissProt' in Dataset:
        cut = table.loc[table.SwissProt >= 8]
    elif 'UniProt' in Dataset:
        cut = table.loc[table.UniProt >= 8]
    name_uni = cut.Substrate
    file = open("Mid_files/dict_C2GO_overall_without_double_and_parent_chebi_inorganic.txt", "r")
    contents = file.read()
    C2GO = ast.literal_eval(contents)
    chebi2name = get_chebi_name(C2GO.keys())


    c2name_dict = {}
    for name in list(name_uni):
        for key, val in chebi2name.items():
            if name == val:
                c2name_dict[key] = name
    for key, val in chebi2name.items():
        if 'ATP' in val:
            c2name_dict[key] = 'ATP(4−)/ADP(-3)'
            break
    delete_list = [key for key in C2GO.keys() if key not in c2name_dict.keys()]
    [C2GO.pop(k) for k in delete_list]
    with open('Mid_files/dict_C2GO_final_'+Dataset+'_inorganic.txt', 'w') as f :
        f.write(str(C2GO))

    C2GO_df = pd.DataFrame([C2GO]).T
    C2GO_df.columns = ['GO_term']
    with open('Dataset/Latex/C2GO_final_' + Dataset + '_latex.txt', 'w') as f:
        with pd.option_context("max_colwidth", 1000):
            f.write(C2GO_df.to_latex(index=True))
    return C2GO

def build_data_label_file(Dataset,acc2GO,C2GO,names):
    data_label_dict = {}
    for acc in acc2GO:
        data_label_dict[acc] = []
    label_list = {}
    i = 0
    for c in C2GO:
        label_list[c] = i
        i += 1
    for item in C2GO.keys():
        for seq in acc2GO.keys():
            if set(C2GO[item]).intersection(acc2GO[seq]):
                data_label_dict[seq].append(label_list[item])

    delete_double_labels = []
    for key in data_label_dict.keys():
        if len(data_label_dict[key]) > 1:
            delete_double_labels.append(key)

    [data_label_dict.pop(k) for k in delete_double_labels]
    records = list(SeqIO.parse(Dataset, "fasta"))
    string_label ={}
    for r in records:
        if str(r.id).split('|')[1] in data_label_dict.keys():
            label = str(data_label_dict[str(r.id).split('|')[1]]).strip('[').strip(']')
            if len(label) > 0:
                string_label[str(r.seq)] = int(label)
    table = pd.DataFrame(string_label.items())
    table.columns = ['seq', 'label']
    table = table.groupby(by='label').count().sort_values(by='seq', ascending=False)

    #find chebi names and ids for each label
    chebi_num = []
    chebi_name = []
    for i in list(table.index):
        for k, val in label_list.items():
            if i == val:
                chebi_num.append(k)
                chebi_name.append(names[k])
    table['CHEBI'] = chebi_num
    table['Substrate'] = chebi_name
    sum = table.sum(numeric_only=True)
    sum.name = 'Total'
    table = table.append(sum, ignore_index=False)
    table = table[['CHEBI', 'Substrate', 'seq']]
    table = table.astype({'seq': int})
    table = table.drop(21)
    with open('temp_table_60_swiss.txt', 'w') as f:
        with pd.option_context("max_colwidth", 1000):
            f.write(table.to_latex(index=True))
    table2 = Table.read('temp_table_60_uni.txt', format='latex').to_pandas()
    table2 = table2.astype({'seq': float})
    table2 = table2.drop(23)
    table2 = table2.drop(0)
    papr = pd.merge(table, table2, how='inner', on='CHEBI')
    proper_table = papr[['CHEBI', 'Substrate_x', 'seq_x', 'seq_y']]
    proper_table.columns = ['CHEBI', 'Substrate', 'SwissProt', 'UniProt']
    sum = proper_table.sum(numeric_only=True)
    sum.name = 'Total'
    proper_table = proper_table.append(sum, ignore_index=False)
    with open('Dataset/Latex/Uni_swiss_paper_latex.txt', 'w') as f:
        with pd.option_context("max_colwidth", 1000):
            f.write(proper_table.to_latex(index=True))

    if 'UniProt' in Dataset:
        with open(Constant.path_uni_ident_60+'UniProt_modified_data_one_label_t10_ident60.txt','w') as f:
            for k,v in string_label.items():
                f.write("%s,%s\n" % (k,v))
        with open(Constant.path_uni_ident_60 + 'Data_label_list_uni_ident60_t10', 'w') as f:
            for k,v, n in zip(list(label_list.keys()), list(label_list.values()), list(names.values())):
                f.write("%s,%s,%s\n" % (k,n,v))
    elif 'SwissProt' in Dataset:
        with open(Constant.path_swiss_ident_60 + 'SwissProt_modified_data_one_label_t10_ident60.txt', 'w') as f:
            for k, v in string_label.items():
                f.write("%s,%s\n" % (k, v))
        with open(Constant.path_swiss_ident_60 + 'Data_label_list_swiss_ident60_t10', 'w') as f:
            for k, v, n in zip(list(label_list.keys()), list(label_list.values()), names):
                f.write("%s,%s,%s\n" % (k, n, v))
    return table

C2GO = C2G_build('UniProt')
names = get_chebi_name(C2GO.keys())
acc2GO_uni = read_acc_GO_term(Constant.path_acc2go_trans_uniprot)
acc2GO_swiss = read_acc_GO_term(Constant.path_acc2go_trans_swissprot)
#dataset_to_latex(C2GO,acc2GO_swiss,acc2GO_uni,'Final')
uni_table = build_data_label_file(Constant.UniProt,acc2GO_uni,C2GO, names)



a=1