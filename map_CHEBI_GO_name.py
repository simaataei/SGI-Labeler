from Methods.Static_methods import get_chebi_name, get_go_name
import ast
import pandas as pd


file = open("Mid_files/dict_C2GO_overall_without_double_and_parent_chebi_organic.txt", "r")
contents = file.read()
C2GO = ast.literal_eval(contents)
C2GO_name = {}

GO_name = get_go_name()
chebi_name = get_chebi_name(C2GO.keys())


for C in C2GO.keys():
    C2GO_name[C+'|'+chebi_name[C]] = []
    for GO in C2GO[C]:
        C2GO_name[C+'|'+chebi_name[C]].append(GO+'|'+GO_name[GO])

C2GO_df = pd.DataFrame([C2GO_name]).T
C2GO_df.columns = ['GO_terms']

with open('Mid_files/CHEBI2GO_names_organic.txt', 'w') as f:
    with pd.option_context("max_colwidth", 5000):
        f.write(C2GO_df.to_latex(index=True, longtable=True))