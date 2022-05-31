import pandas as pd
from Methods.Static_methods import read_dag,find_dag_chebi
import pandas as pd
from owlready2 import *
import ast



def read_leaf_list(dataset_name):
    data = pd.read_csv('./Dataset/'+dataset_name+'/'+dataset_name+'_leaf.txt', sep=",")
    return set(data.GO_term)

def get_leaf_nodes(node,leafs):
    if node is not None:
        if len(node.children) == 0:
            if node.namespace =='molecular_function':
                leafs.append(node)

        for n in node.children:
            if n.namespace == 'molecular_function':
                get_leaf_nodes(n, leafs)
    return leafs
def get_chebi_exact(node):
    '''gets chebi term only from the go term itself according to the go.owl file'''
    chebi = []
    if node != None:
        if len(node.INDIRECT_equivalent_to):
            is_a_list =list(node.INDIRECT_equivalent_to[0].Classes)
            for item in is_a_list:
                if type(item) == Restriction:
                    if 'obo.CHEBI' in str(item.value):
                        chebi.append(item.value.name.replace('_', ':'))
    return chebi
def map_GO2chebi(GO_terms):

    go = get_ontology("./Ref_files/go-plus.owl").load()
    obo = go.get_namespace("http://purl.obolibrary.org/obo/")
    G2C = {}
    for i in GO_terms:
        G2C[i] = []
    for g in GO_terms:
        go_term = g.replace(':', '_')
        chebi = get_chebi_exact(obo[go_term])
        G2C[g] = chebi

    return G2C
transporters = read_dag('GO:0022857')
catalytic = read_dag('GO:0003824')

#inorganic_cation = read_dag('GO:0022890')
#inorganic_anion = read_dag('GO:0015103')

#organic_cation = read_dag('GO:0015101')
#organic_anion = read_dag('GO:0008514')


# Extract transporter leafs
t_leafs = get_leaf_nodes(transporters,[])
#t_leafs = get_leaf_nodes(organic_anion,[]) + get_leaf_nodes(organic_cation,[])
t_leaf_list = []
for l in list(set(t_leafs)):
    t_leaf_list.append(l.id)




# Extract catalytic leafs
c_leafs = get_leaf_nodes(catalytic,[])
c_leaf_list = []
for l in list(set(c_leafs)):
    c_leaf_list.append(l.id)




# Exclude catalytic from the transporter leafs
modified_leaf = [value for value in t_leaf_list if value not in c_leaf_list]
with open("Ref_files/GO_leaf_trans.txt", "w") as f:
    for item in modified_leaf:
        f.write(item + "\n")
GO2C = map_GO2chebi(modified_leaf)



with open('Mid_files/dict_GO2C_organic.txt', 'w') as d:
    d.write(str(GO2C))
C2GO ={}

flat_list = [item for sublist in GO2C.values() for item in sublist]
for chebi in flat_list:
    C2GO[chebi] =[]
for GO in GO2C.keys():
    for item in GO2C[GO]:
        C2GO[item].append(GO)


# Exclude chebis that are not inorganic(organic) from the map
organic_cation = ['CHEBI:25697']
organic_anion = ['CHEBI:25696']
#inorganic_cation = ['CHEBI:36915']
#inorganic_anion = ['CHEBI:24834']
cations = list(find_dag_chebi(organic_cation).values())[0]
anions = list(find_dag_chebi(organic_anion).values())[0]
remove_chebi = []
for C in C2GO.keys():
    if not C in cations and not C in anions:
        remove_chebi.append(C)

[C2GO.pop(key) for key in remove_chebi]


with open('Mid_files/dict_C2GO_from_GO_organic_leaf_CHEBI_filtered.txt', 'w') as d:
    d.write(str(C2GO))
diffrence_list = [value for value in t_leaf_list if value in c_leaf_list]
