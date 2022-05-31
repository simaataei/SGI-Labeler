from Constants import Constants
from goatools.obo_parser import GODag
from goatools.rpt.rpt_lev_depth import RptLevDepth
import obonet
import networkx
import fastobo
from owlready2 import *
import csv
import pandas as pd
import numpy as np



def read_acc_GO_term(dataset_name):
    '''

    :param dataset_name:  SwissProt or UniProt
    :return: Dictionary of sequences and their associated GO terms
    '''
    acc_GO_dict = {}
    f = open(dataset_name, 'r')
    next(f)
    all_seq = f.readlines()
    for item in all_seq:
        item = item.strip('\n')
        item = item.split('\t')
        GO_terms = item[-1].split(';')
        acc_GO_dict[item[0]] = [G.strip(' ') for G in GO_terms]
    return acc_GO_dict

def get_chebi_name(chebi_list):
    '''

    :param chebi_list: a list of chebi terms
    :return: a dictionary of chebi terms and their corresponding names
    '''
    graph = obonet.read_obo('./Ref_files/chebi_lite.obo')
    names = {}
    for chebi in chebi_list:
        names[chebi] = graph.nodes[chebi]['name']
    return names



def get_go_name():
    '''
    :return: a dictionary of GO ids and their name
    '''
    graph = obonet.read_obo("http://geneontology.org/ontology/go-basic.obo")
    id_to_name = {id_: data.get('name') for id_, data in graph.nodes(data=True)}
    return id_to_name

def get_GO_terms(seq_accession, dataset_name):
    '''

    :param seq_accession: accession number of the sequence in the dataset
    :param dataset_name: SwissProt or UniProt
    :return: GO terms associated with the sequence in the dataset
    '''
    acc_GO_dict = read_acc_GO_term(dataset_name)
    return acc_GO_dict[seq_accession]



def read_dag(GO_term):
    '''

    :param GO_term: A GO term
    :return: Directed Acyclic Graph of the GO_term in go-basic.obo file
    '''
    obodag = GODag("Ref_files/go-basic.obo")
    rptobj = RptLevDepth(obodag)
    DAG = rptobj.obo[GO_term]
    return DAG



def dag_to_list(node, dag_list):
    '''

    :param node: Directed Acyclic Graph of the GO_term in go-basic.obo file
    :return: list of GO terms in the DAG with Molucular funtion tag
    '''
    if node is not None:
        if node.namespace == 'molecular_function':
                dag_list.append(node.id)
        for n in node.children:
            if n.namespace == 'molecular_function':
                dag_to_list(n, dag_list)
    return dag_list


def find_dag_chebi(chebi_list):
    chebi_dag = {}
    for che in chebi_list:
        chebi_dag[che] = []
    knowledge_graph = networkx.DiGraph()
    #read obo file
    pato = fastobo.load('./Ref_files/chebi_lite.obo')

    # populate the knowledge graph with is_a relationships
    for frame in pato:
        if isinstance(frame, fastobo.term.TermFrame):
            knowledge_graph.add_node(str(frame.id))
            for clause in frame:
                if isinstance(clause, fastobo.term.IsAClause):
                    knowledge_graph.add_edge(str(frame.id), str(clause.term))

    # find the leaves from all of the nodes
    nodes = knowledge_graph.nodes

    for chebi in chebi_list:
      chebi_dag[chebi] = list(networkx.ancestors(knowledge_graph, chebi))


    #with open('Mid_files/dict_CHEBI_Dag_Ontoclass.txt', 'w') as data:
     #   data.write(str(chebi_dag))

    return chebi_dag
def find_children_chebis(chebi_list, chebi_dag):
    '''

    :param chebi_list: a list of chebi terms
    :param chebi_dag: a dictionary of chebi terms and their childeren in DAG
    :return: a list of chebi terms that have children in the list
    '''
    rm_list = []
    for chebi in chebi_list:
        other_chebis = list(set(chebi_list)- set([chebi]))
        dag_childeren_list = chebi_dag[chebi]
        if set(other_chebis).intersection(dag_childeren_list):
            rm_list.append(chebi)
    return rm_list


def find_GO_from_chebi(node, chebi_list):
    chebi2GO = {}
    if node != None:
        if len(node.INDIRECT_equivalent_to):
            is_a_list = list(node.INDIRECT_equivalent_to[0].Classes)
            for item in is_a_list:
                if type(item) == Restriction:
                    if 'obo.CHEBI' in str(item.value):
                            if item.value.name.replace('_', ':') in chebi_list:
                                chebi2GO[item.value.name.replace('_', ':')] =node.name
    return chebi2GO


def map_chebi2GO(chebi_list):
    C2GO = {}
    for item in chebi_list:
        C2GO[item] = []
    go = get_ontology("./Ref_files/go-plus.owl").load()
    obo = go.get_namespace("http://purl.obolibrary.org/obo/")

    dag_TA_list = list(set(dag_to_list(read_dag('GO:0005215'), [])))
    #dag_catalytic_list = list(set(dag_to_list(read_dag('GO:0003824'), [])))
   # dag_TA_NC = [item for item in dag_TA_list if item not in dag_catalytic_list]

    for child in dag_TA_list:
        child = child.replace(':','_')
        child_go = find_GO_from_chebi(obo[child], chebi_list)
        for chebi in child_go.keys():
            C2GO[chebi].append(child_go[chebi].replace('_',':'))
    rm_list = []
    for item in C2GO.keys():
        if C2GO[item] == []:
            rm_list.append(item)
    [C2GO.pop(key) for key in rm_list]

    return C2GO

def label_data(C2GO, acc2GO):
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
                data_label_dict[seq].append(item)

    return data_label_dict

def delete_multilabel_data(data_label_dict):
    '''

    :param data_label_dict: labeled dictionary of accession with a list of lables
    :return: data_label_dict: a dictionary of accessions with one label on them
    '''
    delete_double_labels = []
    for key in data_label_dict.keys():
        if len(data_label_dict[key]) != 1:
            delete_double_labels.append(key)

    [data_label_dict.pop(k) for k in delete_double_labels]
    for item in data_label_dict.keys():
        data_label_dict[item] = data_label_dict[item][0]
    return data_label_dict
def dataset_to_latex(C2GO, acc2GO_swiss, acc2GO_uni,method):
    data_label_dict_uni = {}
    data_label_dict_swiss = {}
    label_list = {}
    i = 0
    for c in C2GO:
        label_list[c] = i
        i += 1
    for acc in acc2GO_uni.keys():
        data_label_dict_uni[acc] = []
    for acc in acc2GO_swiss.keys():
        data_label_dict_swiss[acc] = []
    C2GO_df = pd.DataFrame([C2GO]).T
    C2GO_df.columns = ['GO_term']
    with open('C2GO_latex.txt', 'w') as f:
        with pd.option_context("max_colwidth", 1000):
            f.write(C2GO_df.to_latex(index=True))

    names = get_chebi_name(C2GO_df.index).values()
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
    C2GO_df = C2GO_df.astype({'UniProt': int})
    C2GO_df = C2GO_df.astype({'SwissProt': int})
    sum = C2GO_df.sum(numeric_only=True)
    sum.name = 'Total'
    C2GO_df = C2GO_df.append(sum, ignore_index=False)
    C2GO_df = C2GO_df.replace(np.nan, '-')
    df1 = C2GO_df[['GO_term']]
    df2 = C2GO_df[['Substrate', 'UniProt', 'SwissProt']]
    with open('Mid_files/CHEBI2GO_'+method+'.txt', 'w') as f:
        with pd.option_context("max_colwidth", 1000):
            f.write(df1.to_latex(index=True))
    with open('Dataset/Latex/CHEBI_from_'+method+'_dataset_latex.txt', 'w') as f:
        with pd.option_context("max_colwidth", 1000):
            f.write(df2.to_latex(index=True))

    with open('./Dataset/UniProt/Data_label_UniProt_'+method+'.txt', 'w') as f:
        for key in data_label_dict_uni.keys():
            f.write("%s\t%s\n" % (key, data_label_dict_uni[key]))
    with open('./Dataset/SwissProt/Data_label_SwissProt_'+method+'.txt', 'w') as f:
        for key in data_label_dict_swiss.keys():
            f.write("%s\t%s\n" % (key, data_label_dict_swiss[key]))
    with open('./Dataset/label_list_' + method + '.txt', 'w') as f:
        for key in label_list.keys():
            f.write("%s\t%s\n" % (C2GO_df.at[key,'Substrate'], label_list[key]))

