import numpy as np
from sklearn.model_selection import train_test_split
import re
import Constants.Constants as Constants
from astropy.table import Table

def test_train(Dataset):
    X = []
    y = []
    with open(Dataset) as f:
        data = f.readlines()
        for d in data:
            d = d.split(',')
            X.append(d[0])
            y.append(d[1].strip('\n'))

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42, stratify=y)
    X_train, X_val, y_train, y_val = train_test_split(X_train, y_train, test_size=0.2, random_state=42,stratify=y_train)
    X_train = [re.sub(r'(\w)', r'\1 ', x) for x in X_train]
    X_train = [re.sub(r"[UZOB]", "X", x) for x in X_train]

    X_test = [re.sub(r'(\w)', r'\1 ', x) for x in X_test]
    X_test = [re.sub(r"[UZOB]", "X", x) for x in X_test]

    X_val = [re.sub(r'(\w)', r'\1 ', x) for x in X_val]
    X_val = [re.sub(r"[UZOB]", "X", x) for x in X_val]

    return X_train,y_train,X_val,y_val,X_test,y_test


X_train, y_train, X_val, y_val, X_test, y_test = test_train('Dataset/UniProt/ident-60/data_label_inorganic_uniProt_ident60.txt')
with open('./Dataset/UniProt/ident-60/inorganic_uniprot_ident60_train.csv', 'w') as f:
    f.write("%s,%s\n" % ('Sequence', 'label'))
    for k, v in zip(X_train, y_train):
        f.write("%s,%s\n" % (k, v))
with open('./Dataset/UniProt/ident-60/inorganic_uniprot_ident60_test.csv', 'w') as f:
    f.write("%s,%s\n" % ('Sequence', 'label'))
    for k, v in zip(X_test, y_test):
        f.write("%s,%s\n" % (k, v))
with open('./Dataset/UniProt/ident-60/inorganic_uniprot_ident60_validation.csv', 'w') as f:
    f.write("%s,%s\n" % ('Sequence', 'label'))
    for k, v in zip(X_val, y_val):
        f.write("%s,%s\n" % (k, v))
a=1