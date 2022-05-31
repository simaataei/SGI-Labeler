import  pandas as pd
from astropy.table import Table
import  ast




df_uni_100 = pd.read_csv ('./Dataset/UniProt/ident-100/UniProt-100-filtered-inorganic_df.csv')
df_uni_100.columns = ['CHEBI','Substrate','UniProt-100']
df_swiss_100 = pd.read_csv ('./Dataset/SwissProt/ident-100/SwissProt-100-filtered-inorganic_df.csv')
df_swiss_100.columns = ['CHEBI','Substrate','SwissProt-100']
df_100 = pd.merge(df_swiss_100,df_uni_100, how='outer')
df_100.columns = ['CHEBI','Substrate','SwissProt-100','UniProt-100']

df_uni_60 = pd.read_csv ('./Dataset/UniProt/ident-60/UniProt-60-filtered-inorganic_df.csv')
df_uni_60.columns = ['CHEBI','Substrate','UniProt-60']
df_swiss_60 = pd.read_csv ('./Dataset/SwissProt/ident-60/SwissProt-60-filtered-inorganic_df.csv')
df_swiss_60.columns = ['CHEBI','Substrate','SwissProt-60']
df_60 = pd.merge(df_swiss_60,df_uni_60, how='outer')
df_60.columns = ['CHEBI','Substrate','SwissProt-60','UniProt-60']

df_merge = pd.merge(df_60, df_100, how='inner', on='CHEBI')
df_merge = df_merge.drop(columns='Substrate_y')
df_merge.columns = ['CHEBI','Substrate','SwissProt-60','UniProt-60','SwissProt-100','UniProt-100']
total =df_merge.loc[df_merge.CHEBI =='Total']
df_merge = df_merge.drop(df_merge[df_merge.CHEBI =='Total'].index)
df_merge = df_merge.sort_values(by='UniProt-100',ascending=False)
df_merge = df_merge.append(total, ignore_index=True)
with open('Dataset/Latex/inorganic_datasets_filtered_latex.txt', 'w') as f:
    with pd.option_context("max_colwidth", 1000):
        f.write(df_merge.to_latex(index=True))
