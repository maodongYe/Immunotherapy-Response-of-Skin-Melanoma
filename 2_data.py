
import mygene
import pandas as pd
import os
os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"


###############################################################################
#                                                                             #
#               change gene ID from ensembl to ENSG                           #
#                                                                             #
###############################################################################
class Preprocess_gene_ID():
    """
    There are two file in need:
       expression_matrix
       gmt:  "pathway_genes_list.txt"

    The expression_matrix file could be downloaded from TCGA, GEO website and so on.
   The gmt file generate from the step change_pathway_to_gmt_format_for_ssgsea.
    """

    def __init__(self, expression_matrix, gmt):
        """
        expression_matrix: gene expression_matrix
        """
        self.expression_matrix = expression_matrix
        self.gmt = gmt
    

    def preprocess_expression(self):
        print('preprocess gene expression files.\n')
        path_ = os.path.join("./Dataset", self.expression_matrix)
        df_expr = pd.read_table(path_, header=None)
        df_expr.columns = df_expr.loc[0,]
        df_expr= df_expr.loc[1:,]
        mg = mygene.MyGeneInfo()
        df_change = mg.querymany(df_expr.id.to_list(),
                                 scopes='symbol',
                                 fields='entrezgene',
                                 species='human',
                                 as_dataframe=True)
        def f1(x):
            return df_change.loc[x, 'entrezgene']
        df_expr['id'] = df_expr['id'].map(f1)
        df_expr= df_expr[~df_expr['id'].isnull()]
        df_expr= df_expr.rename(columns={'id': 'Entrez_Gene_Id'})

        # filter bad line
        bad_idx = []
        for idx, row in df_expr.iterrows():
            eid = row[0]
            try:
                int(eid)
            except:
                bad_idx.append(idx)
        df_expr = df_expr.drop(bad_idx).reset_index(drop=True)
        df_expr = df_expr.drop_duplicates('Entrez_Gene_Id')

        # only keep pathway relative genes
        genes = []
        with open(self.gmt) as F:
            lines = [line.strip() for line in F.readlines()]
            for line in lines:
                genes += line.split('\t')[1:]
        genes = list(set(genes))
        df_expr = df_expr[df_expr['Entrez_Gene_Id'].isin(genes)]

        # output
        df_expr.to_csv('./Dataset/expression_matrix_new_tcga.txt', index=False, sep='\t')


###############################################################################
#                                                                             #
#                              Main function                                  #
#                                                                             #
###############################################################################
if __name__ == '__main__':

    Preprocess_gene_ID(expression_matrix='tcga_expre_matrix.txt',
                   gmt=os.path.join('Pathway', 'pathway_genes_list.txt')
    ).preprocess_expression()
