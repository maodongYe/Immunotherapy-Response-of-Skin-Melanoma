import torch
import glob
import numpy as np
import pandas as pd
from torch_geometric.data import Data
import os
os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"



###############################################################################
#                                                                             #
#                      Preprocess pathway                                     #
#                                                                             #
###############################################################################
class Preprocess_pathway():
    """
    Used gene expression and pathway network to construct graph dataset.

    Args:
        (1) gene_expression_file: Gene expression matrix file with specific format.
            the format should be like this (sep=TAB):
                Entrez_Gene_Id | sample1 | sample2 | sample3 | ...

        (2) pathway_files: pathway network file
            the network file should include two columns named: "scr": sorce node; "dest": dest node

        (3) save_dataset: set graph dataset saved filename.
    """
    def __init__(self, gene_expression_file='./Dataset/expression_matrix_new_tcga.txt',
                       pathway_files='./Pathway/pathways/R*.txt',
                       save_dataset ='./Dataset/tcga_pathway_reactome.pt',
                       keep_pathway='./Pathway/pathways_keep.tsv',
                       pathway_score=None,
                       ):
        self.gene_expression_file = gene_expression_file
        self.pathway_files = pathway_files
        self.save_dataset = save_dataset
        self.pathway_score = pathway_score
        self.keep_pathway = keep_pathway

    def get_score(self, smp, pathway_name):
        """
        get pathway activate level from ssGSEA.

        smp: patient id
        pathway_name: pathway_name
        """
        if self.pathway_score != None:
            df_score = pd.read_csv(self.pathway_score)
            df_score.columns = ['Patient_ID'] + [i for i in df_score.columns[1:]]
            df_score = df_score.set_index('Patient_ID')
            df_score = df_score.apply(lambda x: (x - np.min(x)) / (np.max(x) - np.min(x)))
            try:
                return torch.tensor(round(df_score.loc[pathway_name, smp], 3))
            except:
                return None
        else:
            return None


    def gene2graph(self, df, df_expr, smp, pathway_name):
        """
        construct graph

        df: pathway network info, dataframe.
        df_expr: gene expression info, dataframe.
        smp: patient id
        pathway_name: pathway_name
        """
        # 按照Gene名（entrez_id)升序
        genes = df.src.to_list()
        genes = list(set(genes + df.dest.to_list()))
        genes = sorted(genes)
        
        ######### node feature #########
        # extract row with expression level
        df_tmp = df_expr[df_expr.Entrez_Gene_Id.isin(genes)].loc[:, ['Entrez_Gene_Id', smp]]
        
        # extract gene with expression level
        genes_tmp = df_tmp.Entrez_Gene_Id.to_list()

        # format gene without expression level
        df_1 = pd.DataFrame([list(set(genes).difference(set(genes_tmp))),
                    [df_tmp.iloc[:, 1].mean()]*len(list(set(genes).difference(set(genes_tmp))))
                    ]).T
        df_1.columns = df_tmp.columns

        df_2 = pd.concat([df_tmp, df_1], axis=0, ignore_index=True)
        df_2 = df_2.drop_duplicates('Entrez_Gene_Id', keep='first').sort_values('Entrez_Gene_Id')
        
        node_features = torch.tensor(df_2.loc[:, smp].values.reshape(-1,1))

        # 序号作为gene的图节点
        edge_index = np.array([[genes.index(i) for i in df.src.values],
                            [genes.index(i) for i in df.dest.values]])  
        edge_index = torch.tensor(np.unique(edge_index, axis=1),dtype=torch.long) # 删除重复边

        # 构建图
        if self.pathway_score == None:
            data = Data(edge_index=edge_index, x=node_features.to(torch.float32))
        else:
            score_ = self.get_score(smp=smp, pathway_name=pathway_name)
            data = Data(edge_index=edge_index, x=node_features.to(torch.float32), y=score_)

        return data


    def batch_gene2graph(self):
        ####  read gene expression file: Entrez_Gene_Id | sample_id | sample_id | ...
        df_expr = pd.read_table(self.gene_expression_file)
        df_expr = df_expr[~df_expr.Entrez_Gene_Id.isnull()]  # remove NULL
        df_expr.Entrez_Gene_Id = df_expr.Entrez_Gene_Id.astype(int)  # change Enrez Gene Name to int type
        df_expr = df_expr.sort_values('Entrez_Gene_Id')
        df_expr.iloc[:,1:] = (df_expr.iloc[:,1:] - df_expr.iloc[:,1:].min()) / (df_expr.iloc[:,1:].max() - df_expr.iloc[:,1:].min())    # change to 0-1
        # df_expr.iloc[:,1:] = (df_expr.iloc[:,1:] - df_expr.iloc[:,1:].mean()) / (df_expr.iloc[:,1:].std())    # change to z-score

        #### constrcut graph
        # dict: for save train data: sample  --> pathway
        dataset = {}  # 所有样本的图数据
                
        i = 0
        for smp in df_expr.columns[1:]:
            print('running {}:'.format(i), smp)
            pathway_names = []  #存储pathway名
            i += 1

            dat = []  # 单个患者的所有pathway图， list类型
            files = glob.glob(self.pathway_files)
            keep_pathways = pd.read_table(self.keep_pathway)['id'].to_list()
            for file in files:
                pathway_name = file.split("\\")[-1].split('.')[0]
                if pathway_name not in keep_pathways:
                    # print('{} not in keep pathways.'.format(pathway_name))
                    continue
                df = pd.read_table(file,sep=',')
                df = df[df.direction == 'directed']  ## 基因直接相互作用：在Graph中才有Edge
                data = self.gene2graph(df, df_expr=df_expr, smp=smp, pathway_name=pathway_name)

                # 去除节点小于15的图
                # if df.shape[0] < 15:
                #     continue
                # if data.num_nodes < 15:
                #     continue
                    
                dat.append(data)
                pathway_names.append(pathway_name)
        
            dataset[smp] = dat

        # save dataset
        torch.save(dataset, self.save_dataset)
        # output pathway name
        if os.path.exists(self.save_dataset + '.pathway_names.txt'):
            os.remove(self.save_dataset + '.pathway_names.txt')
        with open(self.save_dataset + '.pathway_names.txt', 'a') as F:
            for line in pathway_names:
                F.writelines(line + '\n')




###############################################################################
#                                                                             #
#                              Main function                                  #
#                                                                             #
###############################################################################
if __name__ == '__main__':

    os.environ["CUDA_VISIBLE_DEVICES"] = "0"

    #  Preprocess_pathway
    project_root_dir= "./"
    data_dir= "Dataset"
    Preprocess_pathway(gene_expression_file=os.path.join(project_root_dir, data_dir, 'expression_matrix_new.txt'),
                       pathway_files=os.path.join(project_root_dir, 'Pathway\\pathways/R*.txt'),
                       pathway_score=os.path.join(project_root_dir, data_dir, 'expression_matrix.ssgsea.csv'),
                       save_dataset =os.path.join(project_root_dir, data_dir, 'pathway_reactome_ssgsea.pt'),
                       keep_pathway = os.path.join(project_root_dir, 'Pathway\\pathways_keep.tsv')
    ).batch_gene2graph()
