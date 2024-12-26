import pandas as pd
import os
os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"



###############################################################################
#                                                                             #
#                        Change pathway to gmt format for ssgsea              #
#                                                                             #
###############################################################################
def generate_gmt(dtail_file, save_file):
    """
    get gmt format pathway file for ssgsea.
    
    dtail_file: keep_pathways_details.tsv
    save_file: gmt pathway file
    """
    df_dtails = pd.read_table(dtail_file,sep=',')

    # remove old file
    with open(save_file, 'w') as F:
        F.writelines('')

    for idx, row in df_dtails.iterrows():
        id_ = row[0]
        df = pd.read_table('./Pathway/pathways/{}.txt'.format(id_),sep=',')
        genes = df.src.to_list()
        genes = list(set(genes + df.dest.to_list()))
        genes = sorted(genes)
        num_genes = len(genes)
        df_dtails.loc[idx, 'Number_of_nodes_update'] = int(num_genes)
        # output
        line = '\t'.join([id_] + [str(g) for g in genes])
        with open(save_file, 'a') as F:
            F.writelines(line + '\n')

    df_dtails.to_csv(dtail_file, index=False, sep='\t')


###############################################################################
#                                                                             #
#                              Main function                                  #
#                                                                             #
###############################################################################
if __name__ == '__main__':


    # step1: get gmt pathway file
    generate_gmt(dtail_file='./Pathway/pathways_keep.tsv',
                 save_file='./Pathway/pathway_genes_list.txt')
