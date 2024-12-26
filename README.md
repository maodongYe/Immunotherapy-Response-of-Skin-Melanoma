The expression matrix needs to be processed by three Python programs in sequence: 1_change_pathway_to_gmt_format_for_ssgsea, 2_data, and 3_data_for_GNN, and then a new pathway-based gene expression data file is obtained.

This file and the clinical information file in txt format are used as inputs for 4_predict.

melamon_predict_model.pt is the trained model for predicting immune therapy response.
