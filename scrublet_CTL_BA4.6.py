import scrublet as scr
import numpy as np

counts_matrix_link = "/data/rusers/sheddn/UCLA-ASD/data/CTL_BA4.6_counts.txt"
genes_list_link = "/data/rusers/sheddn/UCLA-ASD/data/CTL_BA4.6_genes.txt"

counts_matrix = np.array(counts_matrix_link, delimiter='\t')
genes_list = np.array(genes_list_link, delimiter='\t')

scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.15)

doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)

np.savetxt("/data/rusers/sheddn/UCLA-ASD/data/CTL_BA4.6_doublets.txt", predicted_doubles, sep='\t')
