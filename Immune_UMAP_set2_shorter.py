import sys
import os
import umap
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot

def read_matrix():
    matrix = [] # creates an empty array which will hold matrix values when populated
    cells = [] # creates an empty array which will hold cell IDs when populated
    with open("/data/zusers/pratth/sc/atac/GSM3722075_PBMC_Rep3_fragments.tsv.gz.rDHS.matrix.tsv", 'r') as f:
        elements = f.readline().strip().split() # creates an array with element IDs by splitting the first line at each tab
        for line in f: # for the remaining lines...
            fields = line.strip().split() # split the line into fields at each tab
            cells.append(fields[0]) # the cell ID is in the first column
            matrix.append([ float(x) for x in fields[1:] ]) # this converts all the fields from index 1 onward to floating point numbers and adds them to the matrix
    return matrix


def normalize_data(data):
    n_data = []
    for i in data:
        sum = 0
        n_i = []
        for j in i:
            sum = sum + j
        for j in i:
            n_j = j/sum
            n_i.append(n_j)
        n_data.append(n_i)
    return n_data


def main():
    matrix = read_matrix() # reads the matrix from the file

    n_matrix = normalize_data(matrix)

    colors = np.loadtxt(os.path.expanduser("~/set2_top10k_colors.npy"), dtype=str)

    u = umap.UMAP(metric = 'euclidean') # initialize UMAP. different parameters might give better separation
    coordinates = u.fit_transform(n_matrix) # perform the transformation. outputs a list of 2D coordinates, one for each row
    #colors = match_types(elements, t_types)
    matplotlib.pyplot.scatter(
        [ x for x, y in coordinates ], # extract the x-coordinates from the UMAP output
        [ y for x, y in coordinates ], # extract the y-coordinates from the UMAP output
        marker = '.', # make the points small so the plot isn't too crowded
        alpha = 0.1, # make the points semi-transparent so it is easier to tell where points densely cluster together
        c = colors # this makes unstimulated t cells blue and everything else black. TODO: replace with coloring by marker elements
    )
    matplotlib.pyplot.title("<Dataset 2> UMAP, n_neighbors=default, min_dist=default")
    matplotlib.pyplot.savefig(os.path.expanduser("~/umap_colored_set2_default.svg")) # write the plot to "umap.svg" in your home directory
    return 0

if __name__ == "__main__":
    sys.exit(main())
