#!/usr/bin/env python

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
    return elements, cells, matrix

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

def read_cell_types(link):
    with open(link, 'r') as t:
        lines=t.readlines()
        t_types=[]
        for x in lines:
            t_types.append(x.split()[3])
        t.close()  
    return t_types

def match_types(elements, types):
    c=[]
    for i in elements:
        if i in types:
            c.append("blue")
        else:
            c.append("black")
    return c

def create_histograms(matrix, cells, elements, marker, marker_name):
    sums=[]
    indices = []
    for j in elements:
        if j in marker:
            indices.append(elements.index(j))
    for i in range(0, len(cells), 1):
        marker_sum = 0
        for j in indices:
            if elements[j] in marker:
                marker_sum = marker_sum + matrix[i][j]
        if marker_sum > 0.1:
            sums.append(marker_sum)	
    if marker_name=="t_types":
        matplotlib.pyplot.hist(sums, bins=50)
        matplotlib.pyplot.savefig(os.path.expanduser("~/marker_colorsnt1.svg"))
        matplotlib.pyplot.close()
        print("done with t_cells")
    if marker_name=="b_types":
        matplotlib.pyplot.hist(sums, bins=50)
        matplotlib.pyplot.savefig(os.path.expanduser("~/marker_colorsnb1.svg"))
        matplotlib.pyplot.close()
        print("done with b_cells")
    if marker_name=="m_types":
        matplotlib.pyplot.hist(sums, bins=50)
        matplotlib.pyplot.savefig(os.path.expanduser("~/marker_colorsnm1.svg"))
        matplotlib.pyplot.close()
        print("done with m_cells")

def color_graph(matrix, cells, elements, marker, colors, marker_name):
    sum=[]
    indices = []
    for j in elements:
        if j in marker:
            indices.append(elements.index(j))
    for i in range(0, len(cells), 1):
        marker_sum = 0
        for j in indices:
            if elements[j] in marker:
                marker_sum = marker_sum + matrix[i][j]
        if marker_name == "t_types":
            if marker_sum > 0.005:
                colors[i]="red"
        if marker_name == "b_types":
            if marker_sum > 0.01:
                if colors[i] == "red":
                    colors[i]="purple"
                else:
                    colors[i]="blue"
        if marker_name == "m_types":
            if marker_sum > 0.01:
                if colors[i] == "red":
                    colors[i]="orange"
                elif colors[i] == "blue":
                    colors[i]="green"
                elif colors[i] == "purple":
                    colors[i]="brown"
                else:
                    colors[i]="yellow"
    return colors


def main():
    elements, cells, matrix = read_matrix() # reads the matrix from the file

    n_matrix = normalize_data(matrix)

    colors = ["black"] * len(cells)

    t_type_link = "/data/zusers/pratth/ATAC/specific-elements/top-10k/unstimulated_T-cells.bed"
    t_types = read_cell_types(t_type_link) # reads a cell type matrix
    #colors = color_graph(n_matrix, cells, elements, t_types, colors, "t_types") ## later add color_list to the input and just alter the color at an index if it is in a threshold
    create_histograms(n_matrix, cells, elements, t_types, "t_types")

    b_type_link = "/data/zusers/pratth/ATAC/specific-elements/top-10k/B-cell.bed"
    b_types = read_cell_types(b_type_link) # reads a cell type matrix
    #colors = color_graph(n_matrix, cells, elements, b_types, colors, "b_types")
    create_histograms(n_matrix, cells, elements, b_types, "b_types")

    m_type_link = "/data/zusers/pratth/ATAC/specific-elements/top-10k/myeloid_cells.bed"
    m_types = read_cell_types(m_type_link) # reads a cell type matrix
    #colors = color_graph(n_matrix, cells, elements, m_types, colors, "m_types")
    create_histograms(n_matrix, cells, elements, m_types, "m_types")

    u = umap.UMAP(n_neighbors = 10, min_dist = 0.1, metric = 'euclidean') # initialize UMAP. different parameters might give better separation
    coordinates = u.fit_transform(n_matrix) # perform the transformation. outputs a list of 2D coordinates, one for each row
    #colors = match_types(elements, t_types)
    matplotlib.pyplot.scatter(
        [ x for x, y in coordinates ], # extract the x-coordinates from the UMAP output
        [ y for x, y in coordinates ], # extract the y-coordinates from the UMAP output
        marker = '.', # make the points small so the plot isn't too crowded
        alpha = 0.1, # make the points semi-transparent so it is easier to tell where points densely cluster together
        c = colors # this makes unstimulated t cells blue and everything else black. TODO: replace with coloring by marker elements
    )
    matplotlib.pyplot.savefig(os.path.expanduser("~/umap_colored.svg")) # write the plot to "umap.svg" in your home directory
    return 0

if __name__ == "__main__":
    sys.exit(main())
