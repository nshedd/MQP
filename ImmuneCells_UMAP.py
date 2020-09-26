#!/usr/bin/env python

import sys
import os
import umap

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot

def read_matrix():
    matrix = [] # creates an empty array which will hold matrix values when populated
    cells = [] # creates an empty array which will hold cell IDs when populated
    with open("/data/zusers/pratth/sc/atac/PBMC/PBMC.matrix.txt", 'r') as f:
        elements = f.readline().strip().split() # creates an array with element IDs by splitting the first line at each tab
        for line in f: # for the remaining lines...
            fields = line.strip().split() # split the line into fields at each tab
            cells.append(fields[0]) # the cell ID is in the first column
            matrix.append([ float(x) for x in fields[1:] ]) # this converts all the fields from index 1 onward to floating point numbers and adds them to the matrix
    return elements, cells, matrix

def read_cell_types():
    with open("/data/zusers/pratth/ATAC/specific-elements/top-10k/unstimulated_t-cells.bed", 'r') as t:
        elements = t.readline().strip().split() # creates an array with element IDs by splitting the first line at each tab
        lines=t.readlines()
        t_types=[]
        for x in lines:
            t_types.append(x.split(' ')[3])
        f.close()  
    return t_types

def match_types(cells, t_types):
    c=[]
    for i in range(0,len(cells),1):
        if cells[i] in t_types:
            c.append("blue")
        else:
            c.append("red")
    return c

def main():
    elements, cells, matrix = read_matrix() # reads the matrix from the file
    t_types = read_cell_types() # reads a cell type matrix
    u = umap.UMAP(n_neighbors = 10, min_dist = 0.1, metric = 'euclidean') # initialize UMAP. different parameters might give better separation
    coordinates = u.fit_transform(matrix) # perform the transformation. outputs a list of 2D coordinates, one for each row
    matplotlib.pyplot.scatter(
        [ x for x, y in coordinates ], # extract the x-coordinates from the UMAP output
        [ y for x, y in coordinates ], # extract the y-coordinates from the UMAP output
        marker = '.', # make the points small so the plot isn't too crowded
        alpha = 0.1, # make the points semi-transparent so it is easier to tell where points densely cluster together
        c = match_types(cells, t_types) # this makes every point dark red (except for unstimulated t cells). TODO: replace with coloring by marker elements
    )
    matplotlib.pyplot.savefig(os.path.expanduser("~/umap.svg")) # write the plot to "umap.svg" in your home directory
    return 0

if __name__ == "__main__":
    sys.exit(main())
