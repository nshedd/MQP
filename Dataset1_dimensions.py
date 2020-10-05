import sys
import os

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

def count_cells(cells):
    count = 0
    for i in cells:
        count = count + 1
    return count


def count_elements(elements):
    count = 0
    for i in elements:
        count = count + 1
    return count


def count_total(matrix):
    count = 0
    for i in matrix:
        for j in i:
            count = count + 1
    return count

def main():
    elements, cells, matrix = read_matrix() # reads the matrix from the file

    elements_count = count_elements(elements)
    print("Number of elements: ", end = '')
    print(elements_count)

    cells_count = count_cells(cells)
    print("Number of cells: ", end = '')
    print(cells_count)

    total_count = count_total(matrix)
    print("Number of measurements: ", end = '')
    print(total_count)

if __name__ == "__main__":
    sys.exit(main())
