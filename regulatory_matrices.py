
import sys
import os


def read_cell_types(link):
    with open(link, 'r') as t:
        lines=t.readlines()
        t_types=[]
        for x in lines:
            t_types.append(x.split()[3])
        t.close()  
    return t_types


def main():
    t_type_link = "/data/zusers/pratth/ATAC/specific-elements/top-10k/unstimulated_T-cells.bed"
    t_types = read_cell_types(t_type_link) # reads a cell type matrix
    print("Number of regulatory elements in T-cells", end='')
    print(len(t_types))

    b_type_link = "/data/zusers/pratth/ATAC/specific-elements/top-10k/B-cell.bed"
    b_types = read_cell_types(b_type_link) # reads a cell type matrix
    print("Number of regulatory elements in B-cells", end='')
    print(len(b_types))

    m_type_link = "/data/zusers/pratth/ATAC/specific-elements/top-10k/myeloid_cells.bed"
    m_types = read_cell_types(m_type_link) # reads a cell type matrix
    print("Number of regulatory elements in Myeloid cells", end='')
    print(len(m_types))

if __name__ == "__main__":
    sys.exit(main())