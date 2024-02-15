from treelib import Tree, Node
from .Pixel import Pixel
import os

def init_tree(file_path, n_values = 9):
    '''Create a tree from a file containing its topology
    
    This function reads a file containing a tree topology and creates a Tree from it.
    The tree needs to be rooted.
    
    Parameters:
    file_path: string
        the path to the file with the tree topology
    n_values: int
        the number of values that can be stored in the nodes (default 9)
    
    Returns:
    a Tree object with the correct topology
    
    Raises:
    ValueError if n_values is not an integer
    FileNotFoundError if eiher the topology file is not found
    Exception if the tree topology does not respect the expected format
    '''
    if not isinstance(n_values, int):
        raise ValueError("The number of possible values is not an integer")
    tree = Tree()
    # Read file and add the nodes consecutively. If one node should be added but its parent
    # does not exist, raise an Expection
    try:
        f = open(file_path, "r")
        line = f.readline().strip('\n')
        tree.create_node(line, line)
        
        line = f.readline()
        while line != '':
            line = line.strip('\n')
            names = line.split('\t')
            if not tree.contains(names[1]):
                raise Exception("The tree file does not respect the expected format.")
            tree.create_node(names[0], names[0], parent = names[1])
            line = f.readline()
    finally:
        f.close()
    
    # add Pixels to the Nodes' data
    for node in tree.all_nodes():
        node.data = Pixel(len(tree.children(node.tag)), n_values)
    
    return tree
    

