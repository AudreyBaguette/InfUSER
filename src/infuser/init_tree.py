from treelib import Tree, Node
from .Pixel import Pixel
import os
import warnings
import numpy as np

def init_tree(file_path):
    '''Create a tree from a file containing its topology
    
    This function reads a file containing a tree topology and creates a Tree from it.
    The tree needs to be rooted.
    
    Parameters:
    file_path: string
        the path to the file with the tree topology
    
    Returns:
    a Tree object with the correct topology
    
    Raises:
    ValueError if n_values is not an integer
    FileNotFoundError if eiher the topology file is not found
    Exception if the tree topology does not respect the expected format
    '''
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
        node.data = Pixel(len(tree.children(node.tag)))
    
    return tree

def find_new_interval(curr_node_children):
    '''Find the new interval of the node
    
    This function takes a node and fins the interval of values with the lowest combined difference with all children
    
    Parameters
    ----------
    curr_node_children : list of Node
        the nodes from which the new interval should be computed
    
    Returns
    -------
    a tuple containing the new interval
    '''
    # Go through all children and create a vector of possible values based on their intervals
    possible_values = []
    for child in curr_node_children:
        val1, val2 = child.data.get_interval()
        possible_values = possible_values + [val1, val2]
    #
    # For each possible value, compute the distance to each node and store the score in a new list
    scores = []
    possible_values = np.unique(possible_values)
    for pos_val in possible_values:
        total_score = 0
        for child in curr_node_children:
            new_score, x = child.data.compute_distance(pos_val)
            total_score = total_score + new_score
        #
        scores = scores + [total_score]
    #
    # Find what values lead to the minimal score
    candidates = possible_values[np.where(scores == np.min(scores))[0]]
    # TODO remove when dev phase is over
    if len(candidates) > 2:
        warnings.warn("More than 2 candidate values: " + candidates)
    #
    return((min(candidates), max(candidates)))


def process_pixel(array, row_names, new_tree):
    '''Process a single value with a modified Sankoff algo
    
    This function takes as input an array with the values of the initial tree and updates it by following
    a modified Sankoff algorithm
    
    Parameters
    ----------
    array : numpy array
        array containing the values for the current pixel, in the same order as the now_names
    row_names : numpy array
        the name of each node in the tree
    new_tree : Tree
        copy of the differentiation tree
    
    '''
    # Initialize the leaves
    for leaf in new_tree.leaves():
        index = np.where(row_names == leaf.tag)[0][0]
        leaf.data.set_interval(array[index], array[index])
    
    # Bottom-up
    # Iterate throught the non-leaves nodes
    nodes = new_tree.all_nodes()
    nodes.reverse()
    for node in nodes:
        if not node.is_leaf():
            # Update the interval of the current node
            inter_val1, inter_val2 = find_new_interval(new_tree.children(node.tag))
            node.data.set_interval(inter_val1, inter_val2)
    
    # Traceback
    # 1- Start at the root, pick a value in the interval and start scoring.
    min_value, max_value = new_tree.get_node(new_tree.root).data.get_interval()
    #selected_value = np.random.uniform(min_value, max_value) #instead of taking a random value, pick the middle of the interval, to make the algorithm deterministic
    selected_value = min_value + (max_value - min_value)/2
    index = np.where(row_names == new_tree.get_node(new_tree.root).tag)[0][0]
    array[index] = selected_value
    parsimony_score = 0
    # 2- Loop through nodes, top to bottom, set value and update score
    for node in new_tree.all_nodes():
        # For each child, compute the new value and update the score
        children = new_tree.children(node.tag)
        current_value = array[np.where(row_names == node.tag)[0][0]]
        
        for child in children:
            new_score, new_value = child.data.compute_distance(current_value)
            index = np.where(row_names == child.tag)[0][0]
            array[index] = new_value
            parsimony_score = parsimony_score + new_score
    
    return(np.concatenate((array, [parsimony_score])))

