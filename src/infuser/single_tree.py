from treelib import Tree, Node
from .Pixel import Pixel
from .init_tree import *
from .conversions import *
from .utils import *
import warnings
import os
import numpy as np
import time
from joblib import Parallel, delayed
#import pandas as pd

def prepare_bins(chrom_sizes, res, subset = None, chromlist = []):
    '''Prepare the bins to store Hi-C data
    
    This function takes the resolution and chomosome sizes and creates bins along the genome.
    
    Parameters
    ----------
    chrom_sizes : string
        The path to the file containing the chromosome sizes
    res : integer
        The resolution that was used, the width of the bins
    subset : string
        optional, the path to the file containing the regions to subset, if any (default None)
    chromlist : list of string
        optional, the names of the chromosomes to consider. This list is ignored if subset is not null. 
        If set to "all", all autosomes* and chrom X are considered. (default "all") *This default is set to humans,
        if applied to another organism, the chromosome list has to be specified.
    
    Returns
    -------
    the computed bins
    a dictionary with the id of the first bin of each chromosome
    '''
    first_id_dict = {}
    id = 0
    sizes = pd.read_csv(chrom_sizes, sep = "\t", header = None)
    bins = None
    if (subset is not None) and (subset != ""):
        subset_info = pd.read_csv(chrom_sizes, sep = "\t", header = None)
        chromosomes_to_consider = subset_info[0]
    elif len(chromlist) > 0:
        chromosomes_to_consider = chromlist
    else :
        chromosomes_to_consider = sizes[0]
    for chrom in chromosomes_to_consider:
        if (len(chromlist) > 0) and (chrom not in chromlist):
            continue
        first_id_dict[chrom] = id
        max_size = int(sizes[1][sizes[0] == chrom])
        start = range(0, max_size, res)
        end = [i for i in range(res, max_size, res)] + [max_size]
        chr_names = [chrom for i in start]
        df_chrom = pd.DataFrame({'chrom': chr_names, 'start': start, 'end':end})
        id += len(start)
        if bins is None:
            bins = df_chrom
        else :
            bins = pd.concat([bins, df_chrom])
    return bins, first_id_dict



def single_tree(tree_path, sample_file, output_dir, chrom_sizes, chromlist,\
    res = 10000, subset = None, dist = 0,\
    column = 4, transform = ["Z-score"], balance = True, n_jobs = 4):
    '''Find internal nodes of a phylogenic tree
    
    This function takes a tree topology and the data from its leaves to infer the state 
    of internal nodes through a modified Sankoff algorithm
    The function does not return anything but saves the computed internal nodes data in the output_dir.
    
    Parameters
    ----------
    tree_path : string
        the path to the file with the tree topology
    sample_file : string
        the path to the file with the paths to the samples
    output_dir : string
        the path to the output directory
    chrom_sizes : string
        the path to the file containing the size (in bp) of each chromosome
    chromlist : list of string
        the names of the chromosomes to consider. This list is ignored if subset is not null
    res : int
        optional, the resolution to consider
    subset : string
        optional, the path to the file containing the regions to subset, if any (default None)
    dist : int
        optional, the distance to consider. All interactions beyond that distance will be ignored. If set to 0, all interactions 
        are kept. (default 0)
    column : int
        optional, the column conting the score to consider. The first column is column 1. Ignored if the input files are .mcool files (default 4)
    transform : list of string
        the transformation(s) to apply to the matrix. "OE". "log1p" and "Z-score" are supported. (default Z-score)
    balance : boolean
        should the balanced weights be used (default True)
    n_jobs : int
        for paralleliation of pixel computation, how many jobs should be run in parallel (default 4)
    Raises
    ------
    Exception if one sample is of type .mcool and not the others
    
    '''
    start_time = time.time()
    
    # 1- Initialize the tree and create the output directory
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        log = open(output_dir+"/log.txt", "w")
        log.write("Initiation\n")
        log.close()
    else :
        log = open(output_dir+"/log.txt", "w")
        log.write("Initiation\n")
        log.write("The output directory " + output_dir +  " already exists. Overwriting...\n")
        log.close()
        warnings.warn("The output directory " + output_dir +  " already exists. Overwriting...")
    
    tree = init_tree(tree_path)
    
    # 2- Read the sample file and store the data into a dictionary. The keys should match the leaves' names in the tree
    vector_dict = {} # Dictionary containing the Hi-C data, stored as vectors
    means_dict = {} # Dictionary containing the means per diagonal
    vars_dict = {} # Dictionary containing the variances per diagonal
    num_HiC = 0 # Count the number of Hi-C samples to check they are all Hi-C or all other
    num_other = 0
    f = open(sample_file, "r")
    try:
        line = f.readline()
        
        while line != '':
            line = line.strip('\n')
            names = line.split('\t')
            # If the name is not present in the tree, raise a warning
            if not tree.contains(names[0]):
                log = open(output_dir+"/log.txt", "a")
                log.write("The sample " + names[0] +  " is not in the tree. Ignoring...\n")
                log.close()
                warnings.warn("The sample " + names[0] +  " is not in the tree. Ignoring...")
            else :
                log = open(output_dir+"/log.txt", "a")
                log.write("    Reading " + names[0] + "...\n")
                log.close()
                
                sample_is_cool = "cool" in names[1]
                if sample_is_cool:
                    num_HiC += 1
                    v, means, vars = HiC_matrix_to_vector(names[1], res, subset = subset, chromlist = chromlist, balance = balance, transformation = transform, dist = dist)
                else:
                    num_other += 1
                    v, means, vars = get_1D_data(names[1], transformation = transform, col = column)
                
                # Transform the Z-score vector to indices vector
                means_dict[names[0]] = means
                vars_dict[names[0]] = vars
                vector_dict[names[0]] = v
            line = f.readline()
    finally:
        f.close()
    # If we don't have only Hi-C or, only non-Hi-C data, we stop the function
    if (num_HiC != 0 ) and (num_other != 0):
        raise Exception('All samples do not have the same type.')
    init_done = time.time()
    log = open(output_dir+"/log.txt", "a")
    log.write("Initiation done. Time elapsed: %s seconds.\n" % (init_done - start_time))
    log.close()
    
    n_pixels = len(v)
    # Go through the nodes in the tree. For each leave, check that it's in the dictionary. If not, raise an error.
    # For non-leaves nodes, create an empty vector that will store the data
    log = open(output_dir+"/log.txt", "a")
    log.write("Map inference\n")
    log.close()
    for node in tree.all_nodes():
        if node.is_leaf() & (not node.tag in vector_dict.keys()):
            raise Exception("One leaf present in the tree is absent from sample_file.\n")
        elif not node.is_leaf():
            vector_dict[node.tag] = np.zeros(n_pixels)
            
    # 3- Iterate across pixels and save the data sequentially
    percent = np.round(n_pixels/100)

    # Convert the dictionary in an array fo easy parallelization
    matrix = np.array(list(vector_dict.values()))
    row_names = np.array(list(vector_dict.keys()))
    results = Parallel(n_jobs=n_jobs)(delayed(process_pixel)(pixel, row_names, Tree(tree, deep = True)) for pixel in matrix.transpose())
    results = np.array(results).transpose()
    ########

    body_done = time.time()
    log = open(output_dir+"/log.txt", "a")
    log.write("Map inference done. Time elapsed: %s seconds.\n" % (body_done - init_done))
    log.close()

    # 3.5 Stats about the parsimony score
    log = open(output_dir+"/log.txt", "a")
    parsimony_scores = results[-1]
    log.write("Total parsimony score: %s \n" % (parsimony_scores.sum()))
    log.write("Mean parsimony score: %s \n" % (parsimony_scores.mean()))
    log.write("Standard deviation parsimony score: %s \n" % (parsimony_scores.std()))
    log.close()
        
    # 4- Save the results in the right format
    log = open(output_dir+"/log.txt", "a")
    log.write("Saving the files\n")
    log.close()

    if num_HiC > 0:
        bins_df, id_dict = prepare_bins(chrom_sizes, res, subset, chromlist)
    # Save tree structure
    if os.path.exists(output_dir + "/tree_structure.txt"):
        os.remove(output_dir + "/tree_structure.txt")
    tree.save2file(output_dir + "/tree_structure.txt")
    
    # Save the data in matrix format
    if not os.path.isdir(output_dir+"/leaves"): os.mkdir(output_dir+"/leaves")
    if not os.path.isdir(output_dir+"/internal_nodes"): os.mkdir(output_dir+"/internal_nodes")
    for node in vector_dict.keys():
        # Re-convert to Z-score
        row_num = np.where(row_names == node)[0][0]
        vector_dict[node] = results[row_num]
        log = open(output_dir+"/log.txt", "a")
        log.write("    Saving " + node + "...\n")
        log.close()
        
        if tree.get_node(node).is_leaf():
            path = output_dir + "/leaves/"
            samples_to_consider = [node]
        else :
            path = output_dir + "/internal_nodes/"
            samples_to_consider = [c.tag for c in tree.subtree(node).leaves()]
        
        if num_HiC > 0:
            mat_dict = vector_to_matrix(vector_dict[node], res, subset, dist, chrom_sizes, chromlist)
                
            pixel_df = None
            zscores_df = None
            for chrom in mat_dict.keys():
                mat = to_HiC_scores(mat_dict[chrom], chrom, samples_to_consider, transform, means_dict, vars_dict)
                names = range(id_dict[chrom], id_dict[chrom]+len(mat))
                df = pd.DataFrame(mat,columns=names) # Pseudo counts
                z_df = pd.DataFrame(mat_dict[chrom], columns = names) # Z-scores
                # Pseudo-counts
                df['bin1_id'] = names
                df = pd.melt(df,'bin1_id',var_name='bin2_id',value_name='count')
                df = df[[not pd.isna(i) for i in df['count']]]
                df = df.sort_values(by = ['bin1_id', 'bin2_id'])
                df['bin2_id'] = df['bin2_id'].astype('int64')
                if pixel_df is None:
                    pixel_df = df
                else :
                    pixel_df = pd.concat([pixel_df, df])
                # Z-scores
                z_df['bin1_id'] = names
                z_df = pd.melt(z_df,'bin1_id',var_name='bin2_id',value_name='count')
                z_df = z_df[[not pd.isna(i) for i in z_df['count']]]
                z_df = z_df.sort_values(by = ['bin1_id', 'bin2_id'])
                z_df['bin2_id'] = z_df['bin2_id'].astype('int64')
                if zscores_df is None:
                    zscores_df = z_df
                else :
                    zscores_df = pd.concat([zscores_df, z_df])
            # Save pseudo counts
            cooler.create_cooler(path+node+".cool", bins_df, pixel_df, ordered=True, dtypes = {"count": "float"})
        else :
            data = to_scores(vector_dict[node], samples_to_consider, transform, means_dict, vars_dict)
            matrices_dict[node] = vector_dict[node]
            np.savetxt(path+node+".tsv", data, delimiter="\t")

    save_nodes_done = time.time()
    log = open(output_dir+"/log.txt", "a")
    log.write("Nodes' data saved. Time elapsed: %s seconds.\n" % (save_nodes_done - body_done))
    log.close()
    
    # Save the differences between each node and its children in matrix format
    if not os.path.isdir(output_dir+"/edges"): os.mkdir(output_dir+"/edges")
    for new_node in tree.all_nodes():
        children = tree.children(new_node.tag)
        if len(children) > 0:
            for child in children:
                pair = new_node.tag + "_to_" + child.tag
                path = output_dir + "/edges/"
                if not os.path.isdir(path): os.mkdir(path)
                log = open(output_dir+"/log.txt", "a")
                log.write("    Saving " + pair + "...\n")
                log.close()
                
                if num_HiC > 0:
                    pixel_df = None
                    vector_diff = vector_dict[child.tag] - vector_dict[new_node.tag]
                    mat_dict = vector_to_matrix(vector_diff, res, subset, dist, chrom_sizes, chromlist)
                    for chrom in mat_dict.keys():
                        diff = mat_dict[chrom]
                        names = range(id_dict[chrom], id_dict[chrom]+len(diff))
                        df = pd.DataFrame(diff,columns=names)
                        df['bin1_id'] = names
                        df = pd.melt(df,'bin1_id',var_name='bin2_id',value_name='count')
                        df = df[[not pd.isna(i) for i in df['count']]]
                        df = df.sort_values(by = ['bin1_id', 'bin2_id'])
                        df['bin2_id'] = df['bin2_id'].astype('int64')
                        if pixel_df is None:
                            pixel_df = df
                        else :
                            pixel_df = pd.concat([pixel_df, df])
                    cooler.create_cooler(path+pair+".cool", bins_df, pixel_df, ordered=True, dtypes = {"count": "float"})
                else:
                    diff = matrices_dict[child.tag] - matrices_dict[new_node.tag]
                    np.savetxt(path+pair+".tsv", diff, delimiter="\t")
    save_edges_done = time.time()
    log = open(output_dir+"/log.txt", "a")
    log.write("Edges' data saved. Time elapsed: %s seconds.\n" % (save_edges_done - save_nodes_done))
    log.close()
    # end
