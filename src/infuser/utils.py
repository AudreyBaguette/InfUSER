import numpy as np
import pandas as pd
import cooler
import math
from cooltools.lib import numutils


def kth_diag_indices(a, k):
    '''Get the indices of the kth diagonal
    
    This function takes matrix a and a number k to return the indices (rows and columns) of the kth
    diagonal of a.
    
    Parameters
    ----------
    a : numpy array
        a square matrix from which we want the diagonals
    k : int
        the skew, from the normal diagonal
    
    Returns
    -------
    the row indices of the kth diagonal
    the column indices of the kth diagonal

    '''
    rows, cols = np.diag_indices_from(a)
    if k < 0:
        return rows[-k:], cols[:k]
    elif k > 0:
        return rows[:-k], cols[k:]
    else:
        return rows, cols

def get_middle(array, D):
    '''Get the values from the triangular matrix, up to distance D
    
    This function takes matrix and a distance D, then returns the flattened values of the upper triangular 
    separated by at most D elements.
    
    Parameters
    ----------
    array : numpy array
        a square matrix from which we want the values
    D : int
        the maximal distance to consider
    
    Returns
    -------
    numpy array with the values of interest

    '''
    # Get the indices of the upper triangular and the indices to remove
    iu1 = np.triu_indices_from(array)
    iu2 = np.triu_indices_from(array, D)
    # Extract upper triangular
    triu = array[iu1]
    # Convert 2D indices in 1D indices
    triu1_indices = iu1[0]*len(array) + iu1[1]
    triu2_indices = iu2[0]*len(array) + iu2[1]
    # Compute what to remove
    to_remove = np.isin(triu1_indices, triu2_indices)
    vector = np.delete(triu, to_remove)
    
    return vector


def diag_to_matrix(vector, D, size):
    '''Transforms  matrix subset to a matrix
    
    This function takes a flattened vector, as created by get_middle(), and reforms it into
    the initial matrix
    
    Parameters
    ----------
    vector : numpy array
        the values to re-arrange into a matrix
    D : int
        the maximal distance that is to consider between two values
    size: int
        the matrix size (length of each dimension)
    
    Returns
    -------
    a numpy array with the input data converted into a matrix

    '''
    A = np.ones((size, size))
    A[np.triu_indices(size, D)] = np.nan
    A[np.tril_indices(size, -1)] = np.nan
    A[A == 1] = vector
    
    return(A)


def z_transform(matrix):
    '''Applies a z-score transformation on the matrix
    
    This function applies a z-score, stratified by distance, on the counts matrix.
    Z_ij = (B_ij - E_d(i,j)) / sqrt(Var_d(i,j))
    B_ij: (balanced) contact frequency between pair of bins i, j
    E_d(i,j): mean (balanced) contact frequency between all pairs of distance d, that is |i-j| = mean of all B_ij at distance d
    Var_d(I,j): variance of B_ij at distance d
    
    Parameters
    ----------
    matrix : numpy array
        the matrix to transform
    
    Returns
    -------
    a numpy array that is the z-score transformed matrix
    a numpy array that contains the means of the diagonals
    a numpy array that contains the variances of the diagonals

    '''
    means = np.zeros(len(matrix))
    vars = np.zeros(len(matrix))
    for diag in range(0, len(matrix)):
        mean = np.mean([matrix.diagonal(diag), matrix.diagonal(-diag)])
        var = np.var([matrix.diagonal(diag), matrix.diagonal(-diag)])
        matrix[kth_diag_indices(matrix, diag)] = (matrix.diagonal(diag) - mean) / np.sqrt(var)
        matrix[kth_diag_indices(matrix, -diag)] = (matrix.diagonal(-diag) - mean) / np.sqrt(var)
        means[diag] = mean
        vars[diag] = var
    
    return matrix, means, vars


def HiC_matrix_to_vector_chrom(sample, res, chrom, subset = None, balance = True, transformation = ["Z-score"], dist = 0):
    '''Convert a HiC matrix to a vector, for one given chromosome
    
    This function load the Hi-C matrix of a sample, for a given chromosome, and converts it to a vector by flattening it. 
    If a subset file is specified, only regions within the subset are considered. If multiple regions are specified within 
    the subset file, they are considered separately and the interaction across the regions is ignored. 
    Notes:
    - this function assumes that multiple resultions are stored in the .mcool file and they can be accessed with "::/resolutions/". 
    It can read .cool files, but will not check the file's resolution matches the parameter resolution.
    - the start region from the subset file is rounded down to the nearest bin, relative to the resolution.
    - the end region from the subset file is rounded down to the nearest bin, relative to the resolution. This bin is excluded.
    - if the start and end region fall within the exact same bin, the region is considered too small and is ignored.
    
    Parameters
    ----------
    sample : string
        the path to the .mcool file
    res : int
        the resolution to consider
    chrom : string
        the chomosome to consider
    subset : string
        optional, the path to the file containing the regions to subset, if any (default None)
    balance : boolean
        optional, should the balanced matrix be used instead of the raw counts? (default True)
    transformation : list of string
        the transformation(s) to apply to the matrix. "OE". "log1p" and "Z-score" are supported. (default Z-score)
    dist : int
        the distance to consider. All interactions beyond that distance will be ignored. If set to 0, all interactions 
        are kept. (default 0)
    
    Returns
    -------
    a numpy array that is the vectorized (subset of the) matrix
    a numpy array that contains the means of the diagonals (None if not z-transformed)
    a numpy array that contains the variances of the diagonals  (None if not z-transformed)
    
    Raises
    ------
    Warning if the start and end region from the subset file fall in the same bin
    FileNotFoundError if eiher the sample file or the susbet file is not found
    KeyError if the resolution is not present in the mcool file
    ValueError if the chromosome is not present in the mcool file
    NameError if the subset file is specified but does not contain the considered chromosome
    Exception if the vector is empty
    Exception if the vector does not have the expected length

    '''
    if '.mcool' in sample:
        clr = cooler.Cooler(sample+'::/resolutions/'+str(res))
    else:
        clr = cooler.Cooler(sample)
    vector = np.array([])
    
    A = np.array(clr.matrix(balance=balance).fetch(chrom))
    A[~np.isfinite(A)] = 0
    
    means = None
    vars = None
    if "OE" in transformation:
        # get O/E matrix
        mask = A.sum(axis=0) > 0
        OE, _, _, _ = numutils.observed_over_expected(A, mask)
        OE = np.clip(OE, 0, np.percentile(OE[mask, :][:, mask], 99.9))
        A = OE
        
    if "log1p" in transformation:
        A = np.log1p(A)

    if ("Z-score" in transformation) or ("z-score" in transformation) or ("zscore" in transformation) or ("Zscore" in transformation):
        A, means, vars = z_transform(A)

    
    if (not subset is None) and (subset != ""):
        regions = pd.read_csv(subset, sep = "\t")
        regions_chrom = regions[regions.seqnames == chrom]
        if len(regions_chrom) == 0:
            raise NameError('The subset file is not null and the chromosome ' + chrom + ' is not present in the subset file.')
        
        expected_length = 0
        for i in regions_chrom.index:
            start = int(regions_chrom.start[i])
            end = int(regions_chrom.end[i])
            start_i = math.floor(start/res)
            end_i = math.floor(end/res)
            if start_i == end_i:
                raise Warning(regions_chrom.end[i] + ' falls in the same bin as ' + regions_chrom.start[i] + ' at a ' + str(res) + 'bp resolution.')
            
            A_subset = A[start_i:end_i, start_i:end_i]
            vector = np.concatenate((vector, A_subset[np.triu_indices_from(A_subset)]))
            
            added_length = (end_i-start_i)
            expected_length += (added_length*(added_length+1))/2
    elif dist > 0:
        D = np.ceil(dist/res) - 1
        vector = get_middle(A, D)
        expected_length = (len(A)*(len(A)+1)/2) - ((len(A)-D)*((len(A)-D)+1)/2)
    else :
        vector = np.concatenate((vector, A[np.triu_indices_from(A)]))
        expected_length = (len(A)*(len(A)+1))/2
        #vector = np.concatenate((vector, A.flatten()))
        #expected_length = len(A) * len(A)
    
    if len(vector) == 0:
        raise Exception('The resulting vector is empty.')
    if len(vector) != expected_length:
        raise Exception('The resulting vector does not have the expected length.')
    return vector, means, vars


def HiC_matrix_to_vector(sample, res, subset = None, balance = True, transformation = ["Z-score"], dist = 0, chromlist = []):
    '''Convert a HiC matrix to a vector
    
    This function load the Hi-C matrix of a sample, and converts it to a vector by flattening it. 
    If a subset file is specified, only regions within the subset are considered. If multiple regions are specified within 
    the subset file, they are considered separately and the interaction across the regions is ignored. 
    Notes:
    - this function assumes that multiple resultions are stored in the .mcool file and they can be accessed with "::/resolutions/".
    - the start region from the subset file is rounded down to the nearest bin, relative to the resolution.
    - the end region from the subset file is rounded down to the nearest bin, relative to the resolution. This bin is excluded.
    - if the start and end region fall within the exact same bin, the region is considered too small and is ignored.
    
    Parameters:
    sample : string
        the path to the .mcool file
    res : int
        the resolution to consider
    subset : string
        optional, the path to the file containing the regions to subset, if any (default None)
    balance : boolean
        optional, should the balanced matrix be used instead of the raw counts? (default True)
    transformation : list of string
        the transformation(s) to apply to the matrix. "OE". "log1p" and "Z-score" are supported. (default Z-score)
    dist : int
        the distance to consider. All interactions beyond that distance will be ignored. If set to 0, all interactions 
        are kept. (default 0)
    chromlist : list of string
        optional, the names of the chromosomes to consider. This list is ignored if subset is not null. 
        If set to an empty list, all available chromosomes are considered. (default: [])
    
    Returns
    -------
    a numpy array that is the vectorized (subset of the) matrix
    a dictionary that contains one vector of means of the diagonals per chromosome (None if not z-transformed)
    a dictionary that contains oen vector of variances of the diagonals per chromosome  (None if not z-transformed)
    
    Raises
    ------
    EOFError if the subset file is specified but empty.

    '''
    full_vector = np.array([])
    all_means = {}
    all_vars = {}
    if (not subset is None) and (subset != ""):
        regions = pd.read_csv(subset, sep = "\t")
        if len(regions) == 0:
            raise EOFError('The subset file is empty.')
        all_chroms = list(set(regions.seqnames))
    elif len(chromlist) > 0:
        all_chroms = chromlist
    else:
        if '.mcool' in sample:
            clr = cooler.Cooler(sample+'::/resolutions/'+str(res))
        else:
            clr = cooler.Cooler(sample)
        all_chroms = clr.chromnames
    
    for chrom in all_chroms:
        vector, means, vars = HiC_matrix_to_vector_chrom(sample, res, chrom, subset, balance, transformation, dist)
        full_vector = np.concatenate((full_vector, vector))
        all_means[chrom] = means
        all_vars[chrom] = vars
            
    return full_vector, all_means, all_vars

def get_1D_data(sample, transformation = ["Z-score"], col = 4):
    '''Reads a file with binned values into a vector
    
    This function reads a bed-like file for a sample, and converts it to a vector. 
    If a subset file is specified, only regions within the subset are considered.

    Parameters
    ----------
    sample : string
        the path to the file
    transformation : list of string
        the transformation(s) to apply to the matrix. "OE". "log1p" and "Z-score" are supported. (default Z-score)
    col : int
        optional, the column conting the score to consider. The first column is column 1. (default 4)
    
    Returns
    -------
    a numpy array that is the vectorized (subset of the) file
    the mean of the data, if converted to a Z-score (None otherwise)
    the variance of the data, if converted to a Z-score (None otherwise)
    
    Raises
    ------
    EOFError if the subset file is specified but empty.

    '''
    bed = pd.read_csv(sample, sep = "\t", header=None)
    vector = np.array(bed.iloc[:, (col-1)])

    mean = None
    var = None

    if "log10" in transformation:
        vector = np.log1p(vector)
    if "Z-score" in transformation:
        mean = vector.mean()
        var = vector.var()
        vector = (vector - mean) / np.sqrt(var)
        
    return vector, mean, var

def vector_to_matrix_chrom(vector, res, chrom, subset = None, dist = 0, chrom_sizes = "hg38.chrom.sizes"):
    '''Converts a difference or projection vector to a matrix, for a given chromosome
    
    This function reshapes a vector into a matrix, to match the shape of the original inputs and ease visualization.
    Note: A full matrix will be outputted, spawning the full chromosome, regarless of the size of the subset. 
    If a subset has been specified, the cells otside of it will contain zeroes.
    
    Parameters
    ----------
    vector : numpy array
        vector to convert to a matrix
    res : int
        the resolution
    chrom : string
        the chomosome to consider
    subset : string
        optional, the path to the file containing the regions to subset, if any (default None)
    dist : int
        the distance to consider. All interactions beyond that distance will be ignored. If set to 0, all interactions 
        are kept. (default 0)
    chrom_sizes : string
        optional, the path to the file containing the size (in bp) of each chromosome (default "hg38.chrom.sizes")
    
    Returns
    -------
    a np array matrix, representing the difference, per chromosome
    
    Raises
    ------
    FileNotFoundError if eiher the chrom_sizes file is not found
    NameError if the subset fille is specified but does not contain the considered chromosome

    '''
    sizes = pd.read_csv(chrom_sizes, sep = "\t", header = None)
    nbins = math.ceil(sizes[1][sizes[0] == chrom] / res)
    if (not subset is None) and (subset != ""):
        regions = pd.read_csv(subset, sep = "\t")
        regions_chrom = regions[regions.seqnames == chrom]
        if len(regions_chrom) == 0:
            raise NameError('The subset file is not null and the chromosome ' + chrom + ' is not present in the subset file.')
        
        A = np.full((nbins, nbins), np.nan)
        for i in regions_chrom.index:
            start = int(regions_chrom.start[i])
            end = int(regions_chrom.end[i])
            start_i = math.floor(start/res)
            end_i = math.floor(end/res)
            diff = end_i - start_i
            length = (diff*(diff+1))//2
            vector_subset = vector[range(0, length)]
            array = np.full((diff, diff), np.nan)
            array[np.triu_indices_from(array)] = vector_subset
            A[start_i:end_i, start_i:end_i] = array
            vector = np.delete(vector, range(0, length))
        return A
    elif dist > 0:
        D = np.ceil(dist/res) - 1
        length = (nbins*(nbins+1)//2) - ((nbins-D)*((nbins-D)+1)//2)
        vector_subset = vector[range(0, length)]
        vector = np.delete(vector, range(0, length))
        A = diag_to_matrix(vector_subset, D, nbins)
        return A
    else :
        length = (nbins*(nbins+1))//2
        vector_subset = vector[range(0, length)]
        vector = np.delete(vector, range(0, length))
        array = np.full((nbins, nbins), np.nan)
        array[np.triu_indices_from(array)] = vector_subset
        return array


def vector_to_matrix(vector, res, subset = None, dist = 0, chrom_sizes = "hg38.chrom.sizes", chromlist = []):
    '''Converts a difference or projection vector to a matrix
    
    This function reshapes a vector into a matrix, to match the shape of the original inputs and ease visualization.
    Note: A full matrix will be outputted, spawning the full chromosome, regarless of the size of the subset. 
    If a subset has been specified, the cells otside of it will contain zeroes.
    
    Parameters:
    vector : numpy array
        vector to convert to matrices
    res : int
        the resolution
    chrom : string
        the chomosome to consider
    subset : string
        optional, the path to the file containing the regions to subset, if any (default None)
    dist : int
        the distance to consider. All interactions beyond that distance will be ignored. If set to 0, all interactions 
        are kept. (default 0)
    chrom_sizes : string
        optional, the path to the file containing the size (in bp) of each chromosome (default "hg38.chrom.sizes")
    chromlist : list of string
        optional, the names of the chromosomes to consider. This list is ignored if subset is not null. 
        If set to an empty list, all available chromosomes are considered*. (default [])
        *Available chromosomes are determined from the chrom_sizes file, excluding alternative and unknowns

    Returns
    -------
    a dictionary containing one np array matrix, representing the difference, per chromosome (key is the chromosome name)
    
    Raises
    ------
    EOFError if the subset file is specified but empty.

    '''
    sizes = pd.read_csv(chrom_sizes, sep = "\t", header = None)
    matrices_dict = {}
    
    if (not subset is None) and (subset != ""):
        regions = pd.read_csv(subset, sep = "\t")
        if len(regions) == 0:
            raise EOFError('The subset file is empty.')
        all_chroms = list(set(regions.seqnames))
    elif len(chromlist) > 0:
        all_chroms = chromlist
    else:
        to_keep = [("Un" not in i) and ("alt" not in i) and ("rand" not in i) for i in sizes.iloc[:,0]]
        all_chroms = np.array(sizes.iloc[to_keep,0])
    
    for chrom in all_chroms:
        A = vector_to_matrix_chrom(vector, res, chrom, subset, dist, chrom_sizes)
        matrices_dict[chrom] = A
            
    return matrices_dict

###########################




