import numpy as np
from .utils import kth_diag_indices
import copy

def to_indices(vector, step, min_score, max_score):
    '''Convert a Z-scores vector to an indices vector
    
    This function takes a vector containing Z-scores and transforms it to indices.
    
    Parameters
    ----------
    vector : numpy array
        the data to convert
    step : float
        the difference in Z-score between two consecutive indices
    min_score : float
        the minimal value to consider
    max_score : float
        the maximal value to consider
    
    Returns
    -------
    a numpy array containing the converted vector

    '''
    indices = np.zeros(len(vector))
    max_index = (max_score - min_score) / step
    for i in range(0, len(vector)):
        # Convert the scores to indices
        index = (vector[i] - min_score) / step
        # Put everything below 0 to 0 and everything above the max to the max
        if index < 0:
            index = 0
        elif index > max_index:
            index = max_index
        elif np.isnan(index):
            index = 0
        indices[i] = int(index)
    return np.array(indices)

def to_Zscore(vector, step, min_score):
    '''Convert an indices vector to a Z-scores vector
    
    This function takes a vector containing indices and transforms it to Z-scores.
    
    Parameters
    ----------
    vector : numpy array
        the data to convert
    step : float
        the difference in Z-score between two consecutive indices
    min_score : float
        the minimal value to consider
    
    Returns
    -------
    a numpy array containing the converted vector

    '''
    new_vec = (vector * step) + min_score
    return new_vec

def to_HiC_scores(input_matrix, chrom, samples, transformation, means_dict, vars_dict):
    '''Convert a Z-scores matrix to an Hi-C matrix
    
    This function takes a matrix containing Z-scores and transforms it to Hi-C like scores.
    
    Parameters
    ----------
    input_matrix : numpy array
        the data to convert
    chrom : string
        the chromosome represented by the matrix
    samples : array of strings
        the name of samples from which the mean and variance should be considered
    transformation : list of string
        the transformation(s) to apply to the matrix. "log1p" and "Z-score" are supported. (default Z-score)
    means_dict : dictionary
        the diagonal means, per chromosome and per sample
    vars_dict : dictionary
        the diagonal variances, per chromosome and per sample
    
    Returns
    -------
    a numpy array containing the converted matrix

    '''
    matrix = copy.deepcopy(input_matrix)
    if "log1p" in transformation:
        matrix = np.e**matrix - 1
    
    if ("Z-score" in transformation) or ("z-score" in transformation) or ("zscore" in transformation) or ("Zscore" in transformation):
        means_to_consider = means_dict[samples[0]][chrom]
        vars_to_consider = vars_dict[samples[0]][chrom]
        if len(samples) > 1:
            for sample in range(1, len(samples)):
                means_to_consider = np.append(means_to_consider, means_dict[samples[sample]][chrom])
                vars_to_consider = np.append(vars_to_consider, vars_dict[samples[sample]][chrom])
        means_to_consider = means_to_consider.reshape(len(samples), len(means_dict[samples[0]][chrom]))
        vars_to_consider = vars_to_consider.reshape(len(samples), len(means_dict[samples[0]][chrom]))
        
        for diag in range(0, len(matrix)):
            mean = np.mean(means_to_consider[:,diag])
            var = np.mean(vars_to_consider[:,diag])
            matrix[kth_diag_indices(matrix, diag)] = matrix.diagonal(diag) * np.sqrt(var) + mean
            matrix[kth_diag_indices(matrix, -diag)] = matrix.diagonal(-diag) * np.sqrt(var) + mean
    
    return matrix

def to_scores(input_vector, samples, transformation, means_dict, vars_dict):
    '''Convert a Z-scores vector to a vector with "regular" scores
    
    This function takes a vector containing Z-scores and transforms it to bed scores.
    
    Parameters
    ----------
    input_vector : numpy array
        the data to convert
    samples : array of strings
        the name of samples from which the mean and variance should be considered
    transformation : list of string
        the transformation(s) to apply to the matrix. "log1p" and "Z-score" are supported. (default Z-score)
    means_dict : dictionary
        the diagonal means, per chromosome and per sample
    vars_dict : dictionary
        the diagonal variances, per chromosome and per sample
    
    Returns
    -------
    a numpy array containing the converted vector
    
    '''
    vector = copy.deepcopy(input_vector)

    if "log1p" in transformation:
        vector = np.e**vector - 1
        
    if ("Z-score" in transformation) or ("z-score" in transformation) or ("zscore" in transformation) or ("Zscore" in transformation):
        means_to_consider = means_dict[samples[0]]
        vars_to_consider = vars_dict[samples[0]]
        if len(samples) > 1:
            for sample in range(1, len(samples)):
                means_to_consider = np.append(means_to_consider, means_dict[samples[sample]])
                vars_to_consider = np.append(vars_to_consider, vars_dict[samples[sample]])
        #means_to_consider = means_to_consider.reshape(len(samples), len(means_dict[samples[0]]))
        #vars_to_consider = vars_to_consider.reshape(len(samples), len(means_dict[samples[0]]))
        
        mean = np.mean(means_to_consider)
        var = np.mean(vars_to_consider)
        vector = vector * np.sqrt(var) + mean
    
    return vector
