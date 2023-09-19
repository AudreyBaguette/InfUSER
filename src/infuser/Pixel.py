import numpy as np

class Pixel(object):
    """
    A class used to represent de data of a given element in the Hi-C matrix

    ...

    Attributes
    ----------
    n_values : int
        the number of possible values that each pixel can take (default 9).
        the values will range from 0 to n_values-1.
    n_children : int
        the number of children of the nodes.
    scores_matrix : numpy array of floats
        an array storing the score of each possible value. It has n_values x 1 dimensions.
    origins_matrix : numpy array of integers
        an array storing the value of the children producing the scores in scores_matrix. It has n_values x n_children dimensions.
    
    """
    def __init__(self, n_children, n_values = 9):
        """
        Parameters
        ----------
        n_children : int
            The number of children of the node
        n_values : int
            The number of possible values the pixel can take (default is 9)
        
        The constructor fills scores_matrix with np.inf and origins_matrix with np.nan.
        If the node does not have children, the matrix is null
        """
        self.n_children = n_children
        self.n_values = n_values
        self.scores_matrix = np.full((n_values, 1), np.inf)
        if n_children > 0:
            self.origins_matrix = np.full((n_values, n_children), np.nan)
        else :
            self.origins_matrix = None

    def set_score(self, value, score):
        """Updates the scores_matrix

        The matrix is updated such that the line corresponding to the value contain the score.

        Parameters
        ----------
        value : int
            The value for which the score and origins should be updated.
        score : float
            The new score to input in scores_matrix.

        Raises
        ------
        ValueError
            If value is not between 0 and self.n_values.
        """

        if value >= 0 and value < self.n_values:
            self.scores_matrix[value] = score
        else:
            raise ValueError("The value parameter is not within expected values")


    def set_score_and_origins(self, value, score, origins):
        """Updates the scores_matrix and origins_matrix

        The matrices are updated such that the line corresponding to the value contain the score and origins.

        Parameters
        ----------
        value : int
            The value for which the score and origins should be updated.
        score : float
            The new score to input in scores_matrix.
        origins : list of int
            The values from which the score was computed.

        Raises
        ------
        ValueError
            If value is not between 0 and self.n_values.
        """

        assert np.all([isinstance(i, int) for i in origins])
        if value >= 0 and value < self.n_values:
            self.scores_matrix[value] = score
            self.origins_matrix[value] = origins
        else:
            raise ValueError("The value parameter is not within expected values")
        

    def get_min_value(self):
        """Returns the value of resulting in the minimal score

        Returns
        ------
        int, the value corresponding to the minimal score stored in scores_matrix
        float, the score of that value
        """

        min = np.where(self.scores_matrix == np.min(self.scores_matrix))[0]
        # If there are more than one min, select one randomly
        if len(min) > 1:
            select = np.random.randint(len(min))
            min_value = min[select]
        else:
            min_value = min[0]
        return(min_value, np.min(self.scores_matrix))

        
    def get_min_origins(self, value):
        """Returns the value of the children resulting in the value passed as parameter.

        Parameters
        ----------
        value : int
            The value for which the origins' values are retrieved.
        
        Returns
        ------
        numpy array of int, the values corresponding to the origins giving the minimal score stored in scores_matrix
        
        Raises
        ------
        ValueError
            If value is not between 0 and self.n_values.
        """
        if value >= 0 and value < self.n_values:
            min_origin = self.origins_matrix[value]
            return(min_origin.astype(int))
        else:
            raise ValueError("The value parameter is not within expected values")
        

        
    def reset_scores(self):
        """
        Erase the data within scores_matrix and reset the values.
        """
        self.scores_matrix = np.full((self.n_values, 1), np.inf)

