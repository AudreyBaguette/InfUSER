import numpy as np

class Pixel(object):
    """
    A class used to represent de data of a given element in the Hi-C matrix

    ...

    Attributes
    ----------
    values_interval : tuple of floats
        the interval of possible value that the pixel can take.
    n_children : int
        the number of children of the nodes.
    
    """
    def __init__(self, n_children):
        """
        Parameters
        ----------
        n_children : int
            The number of children of the node
        
        The constructor fills the value_interval with a tuple of np.nan.
        """
        self.n_children = n_children
        self.values_interval = (np.nan, np.nan)


    def set_interval(self, val1, val2):
        """Updates the scores_matrix

        The matrix is updated such that the line corresponding to the value contain the score.

        Parameters
        ----------
        val1 : float
            The first end of the interval.
        val2 : float
            The second end of the interval.

        """

        if val1 > val2:
            self.values_interval = (val2, val1)
        else:
            self.values_interval = (val1, val2)

        
    def get_interval(self):
        """Returns the interval of values.

        Returns
        ------
        numpy array of int, the values corresponding to the origins giving the minimal score stored in scores_matrix
        
        """
        return(self.values_interval)
        

    def compute_distance(self, value):
        """Compute the difference between a given value and the interval stored in the Pixel instance.

        Parameters
        ----------
        value : float
            The value to compare to the interval.

        Returns
        ------
        float, the distance (euclidian) between the interval and the value
        float, the value leading to that minimal distance
        
        """
        if value < self.values_interval[0]:
            return (self.values_interval[0] - value, self.values_interval[0])
        elif value > self.values_interval[1]:
            return (value - self.values_interval[1], self.values_interval[1])
        else:
            return (0, value)

