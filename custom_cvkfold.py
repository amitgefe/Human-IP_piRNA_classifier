#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Thu May 12 15:58:28 2022

@author: amitgefe
"""


import numpy as np
from sklearn.model_selection._split import _build_repr

class custom_Repeatedkfold():
    
    def __init__(self, n_repeats=10, n_splits=5,):
        self.n_repeats = int(n_repeats)
        self.n_splits = n_splits
        self.custom_k = custom_kfold(self.n_splits)
    
    def split(self, X, y=None, groups=None):
        """Generates indices to split data into training and test set.
        """
        n_repeats = self.n_repeats

        for idx in range(n_repeats):
            cv = self.custom_k
            for train_index, test_index in cv.split(X, y, groups):
                yield train_index, test_index
    
    def get_n_splits(self, X=None, y=None, groups=None):
        """Returns the number of splitting iterations in the cross-validator    """
        cv = self.custom_k
        return cv.get_n_splits(X, y, groups) * self.n_repeats
    
    def __repr__(self):
        return _build_repr(self)
    
    
class custom_kfold():
    def __init__(self, n_splits=5, *, shuffle=False, random_state=None):
        self.n_splits = int(n_splits)
        self.shuffle = shuffle
        self.random_state = random_state
        
    def split(self, X, y, groups=None):
        """Generate indices to split data into training and test set.
        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Training data, where `n_samples` is the number of samples
            and `n_features` is the number of features.
        y : array-like of shape (n_samples,)
            The target variable for supervised learning problems.
      
        Yields
        ------
        train : ndarray
            The training set indices for that split.
        test : ndarray
            The testing set indices for that split.
        """
        n_samples = len(X)
        if self.n_splits > n_samples :
            raise ValueError(
                ("Cannot have number of splits n_splits={0} greater"
                    " than the number of samples: n_samples={1}.").format(self.n_splits, n_samples)
            )
        elif n_samples//self.n_splits <= y[y==1].shape[0] :
            raise ValueError(
                ("Cannot have fold size {2} from n_splits={0} smaller"
                    " than the number of total positive samples ={1}."
                    ).format(self.n_splits, y[y==1].shape[0], n_samples//self.n_splits)
            )
            

        indices = np.arange(len(X))
        for test_index in self._iter_test_masks(X, y, groups):
            train_index = indices[np.logical_not(test_index)]
            test_index = indices[test_index]
            yield train_index, test_index
    
    def _iter_test_masks(self, X, y, groups=None):
        """Generates boolean masks corresponding to test sets.
        By default, delegates to _iter_test_indices(X, y, groups)
        """
        for test_index in self._iter_test_indices(X, y, groups):
            test_mask = np.zeros(len(X), dtype=bool)
            test_mask[test_index] = True
            yield test_mask
    
            
    def _pos_prob(self, n_splits):
        """ generates n_splits probabilities for positive set folds.
        """
        random = np.random.random_sample(n_splits-1)
        while random.sum() >= 1 :
            random = np.random.random_sample(n_splits-1)
            continue
        return np.append(random , 1- random.sum() )
    
    
    def _iter_test_indices(self, X, y, groups=None):
        n_samples = len(X)
        n_splits = self.n_splits
        indices = np.arange(n_samples)           
        #if self.shuffle:
         #   check_random_state(self.random_state).shuffle(indices)
            
        fold_sizes = np.full(n_splits, n_samples // n_splits, dtype=int)
        fold_sizes[: n_samples % n_splits] += 1
        probs = self._pos_prob(n_splits)
        fold_sizes = np.vstack((fold_sizes,probs)).T
        
        for fold_size, p in fold_sizes[:-1]:
            test_indices = self._test_ind_random_dis(indices, p, fold_size, y)
            indices[test_indices] = False
            yield test_indices 
        yield np.array( [i for i in set(indices)])


    def _test_ind_random_dis(self, indices, p, fold_size , y): 
        try:
            fold_size = int(fold_size)
            pos_indices = np.array([i for i in set(indices) if y[i]==1 and i != 0] )
            pos_size = int( round( p * y[y==1].shape[0] ) )
            if  pos_size >= fold_size:
                raise ValueError(" pos_size >= fold_size.  ")
            
            pos_test_index = np.random.choice(pos_indices, pos_size, replace=False)
            neg_test_index = np.random.choice(  [i for i in set(indices) if y[i]==0 and i != 0],  fold_size-pos_size, replace=False)  
            test_index = np.concatenate( (pos_test_index, neg_test_index) , axis = 0 )
            return test_index
        
        except Exception as err:
            trace = []
            tb = err.__traceback__
            while tb is not None:
                    trace.append({  "filename": tb.tb_frame.f_code.co_filename,
                                  "name": tb.tb_frame.f_code.co_name,
                                  "lineno": tb.tb_lineno
                                  })
                    tb = tb.tb_next
            print(str({ 'type': type(err).__name__,
                        'message': str(err),
                        'trace': trace
                                        }))
            
            print("pos_size =",  pos_size,  " type(pos_size) = " ,  type(pos_size),
                  "\nfold_size = " ,  fold_size, " type(fold_size) = " ,  type(fold_size), 
                   "\nlen( pos_indices) = " , len( pos_indices), " type(pos_indices) = " ,  type(pos_indices),
                  "\n", err.__class__, err )
   
    
    def get_n_splits(self, X=None, y=None, groups=None):
        """Returns the number of splitting iterations in the cross-validator """
        return self.n_splits



