#!/usr/bin/python
## final training
#qsub -cwd -b y -q alal.q -pe ompi 8 -N XGBm_train  /gpfs0/alal/projects/utils/python3.8  XGBm_train.py 
import numpy as np
from xgboost import XGBClassifier
from datetime import datetime
import pickle

# =============================================================================
#     0. load data
# =============================================================================
in_file = "/gpfs0/alal/users/amitgefe/IP_classifier/train_chm13v2.0/final_set/k3-mers_dist_all.csv"
print("in_file = " ,in_file)
data = np.genfromtxt( in_file ,delimiter=',')
print( "data.shape : ", data.shape )
data =  np.unique(data , axis =0)
print( "original data set (no dup) : ",(np.unique(data[:,-1], return_counts=True)))
X_train, y_train = data[:,:-1], data[:,-1].astype("int32")
print( "full train set ", np.unique(y_train, return_counts=True) ) 

# =============================================================================
#       1. define the model 
# =============================================================================
#model 1 - SPW
SPW = XGBClassifier(use_label_encoder=False,  n_jobs = -1 ,verbosity =0, 
                    scale_pos_weight= 12, n_estimators= 330, 
                    colsample_bytree= 0.88, learning_rate= 0.062, 
                    max_depth= 9, min_child_weight= 1, subsample= 0.71)
                    
print("SPW fit start ", datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
SPW.fit( X_train, y_train )
print("SPW fit finish ", datetime.now().strftime("%d/%m/%Y %H:%M:%S"))

# save the model to disk
filename = 'XGB_model.sav'
pickle.dump(SPW , open(filename, 'wb'))
print("\n\nDONE")

