#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 10:39:05 2019

@author: tanerarslan
"""

import os
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import StratifiedKFold
from sklearn.utils import safe_sqr
from collections import Counter

def load_data(path, file_name):
    
    '''
    Load the data
    
    Parameter
    --------
    path: string, path of the data
    
    file_name: string, name of the file
    
    Return
    -------
    df: data.frame
    '''
    
    df_path = os.path.join(path, file_name)
    df = pd.read_table(df_path, sep = '\t')
    
    return df


def feature_sel_CV(step, replicate):
    '''
    Select the most robust protein markers to separate the clusters
    
    Parameter:
    ----------
        step: corresponds to the (integer) number of features to remove at each iteration, (default= 5)
        replicate : numbers of iteration
    
    Return:
    -------
        result_array = 203 X replicate, proteins as array
        counts = frequency of proteins in the list
        train_score = list of each models train score
        test_score = list of each models test score
        
    '''
    
    #define empty lists and numpy arrays
    result_array = np.array([])
    train_score = []
    test_score = []
    pred = {'prediction': [], 'observation': [], 'sample_id':[], 'iteration': []}
    pred_df = pd.DataFrame(data=pred)
    
    for i in range(replicate):
    
        ### main data split into 80% train and 20% validation ##
        train, test = train_test_split(data, test_size=0.2, random_state=i)
            
        clstrs = train.loc[ : ,"Cls"].values
        
        # delete cluster column
        train = train.drop("Cls", axis=1)
        train = train.to_numpy()
    
        #protein names
        features = data.columns[0:1549].tolist()
    
        j = 0
        while  train.shape[1] > 201:
            j += 1
            svc = SVC(kernel="linear")
            Cs = np.array([0.5, 1.0, 10,100])
            
            # get the hyperparameters
            clf = GridSearchCV(estimator=svc, 
                               param_grid=dict(C=Cs), 
                               cv = 5, 
                               return_train_score=True)
            
            clf.fit(train, clstrs)
            
            # do 5-fold cross validation
            cv_test_error = []
            skf = StratifiedKFold(n_splits=5, random_state=j)
            for trn, tst in skf.split(train, clstrs):
                train_train, train_test = train[trn], train[tst]
                train_clstrs, test_clstrs = clstrs[trn], clstrs[tst]
                val_clf = SVC(C = list(clf.best_params_.values())[0], kernel="linear")
                val_clf.fit(train_train, train_clstrs)
                cv_test_error.append(val_clf.score(train_test, test_clstrs))
            mean_cv_test_error = np.array(cv_test_error).mean()
        
            ## train classification for RFE
            
            rfe_clf = SVC(C = list(clf.best_params_.values())[0], kernel="linear")
            rfe_clf.fit(train, clstrs)
            
            # get coeffs
            coefs = rfe_clf.coef_
            
            # get ranks
            if coefs.ndim > 1:
                ranks = np.argsort(safe_sqr(coefs).sum(axis=0))
            else:
                ranks = np.argsort(safe_sqr(coefs))
        
            #remove the X least important features from the array
            to_remove_index = []
            
            for r in range(step):
                to_remove_index.append(ranks[r])
            to_remove_index.sort(reverse = True)
        
            #remove from largest index to smallest
            for f in to_remove_index:
                train = np.delete(train, f, axis = 1)
                del features[f]
            
    
        ##train and test data with the final features
    
        svc.model = SVC(kernel="linear")
        paramC = np.array([0.5, 1.0, 10,100])
    
        # get the hyperparameters
        final_clf = GridSearchCV(estimator=svc.model, 
                                 param_grid=dict(C=paramC), 
                                 cv = 5, 
                                 return_train_score=True)
    
        final_clf.fit(train, clstrs)
    
        ##fix the test set
        test_clstrs = test.loc[ : ,"Cls"].values

        #delete cluster column
        test = test.drop("Cls", axis=1)
        test = test.loc[ : , features]

        #test sample names
        samples = test.index.values

        #convert data into array
        test_data = test.to_numpy()
        
        #get the outputs
        train_score.append(final_clf.best_score_)
        test_score.append(final_clf.score(test_data, test_clstrs))
        
        result_array = np.append(result_array, np.array(features), axis = 0)
        counts = Counter(result_array)

        iterate_pred = {'prediction':final_clf.predict(test_data), 'observation' : test_clstrs, 'sample_id': samples ,'iteration': i}
        iterate_df = pd.DataFrame(data = iterate_pred)
        pred_df = pd.concat([pred_df, iterate_df], ignore_index=True)
        
    return (result_array, counts, train_score, test_score, pred_df)

def dict2txt(counts, fileName):
    '''
    Export dictionary to txt file

    Parameter:
    ---------
    counts: Frequency of the features

    fileName: path

    return:
    -------
    txt file
    '''
    path = fileName + ".txt"
    with open(path, 'w') as f:
        for key, value in counts.items():
            string = "{}\t{}".format(key, value)
            f.write("%s\n" % string)
    return

def list2txt(lsts, fileName):
    '''
    Export list to txt file

    Parameter:
    ---------
    lsts: scores for the training and testing

    fileName: path

    return:
    -------
    txt file
    '''
    path = fileName + ".txt"
    with open(path, 'w') as f:
        for score in lsts:
            f.write("%s\n" % score)
    return


#Proteomics data
data = load_data(path= "./data", file_name= "LC_data_6_Cls.txt")
#Differentially abundant proteins
deqms_features = load_data(path= "./data", file_name= "Most_pvalue_UpDown_sig_prot.txt")        
deqms_feature_list = deqms_features["Proteins"].tolist()
deqms_feature_list.append("Cls")
data = data.loc[ : , deqms_feature_list]  


result_array, counts, train_score, test_score, pred_df = feature_sel_CV(step=1, replicate=100)

#export results
pred_df.to_csv("preds.txt", header=True, index=True, sep='\t', mode='a')
np.savetxt("result_array.txt", result_array, fmt='%s')
dict2txt(counts, "feature_frequency")
list2txt(test_score, "./data")
list2txt(train_score, "./data")   

        
        
