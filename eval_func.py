#!/usr/bin/python
"""
Created on Mon Nov 29 16:13:27 2021

@author: Amit
"""

from sklearn import metrics
import numpy as np
from sklearn.calibration import calibration_curve
import matplotlib.pyplot as plt
import pandas as pd
from xgboost import XGBClassifier
from sklearn.model_selection import learning_curve

# =============================================================================
# # Predict and report metrics:  
# =============================================================================
def eval_by_outhold(X_test, y_test, model, thr=None):
    """ predict X_test set, print auc, f1, accuracy, logloss & confusion_matrix """
    if thr:
        probs = model.predict_proba(X_test)[:, 1]
        y_pred = [1 if x >= thr else 0 for x in probs]
    else: 
        y_pred = model.predict(X_test)
    #y_pred = model.predict(X_test)
    print("y_final: ", np.unique(y_test, return_counts=True))
    print("roc_auc -  ", metrics.roc_auc_score(y_test, y_pred))
    print("f1 -  ", metrics.f1_score(y_test, y_pred))
    print("accuracy -  ", metrics.accuracy_score(y_test, y_pred))
    print("log_loss -  ", metrics.log_loss(y_test, y_pred))
    print("cm -  ", metrics.confusion_matrix(y_test, y_pred))
    return metrics.f1_score(y_test, y_pred)
    
   
    
   
# =============================================================================
# # plot calibrated/ reliability diagram
# =============================================================================
def calib_plot( X_test, y_test, model, save=None):
    fig, ax = plt.subplots()
    # predict probabilities
    probs = model.predict_proba(X_test)[:,1]
    # reliability diagram
    fop, mpv = calibration_curve(y_test, probs, n_bins=10, strategy= 'quantile') 
    #fraction of positives, mean predicted probability
    # plot perfectly calibrated
    ax.plot([0, 1], [0, 1], linestyle='--')
    # plot model reliability
    ax.plot(mpv, fop, marker='.')
    ax.set(xlabel= 'predicted probability frequency' , ylabel='observed', 
           title= "calibration of predicted probabilities")
    if save:
        plt.savefig(save+".png")
    else:
        plt.show()
    plt.close()



# =============================================================================
# #ROC curve
# =============================================================================
def roc_plot( X_test, y_test, model, save=None):
    """ if save plot is needed, save= name of plot.(str)"""
    lr_probs = model.predict_proba(X_test)[:, 1] #model results as probabilities per class 
    lr_auc = metrics.roc_auc_score(y_test, lr_probs)
    print('XGB: ROC AUC=%.3f' % (lr_auc))
    # calculate roc curve
    lr_fpr, lr_tpr, thresholds = metrics.roc_curve(y_test, lr_probs)
    fig1, ax1 = plt.subplots()
    # plot the roc curve for the model
    ax1.plot([0, 1], [0, 1],  linestyle='--') #no skill/random prediction
    ax1.plot(lr_fpr, lr_tpr, color='green', label='ROC curve (area = %0.4f)' % lr_auc)
    # axis labels
    ax1.set(xlabel= 'False Positive Rate', ylabel='True Positive Rate', 
            title= "XGBoost ROC ", xlim= [0.0, 1.0], ylim = [0.0, 1.0])
    ax1.legend(loc="lower right")
    if save:
        plt.savefig(save+".png")
    else:
        plt.show()
    plt.close()
    return lr_auc, best_treshold(lr_fpr, lr_tpr, thresholds)
    
    

def best_treshold(fpr, tpr, thresh):
    #metric = tpr - fpr #cal  Youden’s J statistic
    #OR:
    metric = np.sqrt(tpr * (1-fpr)) #cal Geometric Mean
    
    ix = np.argmax(metric)
    best_thresh = thresh[ix]
    print('Threshold=%f' % (best_thresh))
    return best_thresh


# =============================================================================
#     split data
# =============================================================================
from sklearn.model_selection import train_test_split as tts

def ttsplit_by_Class_ratio(X, y ,test_ratio =None, test_size= None):
    test_ratio = 1 if test_ratio is None else test_ratio
    test_size = 0.3 if test_size is None else test_size
    pos =  X[y ==1]
    pos_train, pos_test = tts(pos, test_size= test_size, random_state=0)
    neg = X[ y==0 ]
    neg_train, neg_test = tts(neg, test_size= int(len(pos_test)*test_ratio), random_state=0)
    train = np.vstack([
        np.hstack([pos_train, np.ones((len(pos_train), 1))]), 
        np.hstack([neg_train, np.zeros((len(neg_train), 1))]) ])
    test = np.vstack([
        np.hstack([pos_test, np.ones((len(pos_test), 1))]), 
        np.hstack([neg_test, np.zeros((len(neg_test), 1))]) ])
    np.random.shuffle(train) , np.random.shuffle(test)
    return train , test
    
    
    
def unique_rows(a):
    a = np.ascontiguousarray(a)
    unique_a = np.unique(a.view([('', a.dtype)]*a.shape[1]))
    return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))
    
    
    
# =============================================================================
# #  precision_recall curve
# =============================================================================
def PreREC_plot( X_test, y_test, model, save=None):
    """ plot precision_recall curve. thresholds in scale up.
    if save plot is needed, save= name of plot.(str)"""
    lr_probs = model.predict_proba(X_test)[:, 1]  #model results as probabilities minority class 
    precision, recall, thresholds = metrics.precision_recall_curve(y_test, lr_probs)
    lr_auc = metrics.average_precision_score(y_test, lr_probs)
    precision = zeroes(precision)[1:]
    recall = zeroes(recall)[1:]
    fscore = (2 * precision * recall) / (precision + recall) 
    # locate the index of the largest f score
    ix = np.argmax(fscore)
    print('PreREC_plot Threshold=%f, F-Score=%.3f' % (thresholds[ix], fscore[ix])) 
    # calculate curve
    fig, ax = plt.subplots()
    # plot the PR curve for the model
    ax.plot(recall, precision, color='green', label="XGB (AUC-PR = %0.4f)" % lr_auc)
    no_skill = len(y_test[y_test==1]) / len(y_test)
    ax.plot([0, 1], [no_skill,no_skill],  color= 'blue', linestyle='--', label='No Skill') #no skill/random prediction
    bad_skill = len(y_test[y_test==0]) / len(y_test) #predict all negative
    ax.plot([0, 1], [bad_skill ,bad_skill ],  color='red', linestyle='--', label='bad Skill')
    ax.axvline(thresholds[np.argmin(np.abs(precision-recall))], color="k", ls = "--") 
    ax.scatter(recall[ix], precision[ix], marker='o', color='black', label='Best f1') 
    # axis labels
    ax.set(xlabel= 'Recall (TPR)', ylabel='Precision (PPV)', 
            title= "XGBoost precision_recall ", xlim= [0.0, 1.0], ylim = [0.0, 1.0])
    ax.legend(loc="best")
    if save:
        plt.savefig(save+".png")
    else:
        plt.show()
    plt.close()
    return thresholds[ix] #change/ add f1
    
    
def zeroes(a):
    a[a==0] = np.finfo(np.float).eps
    return a  
    
# =============================================================================
# #     Optimal Threshold Tuning
# =============================================================================
def Op_Threshold( model, X, y):
    """ for test set on trained model.    """
    lr_probs = model.predict_proba(X)[:, 1]
    thresholds = np.arange(0, 1, 0.001)
    # evaluate each threshold
    scores = [metrics.f1_score(y, (lr_probs >= t).astype('int')) for t in thresholds]
    # get best threshold
    ix = np.argmax(scores)
    print('OP_Threshold=%.3f, F-Score=%.5f' % (thresholds[ix], scores[ix]))
    return thresholds[ix] , scores[ix]
    
    
# =============================================================================
# plot together.  f1, presicion and recall thresholds
# =============================================================================
def threshold_plot(model, X, y, save=None):
    """ for test set on trained model.    plot f1, recall and precision as function of thresholds."""
    lr_probs = model.predict_proba(X)[:, 1]
    precision, recall, thresholds = metrics.precision_recall_curve(y, lr_probs)
    precision = zeroes(precision)[1:]
    recall = zeroes(recall)[1:]
    fscore = (2 * precision * recall) / (precision + recall)
    ix = np.argmax(fscore)
    fig, ax = plt.subplots()
    ax.plot(thresholds, precision , color='green', label='PRE')
    ax.plot(thresholds, recall , color='blue', label='Recall')
    ax.plot(thresholds, fscore , color='red', label='f_score')
    ax.scatter(thresholds[ix], fscore[ix], marker='o', color='black', label='Best')
    ax.axvline(thresholds[np.argmin(np.abs(precision-recall))], color="k", ls = "--") 
    # axis labels
    ax.set(xlabel= 'Thresholds', ylabel="Score",  
                title= "XGBoost thresholds" , #, xlim= [0.0, 1.0], ylim = [0.0, 1.0])
                xlim= ([0.025, thresholds[np.argmin(abs(precision-recall))] + 0.2]))
    ax.legend(loc="best")
    print('threshold_plot - max f1 Threshold=%.3f, F-Score=%.5f, PR-Threshold=%.3f ' % (thresholds[ix], fscore[ix],thresholds[np.argmin(np.abs(precision-recall))] ))
    if save:
        plt.savefig(save+".png")
    else:
        plt.show()
    plt.close() 
    return thresholds[ix]
    
    
# =============================================================================
# 
# =============================================================================   
def best_thresh(lr_probs, y_test, pos_label):
    # by PR:
    precision, recall, thresholds = metrics.precision_recall_curve(y_test, lr_probs, pos_label=pos_label)
    precision = zeroes(precision)[1:]
    recall = zeroes(recall)[1:]
    best = np.argmin(np.abs(precision-recall)) #thresholds[np.argmin(abs(precision-recall))]
    fscore = (2 * precision * recall) / (precision + recall) 
    th_pr, s_pr = thresholds[best], fscore[best]
    print('PR_Threshold=%.3f, F-Score=%.5f' % (th_pr, s_pr))
    #print("PR_Threshold range [%.3f, %.3f] " %(np.min(thresholds), np.max(thresholds)))
    
    #by ROC:
    fpr, tpr, thresholds = metrics.roc_curve(y_test, lr_probs, pos_label=pos_label)
    best = np.argmin(abs(tpr + fpr - 1))  # opposite of  as np.argmax(tpr - fpr)
    th_roc, s_roc  = thresholds[best], np.max(tpr - fpr)
    print('ROC_Threshold=%.3f,  J statistic =%.5f' % (th_roc, s_roc))
    #print("ROC_Threshold range [%.3f, %.3f] " %(np.min(thresholds), np.max(thresholds)))
    
    #by loop & f1 :
    thresholds = np.arange(0, 1, 0.001)
    fscore = [metrics.f1_score(y_test, (lr_probs >= t).astype('int'), pos_label= pos_label) for t in thresholds]
    best = np.argmax(fscore)
    th_f1, s_f1 = thresholds[best], fscore[best]
    print('f1_Threshold=%.3f, F-Score=%.5f' % (th_f1, s_f1))
    #print("f1_Threshold range [%.3f, %.3f] " %(np.min(thresholds), np.max(thresholds)))
  
    return th_pr, s_pr ,th_roc, s_roc , th_f1, s_f1
    
    
def def_threshold( model , xtest):
    y_pred = model.predict(xtest)
    lr_probs = model.predict_proba(xtest)
    preds = [1 if x >= 0.5 else 0 for x in lr_probs[:, 1]]
    print("threshold stays 0.5 - ", np.all(y_pred == preds))
    class_index = np.argmax(lr_probs, axis=1)
    print("decision by max proba - ", np.all(class_index == y_pred))
    
    

# =============================================================================
# # fun1 
# =============================================================================
def threshold_check(X, model, y, clas) :
    metrics_s = []
    probs = model.predict_proba(X)[:, clas]
    th_pr, s_pr , th_roc, s_roc  , th_f1, s_f1  = best_thresh( probs, y,clas )
    for thr in [th_pr, th_roc, th_f1]:
        y_pred = [1 if x >= thr else 0 for x in probs]
        f = metrics.f1_score(y, y_pred , pos_label =clas)
        auc = metrics.roc_auc_score(y, y_pred)
        ap = metrics.average_precision_score(y, y_pred, pos_label= clas)
        metrics_s.append([f,auc,ap])
    
    score = pd.DataFrame([[th_pr] +metrics_s[0], [th_roc]+metrics_s[1] , [th_f1]+ metrics_s[2]], 
                         index = ['PR', 'ROC', 'f1'], 
                         columns=(['th' ,'f1_score', 'auc_score' , 'AP_score']))
    return probs, score

#prob , th_score = ev_f.threshold_check(X_test, clf, y_test, 1)
# =============================================================================
# # fun2
# =============================================================================
def curve_metric(y, p, pos_label):
    lr_auc = metrics.roc_auc_score(y, p)
    fpr, tpr, _ = metrics.roc_curve(y, p ,pos_label=pos_label )
    prec, rec, _ = metrics.precision_recall_curve(y, p, pos_label=pos_label )
    pr_auc = metrics.average_precision_score(y, p, pos_label=pos_label) 
    prec = zeroes(prec)
    rec = zeroes(rec)
    metric_s = pd.DataFrame([[fpr, tpr, lr_auc, prec, rec , pr_auc]],  
                         columns=(['fpr', 'tpr' ,'roc_auc', 'prec', 'rec' , 'pr_auc']))
    return metric_s
  
  
  
#results = pd.DataFrame()
#s = ev_f.curve_metric(y_test, prob)
#results = pd.concat([results, a], axis=0, ignore_index=True)
# =============================================================================
# # fun3 plot means
# =============================================================================  
def tolerant_mean(arrs):
    lens = [len(i) for i in arrs]
    arr = np.ma.empty((np.max(lens),len(arrs)))
    arr.mask = True
    for idx, l in enumerate(arrs):
        arr[:len(l),idx] = l
    return arr.mean(axis = -1)#, arr.std(axis=-1)

  
def plot_curve(metric_arr, test_ratio ,save= None):
    """  if plot train & test sets, metric_arr sould be list of arrays when test 1st, train 2nd
    test_ratio = float(minority/ N)
    plot ROC & PR & thresholds by prec,REC,F1 for 1/2 datasets.
    """
    metric_arr = [metric_arr] if type(metric_arr) != list else metric_arr
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3 ,figsize = (18,6.5) )
    #fig.tight_layout()
    sets = iter([  "test" , 'train']) # "class 1 ", "class 0 "] ) #
    for cals in metric_arr:
        set = next(sets)
        mean_fpr = tolerant_mean(cals.fpr)#.data
        mean_tpr =  tolerant_mean(cals.tpr)#.data #cals.tpr.mean(axis = 0)
        mean_auc = cals.roc_auc.mean(axis = 0)
        mean_prec = tolerant_mean(cals.prec).data #cals.prec.mean(axis = 0)
        mean_rec = tolerant_mean(cals.rec).data #cals.rec.mean(axis = 0)
        pr_auc = cals.pr_auc.mean(axis = 0)
        mean_f = (2 * mean_prec * mean_rec) / (mean_prec + mean_rec)
        ix = np.argmax(mean_f)
    
        # roc curve:
        ax1.plot(mean_fpr, mean_tpr, label='%s AUC = %0.4f'  %( set , mean_auc))
        # PR curve:
        ax2.plot(mean_rec[1:], mean_prec[1:] , label='%s AUC-PR = %0.4f' %( set, pr_auc))
        # thresholds curve:
        thresholds =  np.arange(0, 1, 1/len(mean_f))[1:]
        ax3.plot(thresholds, mean_prec[1:] , color='green')#, label='PRE')
        ax3.plot(thresholds, mean_rec[1:] , color='blue')#, label='Recall')
        ax3.plot(thresholds, mean_f[1:] , color='red')#, label='f_score')
        ax3.scatter(thresholds[ix], mean_f[ix], marker='o' , label=set+' Best f1')
        ax3.axvline(thresholds[np.argmin(np.abs(mean_prec-mean_rec))], color="k", ls = "--")
        print(set+ ": threshold by PR = %.4f , threshold by max f1 = %.4f (%.4f)"%( thresholds[np.argmin(np.abs(mean_prec-mean_rec))], thresholds[ix], mean_f[ix]))
         
    
    # roc curve:  
    ax1.plot([0, 1], [0, 1],  linestyle='--', label='No Skill')#no skill/random prediction
    ax1.set( xlabel= 'False Positive Rate', ylabel='True Positive Rate' ,title = 'ROC curve' ,xlim= [0.0, 1.0], ylim = [0.0, 1.0])
    ax1.legend(loc='upper right')
    # PR curve:
    ax2.plot([0, 1], [test_ratio,test_ratio], linestyle='--', label='No Skill (all 1)') #no skill/predict on;y true
    ax2.plot([0, 1], [1-test_ratio, 1-test_ratio], linestyle='--', label='No Skill (all 0)') #no skill/predict on;y false
    ax2.set( xlabel= 'Recall (TPR)', ylabel= 'Precision (PPV)', title= 'PR curve' , xlim= [0.0, 1.0], ylim = [0.0, 1.0])
    ax2.legend(loc = 'upper right')
    # thresholds curve:
    ax3.set(xlabel= 'Thresholds', ylabel="Score",  title= "XGBoost metrics by thresholds "  ,xlim= [0.0, 1.0], ylim = [0.0, 1.0])
    ax3.legend(["Precision", "Recall", 'F-score'] , loc='upper right')  
    
    if save:
        plt.savefig(save+".png")
    else:
        plt.show()
    plt.close()



# =============================================================================
# #for learning curve
# =============================================================================

def plot_validation_curve(train_scores_arr, test_scores_arr, n_estimators_range, save =None ):
    score = save.split("_")[0].upper()
    fig = plt.figure(figsize=(10, 6), dpi=100)
    for i in range(train_scores_arr.shape[2]):
        train_scores = train_scores_arr[:,:,i]
        test_scores = test_scores_arr[:,:,i]
        train_scores_mean = np.mean(train_scores, axis=1)
        train_scores_std = np.std(train_scores, axis=1)
        test_scores_mean = np.mean(test_scores, axis=1)
        test_scores_std = np.std(test_scores, axis=1)

    
        plt.plot(n_estimators_range, train_scores_mean, color="r")
        plt.plot(n_estimators_range, test_scores_mean, color="g")

        plt.fill_between(n_estimators_range, train_scores_mean - train_scores_std,
                         train_scores_mean + train_scores_std, alpha=0.2, color="r")
        plt.fill_between(n_estimators_range, test_scores_mean - test_scores_std,
                         test_scores_mean + test_scores_std, alpha=0.2, color="g")
        i = np.argmax(test_scores_mean)
        print("Best cross-validation {0} ({1:.2f}) obtained for {2} trees".format(score,test_scores_mean[i], n_estimators_range[i]))
    

    plt.title("Validation Curve with XGBoost")
    plt.xlabel("Number of trees")
    plt.ylabel(score)
    plt.ylim(0.0, 1.1)
    plt.axhline(y=1, color='k', ls='dashed')
    plt.legend(["Training score", "Cross-validation score"], loc="best")
    if save:
        plt.savefig(save+".png")
    else:
        plt.show()
    plt.close()

#ev_f.plot_curve(results, 0.1  ,"nestedCV_plots_k%s_%s_R%s_v5.png" %(k_len,eval_metric ,test_ratio) )



def plot_learning_curve( model , X_train, y_train, cv, alpha=0.1, save=None, scoring='f1'):
    train_sizes, train_scores, test_scores = learning_curve(
        estimator=model , X=X_train, y=y_train,
        train_sizes=np.arange(0.01, 1.0, 0.1), cv=cv, scoring=scoring, n_jobs= -1)
    train_mean = np.mean(train_scores, axis=1)
    train_std = np.std(train_scores, axis=1)
    test_mean = np.mean(test_scores, axis=1)
    test_std = np.std(test_scores, axis=1)
    plt.plot(train_sizes, train_mean, label='Training', color='blue', marker='o')
    plt.fill_between(train_sizes, train_mean + train_std,
                     train_mean - train_std, color='blue', alpha=alpha)
    plt.plot(train_sizes, test_mean, label='Validation', color='red', marker='o')

    plt.fill_between(train_sizes, test_mean + test_std, test_mean - test_std, color='red', alpha=alpha)
    plt.title('Learning curve '+save)
    plt.xlabel('Training set size')
    if scoring == 'f1':
        plt.xlim([0,1])
        plt.ylim([0,1])
    plt.ylabel('Score')
    plt.grid(ls='--')
    plt.legend(loc='best')
    if save:
        plt.savefig(save+"_learning curve.png")
    else:
        plt.show()
    plt.close()
                      


# =============================================================================
# nested CV eval plots.  PR & ROC curves
# =============================================================================

def score_probs(y_true, y_pred):
    """ fun takes y_true, y_pred ( shuld be probabilities)
    #returns auc, fpr, tpr, _ , prec, rec, thresholds, pr_auc
    for ROC & PR & thresholds plot
    from each nestedCV loop.
    for train & test sets."""
    roc_auc = metrics.roc_auc_score(y_true, y_pred) #float per cv
    pr_auc = metrics.average_precision_score(y_true, y_pred) #float per cv
    #f = metrics.f1_score(y_true, y_pred) #float per cv
    
    fpr, tpr = metrics.roc_curve(y_true, y_pred)[:-1]# (fpr, tpr)

    pr_ls = metrics.precision_recall_curve(y_true, y_pred) # tuple(precision, recall, thresholds)
    prec, rec = [arr[:-1] for arr in pr_ls[:-1]]
    thresholds = pr_ls[-1]
    return roc_auc, pr_auc,  fpr, tpr, prec, rec, thresholds #f,
    


def plots_nested_cv(train, test, test_ratio, save=None):
    """ each is pd.DataFrame() with columns ['auc' , 'f1', 'ap',  'roc_array', 'pr_array']
    where the firsts are floats and lasts are arrays.
    each row is CV fold.
    this function plot ROC, PR & threshold by metrics.
    after nestedcv
    """
    fig, ((ax1, ax2), (ax3,ax4)) = plt.subplots(2, 2 ,figsize = (14,14) )
    color = ["blue", "orange"]
    sets = [train, test]
    for i in range(2): #set & color
        for l in range(len(sets[i])): #cv within set
            ax1.plot(sets[i].loc[l].fpr, sets[i].loc[l].tpr , color=color[i]) #ROC curve:
            ax2.plot(sets[i].loc[l].rec, sets[i].loc[l].prec , color=color[i]) #PR curve:
            #threshold by metrics plot:
            prec = zeroes(sets[i].loc[l].prec)
            rec = zeroes(sets[i].loc[l].rec)
            fscore = (2 * prec * rec) / (prec + rec)
            ax_th= ax3 if i==0 else ax4
            ax_th.plot(sets[i].loc[l].thresholds, prec , color='green' , linestyle='-'*(i+1))
            ax_th.plot(sets[i].loc[l].thresholds, rec , color='blue', ls='-'*(i+1))
            ax_th.plot(sets[i].loc[l].thresholds, fscore , color='red', ls='-'*(i+1))
            ax_th.axvline(sets[i].loc[l].best, color="k", ls = (0, (1, 10))) #best threshold of this cv
 
    ax3.set(xlabel= 'Thresholds', ylabel="Score",  title= "XGBoost metrics by thresholds - Train"  ,xlim= [0.0, 1.0], ylim = [0.0, 1.0])
    ax3.legend(["Precision", "Recall", 'F-score', 'CVbest_threshold'] , loc="best")  
    ax4.set(xlabel= 'Thresholds', ylabel="Score",  title= "XGBoost metrics by thresholds - Test (hold-out)"  ,xlim= [0.0, 1.0], ylim = [0.0, 1.0])
    ax4.legend(["Precision", "Recall", 'F-score', 'CVbest_threshold'] , loc="best")
    #ROC curve:
    ax1.plot([0, 1], [0, 1],  linestyle='--', label='No Skill')#no skill/random prediction
    ax1.set( xlabel= 'False Positive Rate', ylabel='True Positive Rate' ,
            title = 'ROC curve' ,xlim= [0.0, 1.0], ylim = [0.0, 1.0])
    ax1.legend(['Train mean AUC = %0.4f(%0.4f)' % (train.auc.mean(), train.auc.std()),
                'Test mean AUC = %0.4f(%0.4f)' % (test.auc.mean(), test.auc.std()) ], loc="best")
    leg = ax1.get_legend()
    leg.legendHandles[0].set_color('blue')
    leg.legendHandles[1].set_color('orange')
    #PR curve:
    no_skill = 1/ (test_ratio+1) #len(y_test[y_test==1]) / len(y_test) #predict all positive
    ax2.plot([0, 1], [no_skill,no_skill], linestyle='--', 
             label='No Skill (all 1)') #no skill/predict on;y true
    ax2.set( xlabel= 'Recall (TPR)', ylabel= 'Precision (PPV)', 
            title= 'PR curve' , xlim= [0.0, 1.0], ylim = [0.0, 1.0])
    ax2.legend(['Train mean PR-AUC = %0.4f (fscore = %0.4f-+%0.4f)' % 
                (train.ap.mean(),train.f1.mean(), train.f1.std()),
                'Test mean PR-AUC = %0.4f (fscore = %0.4f-+%0.4f)' % 
                (test.ap.mean(),test.f1.mean(), test.f1.std())], loc = "best")
    leg = ax2.get_legend()
    leg.legendHandles[0].set_color('blue')
    leg.legendHandles[1].set_color('orange')
    if save:
        plt.savefig(save+".png")
    else:
        plt.show()
    plt.close()

  
# =============================================================================
# get set ratio and size
# =============================================================================
def set_ratio( lable):
    counts = np.unique(lable, return_counts=True)[1]
    if len(counts) != 2:
        return counts[0] , 0.000
    neg, pos = np.unique(lable, return_counts=True)[1]
    ratio = neg/pos
    size = pos+ neg
    return size , ratio

# =============================================================================
# limits of undersample stratgy // threshold (?)			
# =============================================================================
def min_undersample_stratgy(y,n_splits = 5 ):
    TS = len(y) # full train set size =	21900
    TP = len(y[y==1]) # total pos size = 2034
    TN = len(y[y==0]) # total neg size = 19866

    Ste = TS/n_splits # test set size = 4380 (if n_splits = 5)
    Pte = 1 #range( 1, TP*0.9 ) #test pos size : max-min :  TP*0.9	- 1 : 1830 - 1                  #changed 16/5
    Nte = Ste - Pte  #test neg size : max-min : Ste-TP-1 - Ste-TP*0.9	 : 4379 - 2549
    
    Str = TS - Ste # train set size = 	17520
    Ptr = TP - Pte # train pos size:	min-max :427	2032: TP*0.1 - TP-1: 203- 2033
    Ntr = Str - Ptr # train neg size:	min-max :17093	15488: Str-TP-1 -  : 15487 - 17317
    
    return Ptr/Ntr #maximal minimum for sampling_strategy =  pos/neg



# =============================================================================
# cv scorer
# =============================================================================

def f1_ratio_scorer(clf, X, y):
    ''' calculate ratio of f1 to no-skill f1
    for cv scoring.
    '''
    y_pred = clf.predict(X)
    size = len(y)
    f1_score = metrics.f1_score(y, y_pred )
    no_skill_f1 = metrics.f1_score(y, np.ones( size ) )
    
    return f1_score/no_skill_f1



from collections import Counter

def _true_pos_size_scorer(clf, X, y):
    return len(y[y==1])

            
def _true_ratio_scorer(clf, X, y):
    # ratio = neg/pos. 
    target_stats = Counter(y)
    pos_size = target_stats[1]
    neg_size = target_stats[0]        
    return neg_size/pos_size #ev_f.set_ratio(y)[1] same...
