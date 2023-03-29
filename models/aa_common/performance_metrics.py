import sys
sys.path.append("../variant_effect_analysis")

import numpy as np
import pandas as pd
from sklearn import metrics as sklearn_metrics

def get_auc_roc_score(non_nan_result_df:pd.DataFrame):
    df = non_nan_result_df.copy(deep=True)
    auc_roc_score = sklearn_metrics.roc_auc_score(df["class_numeric"], df["pred"])
    if auc_roc_score < 0.5: 
        auc_roc_score = 1 - auc_roc_score
    print(f"\tAUC-ROC: {auc_roc_score:.3f}")
    return auc_roc_score

def get_auc_pr_score(non_nan_result_df:pd.DataFrame):
    df = non_nan_result_df.copy(deep=True)
    precisions, recalls, thresholds = sklearn_metrics.precision_recall_curve(df["class_numeric"], df["pred"])
    auc_pr_score = sklearn_metrics.auc(recalls, precisions)

    if auc_pr_score < 0.5: 
        auc_pr_score = 1 - auc_pr_score
    print(f"\tAUC-PR: {auc_pr_score:.3f}")
    return auc_pr_score, precisions, recalls, thresholds


def get_f1max_and_th(precisions, recalls, thresholds):
    zero_indices = [i for i in range(precisions.shape[0]) if precisions[i]==0. and recalls[i]==0.]
    
    # if precision and recall both are 0, f1=nan
    # so removing those entries where both precision and recall are 0.
    precisions = np.delete(precisions, zero_indices)
    recalls = np.delete(recalls, zero_indices)
    thresholds = np.delete(thresholds, zero_indices)
    
    f1_scores = (2*recalls*precisions)/(recalls+precisions) 
    th_max = thresholds[np.argmax(f1_scores)]
    f1_max = np.max(f1_scores)
    print(f"\tBest F1-Score: {f1_max:.3f} at threshold: {th_max:.3f}")
    return f1_max, th_max

def get_precision_score(non_nan_result_df:pd.DataFrame, th):
    df = non_nan_result_df.copy(deep=True)
    df.loc[df["pred"]>=th, "class_numeric_pred"] = 0
    df.loc[df["pred"]<th, "class_numeric_pred"] = 1
    score = sklearn_metrics.precision_score(df["class_numeric"], df["class_numeric_pred"], zero_division=0)
    
    df.loc[df["pred"]>=th, "class_numeric_pred"] = 1
    df.loc[df["pred"]<th, "class_numeric_pred"] = 0
    score_reverse = sklearn_metrics.precision_score(df["class_numeric"], df["class_numeric_pred"], zero_division=0)
    
    precision = max(score, score_reverse)
    print(f"\tPrecision score: {precision:.3f} at threshold: {th:.3f}")
    return precision

def get_recall_score(non_nan_result_df:pd.DataFrame, th):
    df = non_nan_result_df.copy(deep=True)
    df.loc[df["pred"]>=th, "class_numeric_pred"] = 0
    df.loc[df["pred"]<th, "class_numeric_pred"] = 1
    score = sklearn_metrics.recall_score(df["class_numeric"], df["class_numeric_pred"])
    
    df.loc[df["pred"]>=th, "class_numeric_pred"] = 1
    df.loc[df["pred"]<th, "class_numeric_pred"] = 0
    score_reverse = sklearn_metrics.recall_score(df["class_numeric"], df["class_numeric_pred"])
    
    recall = max(score, score_reverse)
    print(f"\tRecall score: {recall:.3f} at threshold: {th:.3f}")
    return recall

def get_accuracy_score(non_nan_result_df:pd.DataFrame, th):
    df = non_nan_result_df.copy(deep=True)
    df.loc[df["pred"]>=th, "class_numeric_pred"] = 0
    df.loc[df["pred"]<th, "class_numeric_pred"] = 1
    score = sklearn_metrics.accuracy_score(df["class_numeric"], df["class_numeric_pred"])
    
    df.loc[df["pred"]>=th, "class_numeric_pred"] = 1
    df.loc[df["pred"]<th, "class_numeric_pred"] = 0
    score_reverse = sklearn_metrics.accuracy_score(df["class_numeric"], df["class_numeric_pred"])
    
    accuracy = max(score, score_reverse)
    print(f"\tAccuracy score: {accuracy:.3f} at threshold: {th:.3f}")
    return accuracy

def get_matthews_corrcoef(non_nan_result_df:pd.DataFrame, th):
    df = non_nan_result_df.copy(deep=True)
    df.loc[df["pred"]>=th, "class_numeric_pred"] = 0
    df.loc[df["pred"]<th, "class_numeric_pred"] = 1
    score = sklearn_metrics.matthews_corrcoef(df["class_numeric"], df["class_numeric_pred"])
    
    df.loc[df["pred"]>=th, "class_numeric_pred"] = 1
    df.loc[df["pred"]<th, "class_numeric_pred"] = 0
    score_reverse = sklearn_metrics.matthews_corrcoef(df["class_numeric"], df["class_numeric_pred"])
    
    mcc = max(score, score_reverse)
    print(f"\tMCC score: {mcc:.3f} at threshold: {th:.3f}")
    return mcc
    
def get_balanced_accuracy_score(non_nan_result_df:pd.DataFrame, th):
    df = non_nan_result_df.copy(deep=True)
    df.loc[df["pred"]>=th, "class_numeric_pred"] = 0
    df.loc[df["pred"]<th, "class_numeric_pred"] = 1
    score = sklearn_metrics.balanced_accuracy_score(df["class_numeric"], df["class_numeric_pred"])
    
    df.loc[df["pred"]>=th, "class_numeric_pred"] = 1
    df.loc[df["pred"]<th, "class_numeric_pred"] = 0
    score_reverse = sklearn_metrics.balanced_accuracy_score(df["class_numeric"], df["class_numeric_pred"])
    
    balanced_accuracy = max(score, score_reverse)
    print(f"\tBalanced accuracy score: {balanced_accuracy:.3f} at threshold: {th:.3f}")
    return balanced_accuracy