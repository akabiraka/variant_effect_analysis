import sys
sys.path.append("../variant_effect_analysis")

import numpy as np
import pandas as pd
from sklearn import metrics as sklearn_metrics

def get_auc_roc_score(non_nan_result_df:pd.DataFrame):
    larger_means_positive_class = True
    df = non_nan_result_df.copy(deep=True)
    auc_roc_score = sklearn_metrics.roc_auc_score(df["class_numeric"], df["pred"])
    if auc_roc_score < 0.5: 
        larger_means_positive_class = False
    print(f"\tAUC-ROC: {auc_roc_score:.3f}")
    return auc_roc_score, larger_means_positive_class

def get_auc_pr_score(non_nan_result_df:pd.DataFrame):
    df = non_nan_result_df.copy(deep=True)
    precision, recall, thresholds = sklearn_metrics.precision_recall_curve(df["class_numeric"], df["pred"], pos_label=1)
    auc_pr_score = sklearn_metrics.auc(recall, precision)

    # if auc_pr_score < 0.5: 
    #     auc_pr_score = 1 - auc_pr_score
    print(f"\tAUC-PR: {auc_pr_score:.3f}")
    return auc_pr_score, precision, recall, thresholds

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
    return f1_max, th_max, precisions, recalls, thresholds

def get_f1max_and_th_(non_nan_result_df:pd.DataFrame):
    df = non_nan_result_df.copy(deep=True)
    f1_scores, precisions, recalls = [], [], []
    thresholds = np.arange(0, 1, .01)

    for th in thresholds:
        df.loc[df["pred"]>=th, "class_numeric_pred"] = 1
        df.loc[df["pred"]<th, "class_numeric_pred"] = 0
        prec = sklearn_metrics.precision_score(df["class_numeric"], df["class_numeric_pred"], pos_label=1, average="weighted")
        rec = sklearn_metrics.recall_score(df["class_numeric"], df["class_numeric_pred"], pos_label=1, average="weighted")
        f1 = sklearn_metrics.f1_score(df["class_numeric"], df["class_numeric_pred"], pos_label=1, average="weighted")

        precisions.append(prec)
        recalls.append(rec)
        f1_scores.append(f1)

    f1_max = np.max(f1_scores)
    th_max = thresholds[np.argmax(f1_scores)]
    print(f"\tcBest F1-Score: {f1_max:.3f} at threshold: {th_max:.3f}")
    return f1_max, th_max, precisions, recalls, thresholds

    
# def get_f1max_and_th_(precision, recall, thresholds):
#     numerator = 2 * recall * precision
#     denom = recall + precision
#     f1_scores = np.divide(numerator, denom, out=np.zeros_like(denom), where=(denom!=0))
#     f1_max = np.max(f1_scores)
#     th_max = thresholds[np.argmax(f1_scores)]
#     print(f"\tBest F1-Score: {f1_max:.3f} at threshold: {th_max:.3f}")
#     return f1_max, th_max

def get_precision_score(non_nan_result_df:pd.DataFrame, th):
    df = non_nan_result_df.copy(deep=True)
    df.loc[df["pred"]>=th, "class_numeric_pred"] = 1
    df.loc[df["pred"]<th, "class_numeric_pred"] = 0
    score = sklearn_metrics.precision_score(df["class_numeric"], df["class_numeric_pred"], zero_division=0)
    
    # df.loc[df["pred"]>=th, "class_numeric_pred"] = 1
    # df.loc[df["pred"]<th, "class_numeric_pred"] = 0
    # score_reverse = sklearn_metrics.precision_score(df["class_numeric"], df["class_numeric_pred"], zero_division=0)
    
    # precision = max(score, score_reverse)
    precision = score
    print(f"\tPrecision score: {precision:.3f} at threshold: {th:.3f}")
    return precision

def get_recall_score(non_nan_result_df:pd.DataFrame, th):
    df = non_nan_result_df.copy(deep=True)
    df.loc[df["pred"]>=th, "class_numeric_pred"] = 1
    df.loc[df["pred"]<th, "class_numeric_pred"] = 0
    score = sklearn_metrics.recall_score(df["class_numeric"], df["class_numeric_pred"])
    
    # df.loc[df["pred"]>=th, "class_numeric_pred"] = 1
    # df.loc[df["pred"]<th, "class_numeric_pred"] = 0
    # score_reverse = sklearn_metrics.recall_score(df["class_numeric"], df["class_numeric_pred"])
    
    # recall = max(score, score_reverse)
    recall = score
    print(f"\tRecall score: {recall:.3f} at threshold: {th:.3f}")
    return recall

def get_accuracy_score(non_nan_result_df:pd.DataFrame, th):
    df = non_nan_result_df.copy(deep=True)
    df.loc[df["pred"]>=th, "class_numeric_pred"] = 1
    df.loc[df["pred"]<th, "class_numeric_pred"] = 0
    score = sklearn_metrics.accuracy_score(df["class_numeric"], df["class_numeric_pred"])
    
    # df.loc[df["pred"]>=th, "class_numeric_pred"] = 1
    # df.loc[df["pred"]<th, "class_numeric_pred"] = 0
    # score_reverse = sklearn_metrics.accuracy_score(df["class_numeric"], df["class_numeric_pred"])
    
    # accuracy = max(score, score_reverse)
    accuracy = score
    print(f"\tAccuracy score: {accuracy:.3f} at threshold: {th:.3f}")
    return accuracy

def get_matthews_corrcoef(non_nan_result_df:pd.DataFrame, th):
    df = non_nan_result_df.copy(deep=True)
    df.loc[df["pred"]>=th, "class_numeric_pred"] = 1
    df.loc[df["pred"]<th, "class_numeric_pred"] = 0
    score = sklearn_metrics.matthews_corrcoef(df["class_numeric"], df["class_numeric_pred"])
    
    # df.loc[df["pred"]>=th, "class_numeric_pred"] = 1
    # df.loc[df["pred"]<th, "class_numeric_pred"] = 0
    # score_reverse = sklearn_metrics.matthews_corrcoef(df["class_numeric"], df["class_numeric_pred"])
    
    # mcc = max(score, score_reverse)
    mcc = score
    print(f"\tMCC score: {mcc:.3f} at threshold: {th:.3f}")
    return mcc
    
def get_balanced_accuracy_score(non_nan_result_df:pd.DataFrame, th):
    df = non_nan_result_df.copy(deep=True)
    df.loc[df["pred"]>=th, "class_numeric_pred"] = 1
    df.loc[df["pred"]<th, "class_numeric_pred"] = 0
    score = sklearn_metrics.balanced_accuracy_score(df["class_numeric"], df["class_numeric_pred"])
    
    # df.loc[df["pred"]>=th, "class_numeric_pred"] = 1
    # df.loc[df["pred"]<th, "class_numeric_pred"] = 0
    # score_reverse = sklearn_metrics.balanced_accuracy_score(df["class_numeric"], df["class_numeric_pred"])
    
    # balanced_accuracy = max(score, score_reverse)
    balanced_accuracy = score
    print(f"\tBalanced accuracy score: {balanced_accuracy:.3f} at threshold: {th:.3f}")
    return balanced_accuracy

from scipy.stats import ks_2samp
def get_KS_test_score(df:pd.DataFrame):
    pos_cls_preds = df[df["class_numeric"]==1]["pred"]
    neg_cls_preds = df[df["class_numeric"]==0]["pred"]
    
    res = ks_2samp(pos_cls_preds, neg_cls_preds)
    print(f"\tKS-test score. statistic: {res.statistic:.3f}, p-value: {res.pvalue:.3f}")
    return res.statistic, res.pvalue