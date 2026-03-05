#################################################################################
# Version: 1.0
# Date: 08.02.2021
# Tong Wu, tong.wu@helmholtz-muenchen.de
# Python version: 3.8
# python EV_estimationLightGBM_genetics.py -path_to_X Genetics_data_FINAL_dec2023.csv -path_to_Y metabolomics_data_for_log2_FINAL_dec2023.csv -model lightgbam 'Pathtooutput/' -path_to_link linkfile2.csv
##################################################################################

import numpy as np
import pandas as pd
import lightgbm as lgb
import os
import pickle
import argparse
from sklearn.svm import SVR
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import Lasso
from sklearn.metrics import r2_score, precision_recall_curve, explained_variance_score
from sklearn.model_selection import GroupKFold, RandomizedSearchCV
from scipy.stats.stats import spearmanr, pearsonr


# constant hyper parameters for all runs and permutations - Suggested by Eran 5.12.2018
# The lgb model the params referenced by the A reference map of potential determinants for the human serum metabolome
lgb = lgb.LGBMRegressor(learning_rate=0.01, max_depth =5, feature_fraction=0.8, num_leaves = 25, min_data_in_leaf=15,
                        metric='l2', early_stopping_rounds=None, n_estimators=200, bagging_fraction=0.9,
                        bagging_freq=5, num_threads=1, verbose=-1, silent=True)

lasso = Lasso()


def cross_validation(command_args, Y, idx):
    X_genetics = pd.read_csv(command_args.path_to_X, index_col=0)
    link = pd.read_csv(command_args.path_to_link,index_col=0)
    df_link = pd.DataFrame(link)
    results_df = pd.DataFrame(index=Y.columns)
    predictions_df = pd.DataFrame(index=Y.index, columns=Y.columns)

    for y_name in Y.columns:
        subset=df_link.loc[df_link['phenotype'] == y_name]
        snps = pd.DataFrame(subset, columns=['ID'])
        snp_list=np.array(snps.index.array)
        snplist=snp_list.tolist()
        X=X_genetics.loc[:, snplist]
        y = Y[y_name].dropna().astype(float).copy()
        X_temp = X.loc[y.index].dropna(how="all").copy()
        y = y.loc[X_temp.index]

        groups = np.array(range(X.shape[0]))
        group_kfold = GroupKFold(n_splits=5)
        final_pred = pd.DataFrame(index=X_temp.index, columns=[y_name])

        try:
            for train_index, test_index in group_kfold.split(X, y, groups):
                X_train, X_test = X.iloc[train_index,
                                         :].copy(), X.iloc[test_index, :].copy()
                y_train, y_test = y.iloc[train_index].copy(
                ), y.iloc[test_index].copy()
                model = _choose_model(command_args.model)
                model.fit(X_train, y_train)
                y_pred = model.predict(X_test)
                final_pred.loc[X_test.index, :] = np.expand_dims(y_pred, 1)
            results_df = _evaluate_performance(
                y_name, final_pred.values.ravel(), y, results_df)
            predictions_df.loc[final_pred.index,
                               y_name] = final_pred.values.ravel()
        except:
            continue
    _save_temporary_files(command_args, idx, results_df, predictions_df)
    return

#chose which model would be used 
def _choose_model(model):
    if model == 'lightgbam':
        return lgb
    elif modle == 'lasso':
        return lasso


#svae the temporay result predicted by the models
def _save_temporary_files(command_args, idx, results_df, predictions_df):
    with open(command_args.output_dir + '/temp_resdf_' + str(idx) + '.pkl', 'wb') as fout:
        pickle.dump(results_df, fout)
    with open(command_args.output_dir + '/temp_pred_' + str(idx) + '.pkl', 'wb') as fout:
        pickle.dump(predictions_df, fout)
    return


def _evaluate_performance(y_name, y_pred, y_test, results_df):
    results_df.loc[y_name, 'Size'] = y_pred.shape[0]
    results_df.loc[y_name, 'Coefficient_of_determination'] = r2_score(
        y_true=y_test, y_pred=y_pred)
    results_df.loc[y_name, 'explained_variance_score'] = explained_variance_score(
        y_true=y_test, y_pred=y_pred)
    results_df.loc[y_name, 'pearson_r'], results_df.loc[y_name,
                                                        'pearson_p'] = pearsonr(y_pred, y_test)
    results_df.loc[y_name, 'spearman_r'], results_df.loc[y_name,
                                                         'spearman_p'] = spearmanr(y_pred, y_test)
    return results_df


def concat_outputs(command_args):
    all_temp_files = os.listdir(command_args.output_dir)
    resdf_files = [command_args.output_dir +
                   f for f in all_temp_files if f.startswith('temp_resdf_')]
    pred_files = [command_args.output_dir +
                  f for f in all_temp_files if f.startswith('temp_pred_')]
    _concat_files(resdf_files, command_args.output_dir +
                  '/results.pkl', how='dataframe')
    _concat_files(pred_files, command_args.output_dir +
                  '/predictions_df.pkl', how='dataframe', axis=1)
    return


def _concat_files(files, final_path, how='dataframe', axis=0):
    if how == 'dataframe':
        final_file = pd.DataFrame()
        for f in files:
            final_file = pd.concat((final_file, pd.read_pickle(f)), axis=axis)
            os.remove(f)
        with open(final_path, 'wb') as fout:
            pickle.dump(final_file, fout)
        final_file.to_csv(
            ((final_path.split('.pkl')[0]).split('.dat')[0]) + '.csv')
    elif how == 'dic':
        final_file = {}
        for f in files:
            final_file.update(pd.read_pickle(f))
            os.remove(f)
        with open(final_path, 'wb') as fout:
            pickle.dump(final_file, fout)
    return


def upload_these_jobs(command_args):
    if command_args.path_to_Y.endswith('.csv'):
        Y = pd.read_csv(command_args.path_to_Y, index_col=0)
    else:
        Y = pd.read_pickle(command_args.path_to_Y)

    for idx in range(0, Y.shape[1], command_args.n_cols_per_job):
        cross_validation(
            command_args, Y.iloc[:, idx:idx + command_args.n_cols_per_job], idx)

    # merge the temp results files
    concat_outputs(command_args)
    return


def make_dir_if_not_exists(path):
    if not os.path.exists(path):
        os.makedirs(path)


def main():
    print('main')
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'output_dir', help='Path to output directory', type=str, default=None)
    parser.add_argument(
        '-model', help='Which prediction model to use', type=str, default='lightgbm')
    parser.add_argument('-n_cols_per_job',
                        help='Number of columns per job', type=int, default=10)
    parser.add_argument('-path_to_X', '--path_to_X',
                        help='Path to features data - X', type=str, default='')
    parser.add_argument(
        '-path_to_Y', help='Path to labels - Y', type=str, default='')
    parser.add_argument(
        '-only_concat', help='Whether to only concatenate the output files', type=bool, default=False)
    parser.add_argument(
        '-path_to_link', help='Path to links of X to Y', type=str, default='')
    command_args = parser.parse_args()

    if (not os.path.exists(command_args.path_to_X)) or (not os.path.exists(command_args.path_to_Y)):
        print("X or Y doesn't exist!")
        return
    if command_args.n_cols_per_job < 1 or command_args.n_cols_per_job > 1000:
        print("n_cols_per_job must be between 1 and 1000")
        return

    if command_args.only_concat:
        concat_outputs(command_args)
        return

    make_dir_if_not_exists(command_args.output_dir)
    upload_these_jobs(command_args)


if __name__ == "__main__":
    main()

