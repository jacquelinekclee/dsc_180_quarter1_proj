import pandas as pd
import logging
import statsmodels.api as sm
from datetime import datetime
import logging
import numpy as np

def run_all_linear_regressions(data_dict):
    final_results_fps = {}
    
    total_num_tests_pops = 0
    total_num_tests_all = 0
    
    date = datetime.today().strftime('%Y_%m_%d')
    
    for key, val in data_dict.items():
        final_gene_expr_df = val[0]
        num_tests = final_gene_expr_df.snps.str.len().sum()
        if key != 'all':
            total_num_tests_pops += num_tests
        else:
            total_num_tests_all += num_tests
    final_results_fps['total_num_tests_all'] = total_num_tests_all
    final_results_fps['total_num_tests_pops'] = total_num_tests_pops
    # perform analysis for each population        
    
    for key, val in data_dict.items():
        logging.info('start analysis for {}'.format(key))
        if key == 'all':
            total_num_tests = total_num_tests_all
        else:
            total_num_tests = total_num_tests_pops
            
        logging.info('total number of tests for {}: {}'.format(key, str(total_num_tests)))
        final_gene_expr_df = val[0]
        genotype_df = val[1]
        final_gene_expr_df['results'] = final_gene_expr_df.apply(lambda row: linear_regressions_upd(row, genotype_df, total_num_tests), axis = 1)
        final_results = process_results(final_gene_expr_df)
        # key = population
        fp = 'data/out/{}_eqtl_analysis_results_chr22_{}.csv'.format(date, key)
        final_results.to_csv(fp)
        final_results_fps[key] = fp
        logging.info('finished analysis for {}'.format(key))
    print(final_results_fps.keys())
    return final_results_fps
    
def linear_regressions_upd(row, genotype_df, total_num_tests):
    snps = row['snps']
    gene = row.Gene_Symbol
    if len(snps) == 0:
        return {'no matches':[-1., -1., -1., -1., gene, -1., False]}
    # get genotype for given snps
    genotypes = genotype_df.loc[snps]
    genotypes_t = genotypes.iloc[:, 4:].transpose()
    gene_expr = row[7:].to_frame(name = 'y')
    # merge on sample_id
    merged = genotypes_t.merge(gene_expr, left_index=True, right_index=True)
    # gene expr data
    y = merged.y.astype(float)
    # Xs = genotype data
    Xs = merged.drop(columns = ['y'])
    results = {}
    p_val_thresh = 0.05 / total_num_tests
    for snp in Xs:
        X = Xs[snp].values
        if len(X) == 0:
            results[snp] = [-1., -1., -1., -1., gene, -1., False]
            continue
        pos = genotype_df.loc[snp].pos
        if type(pos) == pd.Series:
            pos = pos.iloc[0]
        # run linear regression
        X = sm.add_constant(X)
        model = sm.OLS(y, X)
        fit_model = model.fit()
        if 'x1' not in fit_model.params.index:
            results[snp] = [-1., -1., -1., -1., gene, -1., False]
            continue
        coeffs = fit_model.params['x1']
        if np.isnan(coeffs):
            results[snp] = [-1., -1., -1., -1., gene, -1., False]
            continue
        p_val = fit_model.pvalues.x1
        std_error = fit_model.bse['x1']
        sig_bool = p_val < p_val_thresh
        results[snp] = [p_val, fit_model.rsquared, std_error, pos, gene, coeffs, sig_bool]
    return results

def process_results(final_gene_expr_df):
    # extract results from dicts
    final_results = pd.concat([pd.DataFrame(final_gene_expr_df['results'].iloc[i]).transpose()
                               for i in range(len(final_gene_expr_df))])
    final_results.columns = ['pvalue', 'rsquared', 'stderror', 'pos', 'gene', 'coeffs', 'is_sig']
    final_results.fillna(0, inplace=True)
    final_results['minuslog10pvalue'] = -np.log10(final_results.pvalue.astype(float))
    final_results.sort_values('minuslog10pvalue', ascending=False, inplace=True)
    return final_results
    
