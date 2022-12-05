import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import random

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

def generate_results(final_results_fps, data_dict, total_num_tests_all, total_num_tests_pops):
    genotype_df = data_dict['all'][-1]
    genotype_df = genotype_df.reset_index()
    gene_expr = data_dict['all'][0]
    all_results_fps = generate_all_results_viz(final_results_fps['all'], genotype_df, gene_expr, total_num_tests_all)
    population_fps = {key:val for key, val in final_results_fps.items() if key not in ['all', 'total_num_tests_pops', 'total_num_tests_all']}
        
    population_results_fps = generate_population_results_viz(population_fps, total_num_tests_pops)
    # return all_results_fps + population_results_fps
    return population_results_fps
    
def generate_all_results_viz(fp, genotype_df, gene_expr, total_num_tests_all):
    full_results = pd.read_csv(fp)
    full_results.rename(columns = {'Unnamed: 0':'ref_snp'}, inplace=True)
    full_results.loc[full_results.coeffs == 0, 'pvalue'] = 1
    full_results_merged = merged_results(full_results, genotype_df)
    full_results_merged_filt = full_results_merged.loc[full_results_merged.ref_snp != 'no matches']
    full_results_merged_filt_sig = full_results_merged_filt.loc[full_results_merged_filt.pvalue < (0.05 / total_num_tests_all)]
    # viz 1: bar chart 
    bar_chart = all_results_bar_chart(full_results_merged_filt_sig)
    # viz 2: box plot
    box_plot = all_results_most_sig_boxplot(full_results_merged_filt, genotype_df, gene_expr)
    return {'bar chart file path':bar_chart, 'box plot file path':box_plot}
    
def generate_population_results_viz(fps, total_num_tests_pops):
    final_results_fps = {}
    pop_results = {pop:pd.read_csv(fp) for pop, fp in fps.items()}
    for pop in pop_results:
        df = pop_results[pop]
        df.loc[df.coeffs == 0, 'pvalue'] = 1
        pop_results[pop] = df.copy()
    for pop in pop_results:
        df = pop_results[pop]
        ps = []
        for p in df.pos:
            try:
                x = int(p)
                ps.append(x)
            except:
                ps.append(int(p.split()[2]))
        df['pos'] = ps
        pop_results[pop] = df.copy()
    ren = [df.rename(columns = {'Unnamed: 0': 'ref_snp'}, inplace = True) for df in pop_results.values()]
    sig_results = {pop:val.loc[val.is_sig] for pop,val in pop_results.items()}
    try:
        counts = generate_population_counts(sig_results, pop_results)
    except:
        pass
    try:
        # viz 3: bar plot for number of cis-eQTLs for each population
        bar_chart = population_results_bar_chart(sig_results)
        final_results_fps['population bar chart file path'] = bar_chart
    except:
        pass
    try:
        # viz 4: stacked bar plot for number of genes that found x number of sig cis-eQTLs
        stacked_bar_chart = population_results_stacked_bar_chart(counts)
        final_results_fps['population stacked bar chart file path'] = stacked_bar_chart
                         
    except:
        pass
    try:
        merged = generate_population_merged(sig_results)
        all_overlap = population_all_overlap(merged)
        eu_overlap = population_eu_overlap(sig_results)
        results_df_dict = {'Metric':['Shared Across All Populations', 'Shared Across European Populations'],
                          'Number of Significant SNP-Gene Pairs':[all_overlap, eu_overlap]}
        # table for overlap summary
        results_df = pd.DataFrame(data = results_df_dict)
        results_df_fp = 'results/overlap_results_summary.csv'
        results_df.to_csv(results_df_fp, index = False)
        # table for overlap counts
        overlap_matrix = population_overlap_matrix(sig_results)
        overlap_matrix_fp = 'results/population_overlap_matrix.csv'
        overlap_matrix.to_csv(overlap_matrix_fp, index_label = 'population')
        final_results_fps['overlap results summary table file path'] = results_df_fp
        final_results_fps['population overlap matrix file path'] = overlap_matrix_fp
                         
    except:
        pass
    try:
        # viz 5: dot plot
        samp_for_dot_plot = merged.sample().iloc[0]
        dot_plot = dot_plot_ci(samp_for_dot_plot.ref_snp, samp_for_dot_plot.gene, merged)
        final_results_fps['dot plot for e-QTL across all populations file path'] = dot_plot,
                         
    except:
        pass
    try:
        # viz 5a: dot plot for unique YRI
        yri_pairs_filt = generate_yri_pairs(sig_results)
        samp = random.sample(yri_pairs_filt, 1)[0].split('; ')
        ref = samp[0]
        g = samp[1]
        yri_dot_plot = dot_plot_ci_yri(ref, g, pop_results)
        final_results_fps['dot plot for e-QTL exclusive to YRI population file path'] = yri_dot_plot
                         
    except:
        pass
    try:
        # viz 6: scatter plot
        highest_overlap_pop1 = overlap_matrix.max().idxmax()
        highest_overlap_pop2 = overlap_matrix.loc[highest_overlap_pop1].idxmax()
        pop_scatter_plot = two_populations_scatter_plot(merged, highest_overlap_pop1, highest_overlap_pop2)
        final_results_fps['scatter plot comparing 2 populations effect sizes'] = pop_scatter_plot
    except:
        pass
    return final_results_fps
    
def merged_results(results, genotype_df):
    genotype_df_upd_filt = genotype_df.loc[(genotype_df.ref_allele.str.len() == 1) & (genotype_df.alt_allele.str.len() == 1)][['ref_snp', 'pos', 'ref_allele', 'alt_allele']]
    genotype_df_upd_filt['pos'] = genotype_df_upd_filt.pos.astype(int)
    results_merged = results.merge(genotype_df_upd_filt, on = ['ref_snp', 'pos'], how = 'left')
    return results_merged
    
def all_results_bar_chart(full_results_merged_filt_sig):
    fig, ax = plt.subplots()
    x = full_results_merged_filt_sig.groupby('gene')['ref_snp'].count().value_counts().sort_index().index[:25]
    y = full_results_merged_filt_sig.groupby('gene')['ref_snp'].count().value_counts().sort_index().values[:25]
    plt.bar(x, y, edgecolor = 'black')
    plt.title('Number of Genes With X Significant SNPs')
    plt.xlabel('# significant SNPs')
    plt.ylabel('# Genes')
    date = datetime.today().strftime('%Y_%m_%d')
    fp = 'results/{}_genes_snp_count_bar_chart.png'.format(date)
    plt.savefig(fp,bbox_inches = "tight")
    return fp
    
def all_results_most_sig_boxplot(full_results_merged_filt, genotype_df, gene_expr):
    fig, ax = plt.subplots()
    most_sig = full_results_merged_filt.sort_values("pvalue").iloc[0]
    ref_snp_vals = genotype_df.loc[genotype_df.ref_snp == most_sig.ref_snp].melt(list(genotype_df.columns[:5]), var_name='sample_id', value_name='Value')
    gene_expr_vals = gene_expr.loc[gene_expr.Gene_Symbol == most_sig.gene].melt(list(gene_expr.columns[:7]), var_name='sample_id', value_name='gene_expression')
  
    ref_snp_vals_filt = ref_snp_vals[['sample_id', 'Value']]
    gene_expr_vals_filt = gene_expr_vals[['sample_id', 'gene_expression']]
    one_gene_snp_pair = ref_snp_vals_filt.merge(gene_expr_vals_filt, on = 'sample_id', how = 'inner')
    one_gene_snp_pair_filt = one_gene_snp_pair[['gene_expression', 'Value']]
    genotype_dict = {}
    genotype_dict[0] = most_sig.ref_allele * 2
    genotype_dict[1] = most_sig.ref_allele + most_sig.alt_allele
    genotype_dict[2] = most_sig.alt_allele * 2
    one_gene_snp_pair_filt['Value'] = one_gene_snp_pair_filt['Value'].replace(genotype_dict)
    one_gene_snp_pair_filt.rename(columns = {'Value':'genotype'}, inplace = True)
    one_gene_snp_pair_filt['gene_expression'] = one_gene_snp_pair_filt['gene_expression'].astype(float)
    one_gene_snp_pair_filt.boxplot('gene_expression', by = 'genotype')
    plt.suptitle('')
    plt.title('Gene: {}. SNP: {}.'.format(most_sig.gene, most_sig.ref_snp))
    plt.xlabel('genotype')
    date = datetime.today().strftime('%Y_%m_%d')
    fp = 'results/{}_most_sig_pair_boxplot.png'.format(date)
    plt.savefig(fp,bbox_inches = "tight")
    return fp
    
def population_results_bar_chart(sig_results):
    fig, ax = plt.subplots()
    pd.Series({pop:sig_results[pop].shape[0] for pop in sig_results}).plot(kind = 'bar', color = colors[:5], edgecolor = 'black')
    plt.xlabel('Population')
    plt.ylabel('# Significant SNPs')
    plt.title('# Significant SNPs for Each Population')
    date = datetime.today().strftime('%Y_%m_%d')
    fp = 'results/{}_num_sig_snps_by_population.png'.format(date)
    plt.savefig(fp,bbox_inches = "tight")
    return fp

def generate_population_counts(sig_results, pop_results):
    counts = {}
    indexes = {}
    populations = sig_results.keys()
    for pop in populations:
        c = sig_results[pop].groupby('gene')['pvalue'].count().value_counts()
        c[0] = pop_results[pop].shape[0] - c.sum()
        for i in c.index:
            if i not in indexes:
                indexes[i] = 1
            else:
                indexes[i] += 1
        counts[pop] = c.sort_index()

    for c in counts: 
        series = counts[c]
        for i in indexes:
            if i in series:
                continue
            else:
                series[i] = np.nan
        series.sort_index(inplace = True)
    return counts

def population_results_stacked_bar_chart(counts):
    fig, ax = plt.subplots()
    counts_dfs = []
    for pop in counts:
        counts_df = pd.DataFrame({'num_genes':counts[pop].iloc[1:16].values, 'num_sig_snps':counts[pop].iloc[1:16].index})
        counts_df['population'] = pop
        counts_dfs.append(counts_df)
    counts_dfs_conc = pd.concat(counts_dfs)
    counts_dfs_conc_plt = counts_dfs_conc.set_index('num_sig_snps').pivot(columns='population')
    p = counts_dfs_conc_plt.plot(kind = 'bar', stacked = True)
    plt.title('Number of Genes based on Number of cis-eQTLs')
    plt.ylabel('Number of Genes')
    plt.legend(['CEU', 'FIN', 'GBR', 'TSI', 'YRI'])
    date = datetime.today().strftime('%Y_%m_%d')
    fp = 'results/{}_genes_snp_count_by_pop_stacked_bar_chart.png'.format(date)
    plt.savefig(fp,bbox_inches = "tight")
    return fp
    
def generate_population_merged(sig_results):
    merged = sig_results['CEU']
    merged = merged.rename(columns = {col:col + '_CEU' for col in merged.columns if col not in ['ref_snp', 'gene']})
    prev = ''
    for pop in sig_results:
        if pop == 'CEU':
            continue
        curr = sig_results[pop]
        if prev != '':
            prev = '_' + prev
        merged = merged.merge(curr, on = ['ref_snp', 'gene'], how = 'inner', suffixes = [prev, '_' + pop])
        prev = pop
    return merged
    
def population_all_overlap(merged):
    return merged.shape[0]

def population_eu_overlap(sig_results):
    merged_eu = sig_results['CEU']
    merged_eu = merged_eu.rename(columns = {col:col + '_CEU' for col in merged_eu.columns if col not in ['ref_snp', 'gene']})
    prev = ''
    for pop in sig_results:
        if pop == 'CEU' or pop == 'YRI':
            continue
        curr = sig_results[pop]
        if prev != '':
            prev = '_' + prev
        merged_eu = merged_eu.merge(curr, on = ['ref_snp', 'gene'], how = 'inner', suffixes = [prev, '_' + pop])
        prev = pop

    merged_eu = merged_eu.rename(columns = {col:col + '_TSI' for col in merged_eu.columns if col not in ['ref_snp', 'gene'] and '_' not in col})
    return merged_eu.shape[0]

def population_overlap_matrix(sig_results):
    yri_gbr, yri_tsi, yri_ceu, yri_fin = yri_overlap(sig_results)
    gbr_tsi, gbr_ceu, gbr_fin = gbr_overlap(sig_results)
    tsi_ceu, tsi_fin = tsi_overlap(sig_results)
    merged_fin_ceu = sig_results['FIN'].merge(sig_results['CEU'], on = ['ref_snp', 'gene'], how = 'inner', suffixes = ['_fin', '_ceu'])
    ceu_fin = merged_fin_ceu.shape[0]
    yri_column = [np.nan, yri_gbr, yri_tsi, yri_ceu, yri_fin]
    gbr_column = [yri_gbr, np.nan, gbr_tsi, gbr_ceu, gbr_fin]
    tsi_column = [yri_tsi, yri_gbr, np.nan, tsi_ceu, tsi_fin]
    ceu_column = [yri_ceu, gbr_ceu, tsi_ceu, np.nan, ceu_fin]
    fin_column = [yri_fin, gbr_fin, tsi_fin, ceu_fin, np.nan]
    pops = ['YRI', 'GBR', 'TSI', 'CEU', 'FIN']
    matrix = pd.DataFrame(data = [yri_column, gbr_column, tsi_column, ceu_column, fin_column], columns = pops, index = pops)
    return matrix
    

def yri_overlap(sig_results):
    merged_yri_gbr = sig_results['YRI'].merge(sig_results['GBR'], on = ['ref_snp', 'gene'], how = 'inner', suffixes = ['_yri', '_gbr'])
    merged_yri_tsi = sig_results['YRI'].merge(sig_results['TSI'], on = ['ref_snp', 'gene'], how = 'inner', suffixes = ['_yri', '_tsi'])
    merged_yri_ceu = sig_results['YRI'].merge(sig_results['CEU'], on = ['ref_snp', 'gene'], how = 'inner', suffixes = ['_yri', '_ceu'])
    merged_yri_fin = sig_results['YRI'].merge(sig_results['FIN'], on = ['ref_snp', 'gene'], how = 'inner', suffixes = ['_yri', '_fin'])
    return merged_yri_gbr.shape[0], merged_yri_tsi.shape[0], merged_yri_ceu.shape[0], merged_yri_fin.shape[0]

def gbr_overlap(sig_results):

    merged_gbr_tsi = sig_results['GBR'].merge(sig_results['TSI'], on = ['ref_snp', 'gene'], how = 'inner', suffixes = ['_gbr', '_tsi'])
    merged_gbr_ceu = sig_results['GBR'].merge(sig_results['CEU'], on = ['ref_snp', 'gene'], how = 'inner', suffixes = ['_gbr', '_ceu'])
    merged_gbr_fin = sig_results['GBR'].merge(sig_results['FIN'], on = ['ref_snp', 'gene'], how = 'inner', suffixes = ['_gbr', '_fin'])
    return merged_gbr_tsi.shape[0], merged_gbr_ceu.shape[0], merged_gbr_fin.shape[0]

def tsi_overlap(sig_results):
    merged_tsi_ceu = sig_results['TSI'].merge(sig_results['CEU'], on = ['ref_snp', 'gene'], how = 'inner', suffixes = ['_tsi', '_ceu'])
    merged_tsi_fin = sig_results['TSI'].merge(sig_results['FIN'], on = ['ref_snp', 'gene'], how = 'inner', suffixes = ['_tsi', '_fin'])
    return merged_tsi_ceu.shape[0], merged_tsi_fin.shape[0]


def dot_plot_ci(ref, g, merged):

    row = merged.loc[(merged.gene == g) & (merged.ref_snp == ref)].iloc[0]
    
    fig, ax = plt.subplots()

    plt.xticks([1, 2, 3, 4, 5], ['GBR', 'FIN', 'CEU', 'TSI', 'YRI'])
    # 1
    x = row['coeffs_GBR']
    se = row['stderror_GBR']
    horizontal_line_width=0.25
    color='#2187bb'
    confidence_interval = 1.96 * se
    y = 1

    left = y - horizontal_line_width / 2
    top = x - confidence_interval
    right = y + horizontal_line_width / 2
    bottom = x + confidence_interval

    ax.plot([y, y], [top, bottom], color=color)
    ax.plot([left, right], [top, top], color=color)
    ax.plot([left, right], [bottom, bottom], color=color)

    ax.plot(y, x, 'o', color='#f44336')

    # 2
    x = row['coeffs_FIN']
    se = row['stderror_FIN']
    horizontal_line_width=0.25
    color='#2187bb'
    confidence_interval = 1.96 * se
    y = 2

    left = y - horizontal_line_width / 2
    top = x - confidence_interval
    right = y + horizontal_line_width / 2
    bottom = x + confidence_interval

    ax.plot([y, y], [top, bottom], color=color)
    ax.plot([left, right], [top, top], color=color)
    ax.plot([left, right], [bottom, bottom], color=color)

    ax.plot(y, x, 'o', color='#f44336')

    # 3: ceu
    x = row['coeffs_CEU']
    se = row['stderror_CEU']
    horizontal_line_width=0.25
    color='#2187bb'
    confidence_interval = 1.96 * se
    y = 3

    left = y - horizontal_line_width / 2
    top = x - confidence_interval
    right = y + horizontal_line_width / 2
    bottom = x + confidence_interval

    ax.plot([y, y], [top, bottom], color=color)
    ax.plot([left, right], [top, top], color=color)
    ax.plot([left, right], [bottom, bottom], color=color)

    plt.plot(y, x, 'o', color='#f44336')

    # 4
    x = row['coeffs_TSI']
    se = row['stderror_TSI']
    horizontal_line_width=0.25
    color='#2187bb'
    confidence_interval = 1.96 * se
    y = 4

    left = y - horizontal_line_width / 2
    top = x - confidence_interval
    right = y + horizontal_line_width / 2
    bottom = x + confidence_interval

    ax.plot([y, y], [top, bottom], color=color)
    ax.plot([left, right], [top, top], color=color)
    ax.plot([left, right], [bottom, bottom], color=color)

    plt.plot(y, x, 'o', color='#f44336')

    # 5
    x = row['coeffs_YRI']
    se = row['stderror_YRI']
    horizontal_line_width=0.25
    color='#2187bb'
    confidence_interval = 1.96 * se
    y = 5

    left = y - horizontal_line_width / 2
    top = x - confidence_interval
    right = y + horizontal_line_width / 2
    bottom = x + confidence_interval

    ax.plot([y, y], [top, bottom], color=color)
    ax.plot([left, right], [top, top], color=color)
    ax.plot([left, right], [bottom, bottom], color=color)

    ax.plot(y, x, 'o', color='#f44336')

    plt.title('Gene: {}; SNP: {}'.format(g, ref))
    date = datetime.today().strftime('%Y_%m_%d')
    fp = 'results/{}_dot_plot_with_confidence_intervals_by_pop.png'.format(date)
    plt.savefig(fp,bbox_inches = "tight")
    return fp

def generate_yri_pairs(sig_results):
    yri_non_zero_pval = sig_results['YRI']
    yri_pairs = yri_non_zero_pval.ref_snp + '; ' + yri_non_zero_pval.gene
    all_pairs = {}
    for pop, df in sig_results.items():
        pairs = df.ref_snp + '; ' + df.gene
        all_pairs[pop] = pairs


    all_pairs_filt = {k:v for k,v in all_pairs.items() if k != 'YRI'}

    def find_sim_1121(pair):
        for pop, pairs in all_pairs_filt.items():
            for p in pairs:
                if pair == p:
                    return False
        return True

    yri_pairs_filt = list(filter(find_sim_1121, yri_pairs))
    return yri_pairs_filt

def dot_plot_ci_yri(ref, g, pop_results):
    
    fig, ax = plt.subplots()

    plt.xticks([1, 2, 3, 4, 5], ['GBR', 'FIN', 'CEU', 'TSI', 'YRI'])
    # 1
    row = pop_results['GBR'].loc[(pop_results['GBR'].ref_snp == ref) & (pop_results['GBR'].gene == g)].iloc[0]
    x = row['coeffs']
    se = row['stderror']
    horizontal_line_width=0.25
    color='#2187bb'
    confidence_interval = 1.96 * se
    y = 1

    left = y - horizontal_line_width / 2
    top = x - confidence_interval
    right = y + horizontal_line_width / 2
    bottom = x + confidence_interval

    ax.plot([y, y], [top, bottom], color=color)
    ax.plot([left, right], [top, top], color=color)
    ax.plot([left, right], [bottom, bottom], color=color)

    ax.plot(y, x, 'o', color='#f44336')

    # 2
    row = pop_results['FIN'].loc[(pop_results['FIN'].ref_snp == ref) & (pop_results['FIN'].gene == g)].iloc[0]
    x = row['coeffs']
    se = row['stderror']
    horizontal_line_width=0.25
    color='#2187bb'
    confidence_interval = 1.96 * se
    y = 2

    left = y - horizontal_line_width / 2
    top = x - confidence_interval
    right = y + horizontal_line_width / 2
    bottom = x + confidence_interval

    ax.plot([y, y], [top, bottom], color=color)
    ax.plot([left, right], [top, top], color=color)
    ax.plot([left, right], [bottom, bottom], color=color)

    ax.plot(y, x, 'o', color='#f44336')

    # 3: ceu
    row = pop_results['CEU'].loc[(pop_results['CEU'].ref_snp == ref) & (pop_results['CEU'].gene == g)].iloc[0]
    x = row['coeffs']
    se = row['stderror']
    horizontal_line_width=0.25
    color='#2187bb'
    confidence_interval = 1.96 * se
    y = 3

    left = y - horizontal_line_width / 2
    top = x - confidence_interval
    right = y + horizontal_line_width / 2
    bottom = x + confidence_interval

    ax.plot([y, y], [top, bottom], color=color)
    ax.plot([left, right], [top, top], color=color)
    ax.plot([left, right], [bottom, bottom], color=color)

    ax.plot(y, x, 'o', color='#f44336')

    # 4
    row = pop_results['TSI'].loc[(pop_results['TSI'].ref_snp == ref) & (pop_results['TSI'].gene == g)].iloc[0]
    x = row['coeffs']
    se = row['stderror']
    horizontal_line_width=0.25
    color='#2187bb'
    confidence_interval = 1.96 * se
    y = 4

    left = y - horizontal_line_width / 2
    top = x - confidence_interval
    right = y + horizontal_line_width / 2
    bottom = x + confidence_interval

    ax.plot([y, y], [top, bottom], color=color)
    ax.plot([left, right], [top, top], color=color)
    ax.plot([left, right], [bottom, bottom], color=color)

    plt.plot(y, x, 'o', color='#f44336')

    # 5
    row = pop_results['YRI'].loc[(pop_results['YRI'].ref_snp == ref) & (pop_results['YRI'].gene == g)].iloc[0]
    x = row['coeffs']
    se = row['stderror']
    horizontal_line_width=0.25
    color='#2187bb'
    confidence_interval = 1.96 * se
    y = 5

    left = y - horizontal_line_width / 2
    top = x - confidence_interval
    right = y + horizontal_line_width / 2
    bottom = x + confidence_interval

    ax.plot([y, y], [top, bottom], color=color)
    ax.plot([left, right], [top, top], color=color)
    ax.plot([left, right], [bottom, bottom], color=color)

    ax.plot(y, x, 'o', color='#f44336')

    plt.title('Gene: {}; SNP: {}'.format(g, ref))
    date = datetime.today().strftime('%Y_%m_%d')
    fp = 'results/{}_dot_plot_with_confidence_intervals_by_pop_yri.png'.format(date)
    plt.savefig(fp,bbox_inches = "tight")
    return fp

def two_populations_scatter_plot(merged, pop1, pop2):
    x = merged['coeffs_{}'.format(pop1)]
    y = merged['coeffs_{}'.format(pop2)]
    fig, ax = plt.subplots()
    axis_range = range(int(min(x.min(), y.min())), int(max(x.max(), y.max())))
    ax.plot(axis_range, axis_range, color = 'black')
    plt.scatter(x, y)
    plt.title('{} vs. {}'.format(pop1, pop2))
    plt.xlabel('GBR')
    plt.ylabel('CEU')
    date = datetime.today().strftime('%Y_%m_%d')
    fp = 'results/{}_{}_vs_{}_scatter_plot.png'.format(date, pop1, pop2)
    plt.savefig(fp,bbox_inches = "tight")
    return fp

