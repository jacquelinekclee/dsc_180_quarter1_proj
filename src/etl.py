import pandas as pd
import json
import subprocess
from bed_reader import open_bed, sample_file
import numpy as np

def get_all_data(gene_expression_fp, genotype_fps, plink_output_fps, population_mapping_fp, populations = ['CEU', 'FIN', 'GBR', 'TSI', 'YRI']):
    # gene expression data
    population_mapping_df = get_population_sample_mapping(population_mapping_fp)
    # keep pops in the 5 populations
    population_mapping_df = population_mapping_df.loc[population_mapping_df.population.isin(populations)]
     
    gene_expr = process_gene_expression_data(gene_expression_fp)
    # gene expr sample ids
    gene_expr_samps = list(gene_expr.columns[6:])
    
    # keep sample ids in gene expr
    population_mapping_df = population_mapping_df.loc[population_mapping_df['sample_id'].isin(gene_expr_samps)]
    
    five_pop_samps = population_mapping_df['sample_id'].to_list()
    gene_expr = gene_expr[list(gene_expr.columns[:6]) + five_pop_samps].copy()
    
    # dict pop:samp ids
    pop_sample_ids = population_mapping_df.groupby('population')['sample_id'].agg('unique').to_dict()
    
    # finish gene expr
    gene_expr_transposed = get_gene_expr_transposed(gene_expr)
    gene_expr_info = get_gene_expression_info(gene_expr_transposed)
    
    # genotype data
    genotype_df = create_large_genotype_df(genotype_fps, plink_output_fps, five_pop_samps)
    
    # nearby snps 
    final_gene_expr_df = get_final_gene_expr_df(genotype_df, gene_expr, gene_expr_info)
    
    samps_genotype = list(genotype_df.columns[5:])
    
    pop_data = {'all':[final_gene_expr_df, genotype_df]}
    
    
    for pop in populations:
        pop_samps = list(pop_sample_ids[pop])
        gene_expr_pop = gene_expr[list(gene_expr.columns[:6]) + pop_samps]
        gene_expr_transposed_pop = get_gene_expr_transposed(gene_expr_pop)
        gene_expr_info_pop = get_gene_expression_info(gene_expr_transposed_pop)
        gene_expr_info_pop.drop_duplicates(inplace = True)
        
        pop_samps_genotype = list(set(samps_genotype).intersection(set(pop_samps)))
        genotype_df_pop = genotype_df[list(genotype_df.columns[:5]) + pop_samps_genotype]
        
        final_gene_expr_df_pop = get_final_gene_expr_df(genotype_df_pop, gene_expr_pop, gene_expr_info_pop)
        
        pop_data[pop] = [final_gene_expr_df_pop, genotype_df_pop]

    
    return pop_data
    
def get_population_sample_mapping(population_mapping_fp):
    population_mapping_df = pd.read_csv(population_mapping_fp, sep = ' ')
    population_mapping_df.rename(columns = {'sample':'sample_id'}, inplace=True)
    return population_mapping_df
    
def process_gene_expression_data(gene_expression_fp):
    gene_expr = pd.read_csv(gene_expression_fp, sep="\t", low_memory = False)
    if 0 in gene_expr.columns or '0' in gene_expr.columns:
        gene_expr.columns = gene_expr.loc[0]
        gene_expr = gene_expr.drop([0])
    gene_expr.Coord = gene_expr.Coord.astype(int)
    gene_expr['Coord_min'] = gene_expr.Coord - 500000
    gene_expr['Coord_max'] = gene_expr.Coord + 500000
    reorder_cols = list(gene_expr.columns[:4]) + list(gene_expr.columns[-2:]) + list(gene_expr.columns[4:-2])
    gene_expr = gene_expr[reorder_cols]
    return gene_expr
    
def get_gene_expr_transposed(gene_expr):
    gene_expr_transposed = gene_expr.melt(list(gene_expr.columns[:6]), var_name='sample_id', value_name='Value')
    return gene_expr_transposed


def get_gene_expression_info(gene_expr_transposed):
    gene_expr_info = gene_expr_transposed.drop(columns = ['sample_id', 'Value','TargetID'])
    gene_expr_info.drop_duplicates(inplace = True)
    return gene_expr_info

    
def process_one_genotype_data(genotype_fp, plink_output_fp):
    plink_process = "plink --vcf {} --biallelic-only --maf 0.05 --make-bed --out {}"
    result = subprocess.run(plink_process.format(genotype_fp, plink_output_fp),
                            check=True,
                            capture_output=True,
                            shell=True)
    output_files = {ext:plink_output_fp + ext for ext in ['.bim', '.bed', '.fam']}
    return output_files

def create_large_genotype_df(genotype_fps, plink_output_fps, gene_expr_samps):
    all_output_files = [process_one_genotype_data(genotype_fps[i], plink_output_fps[i])
                       for i in range(len(plink_output_fps))]
    dfs = [create_one_genotype_df(files, gene_expr_samps) for files in all_output_files]
    return pd.concat(dfs)
    
def create_one_genotype_df(output_files, gene_expr_samps):
    bed = open_bed(output_files['.bed'])
    vals = bed.read(index=np.s_[::,::])
    samples = bed.iid
    df = pd.read_csv(output_files['.bim'], sep = '\t', header = None)
    df.columns = ['Chr', 'ref_snp', 'gene_expression', 'pos', 'ref_allele', 'alt_allele']
    df.drop(columns = ['gene_expression'], inplace=True) 
    samples_df = pd.DataFrame(data = {samples[i]:vals[i] for i in range(len(samples))})
    genotype_df = pd.concat([df, samples_df], axis=1)
    cols = ['Chr', 'ref_snp', 'pos', 'ref_allele', 'alt_allele']
    # keep only biallelic
    genotype_df = genotype_df.loc[(genotype_df.ref_allele.str.len() == 1) & (genotype_df.alt_allele.str.len() == 1)].copy()
    return genotype_df.loc[genotype_df.ref_snp != '.'][cols + list(gene_expr_samps)].set_index('ref_snp')

def nearby_snps(row, genotype_df):
    matches = genotype_df.loc[genotype_df.pos.between(row.Coord_min, row.Coord_max)].index
    if len(matches) == 0:
        return []
    else:
        return matches.values
    
def get_final_gene_expr_df(genotype_df, gene_expr, gene_expr_info):
    genotype_df['Chr'] = genotype_df['Chr'].astype(str)
    gene_expr_info['snps'] = gene_expr_info.apply(lambda row: nearby_snps(row, genotype_df.loc[genotype_df.Chr == row.Chr]), axis = 1)
    final_gene_expr_df = gene_expr.merge(gene_expr_info, on = list(gene_expr_info.columns[:5]), how = 'left')
    return final_gene_expr_df

    
