#!/usr/bin/env python

# import argparse
# parser = argparse.ArgumentParser()
# parser.add_argument("data", type=str)
# args = parser.parse_args()
# print(args.square**2)

import sys
import json

sys.path.insert(0, 'src')

from etl import get_all_data
from eqtl_analysis import run_all_linear_regressions
from visualize import generate_results
from datetime import datetime
import logging

import os
import os.path
import shutil

date = datetime.today().strftime('%Y_%m_%d')
logging.basicConfig(filename='log_{}.txt'.format(date), 
		    filemode='a', 
		    level=logging.INFO,
		    datefmt='%H:%M:%S',
		    format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s')

def main(targets):
    if 'test' in targets:
        with open('config/data-params.json') as fh:
            file_paths = json.load(fh)['test']
            logging.info('starting test process')
    elif 'data' in targets:
        with open('config/data-params.json') as fh:
            file_paths = json.load(fh)['data']
            logging.info('starting all process')
    logging.info('extract file paths')
    gene_expression_fp, genotype_fps, plink_output_fps, population_mapping_fp = file_paths.values()
    logging.info('process all data')
    data_dict = get_all_data(gene_expression_fp, genotype_fps, plink_output_fps, population_mapping_fp)
    logging.info('run linear regressions')
    final_results_fps = run_all_linear_regressions(data_dict)
    total_num_tests_all = final_results_fps['total_num_tests_all']
    total_num_tests_pops = final_results_fps['total_num_tests_pops']
    final_results_fps = {k:v for k, v in final_results_fps.items() if k not in ['total_num_tests_all', 'total_num_tests_pops']}
    viz_results_fps = generate_results(final_results_fps, data_dict, total_num_tests_all, total_num_tests_pops)
    return viz_results_fps
        

if __name__ == '__main__':
    targets = sys.argv[1:]
    final_results_fps = main(targets)
    print(len(final_results_fps) > 0)
    
