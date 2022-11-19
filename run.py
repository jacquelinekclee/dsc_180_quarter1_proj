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
from datetime import datetime
import logging

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
        logging.info('extract file paths')
        gene_expression_fp, genotype_fps, plink_output_fps, population_mapping_fp = file_paths.values()
        logging.info('process all data')
        data_dict = get_all_data(gene_expression_fp, genotype_fps, plink_output_fps, population_mapping_fp)
        logging.info('run linear regressions')
        final_results_fps = run_all_linear_regressions(data_dict)
        return final_results_fps
        

if __name__ == '__main__':
    targets = sys.argv[1:]
    final_results_fps = main(targets)