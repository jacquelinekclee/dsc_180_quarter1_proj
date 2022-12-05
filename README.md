# DSC 180A Quarter 1 Project: Identifying cis-eQTLs in European and African Populations Using Linear Regression

## How to Access Raw Data

- Gene expression data: Navigate to [this link](https://zenodo.org/record/6998955#.Y4217C-B00c) and download the file with path "GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz". Unzip the file and save it as "GD462.GeneQuantRPKM.50FN.samplename.resk10.txt" in the data/raw directory.
- Genotype data: Download the file [linked here](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz). Save this file as "ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz" in the data/raw directory.
- Sample ID Mapping: Navigate to [this link](https://drive.google.com/file/d/17YcKfMtpKItXDuC5kJjUhth6aA7XycXA/view?usp=sharing) and download the file. Save it as "ALL_1000G_phase1integrated_v3.sample" in the data/raw directory. 

After completing the above, you should have the following files:

```
├── data             
    ├── raw           
    │   └── GD462.GeneQuantRPKM.50FN.samplename.resk10.txt
    │   └── ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz
    │   └── ALL_1000G_phase1integrated_v3.sample
  
```
