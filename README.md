# code-examples

This repositories contains 3 scripts with three different examples of code.
1) ukb_dosages_newrelease.R is uses a .gen file generated using QCTOOL from a UK Biobank .bgen file. The file contains genotype data from a subset of SNPs from the UK Biobank dataset. The data matrix has SNPs in the rows and subjects in the columns and there are 3 columns per subject. This script efficiently converts the three columns to one column of dosage data and transposes the data matrix in order to obtain a final matrix with subject IDs in the row and SNP RSID's in the columns. This is saved to a .RData file (there is also an option to save to a .csv file).

Data files used by this script:
- ukb1913_imp_chr1_v2_s487398.sample - template sample file with the subject IDs. These are in the same order as the subjects in the columns of the subset.gen file (which does not contain meaningful column names). The subject IDs will be used as the ID column in the final output matrix
- subset.gen


--------------------------------------------------------------------------------------------------------------------------------------
2) rank and var ratios - empirical.R - this script processes and prepares genotype data from the ECRHS dataset to be analysed and performs the standard statistical analysis (with and without adjusting for other covariates) and a bayesian joint model with and without the inclusion of prior bioligical knowledge. Prior biological knowledge corresponds a set of 10 binary questions (0=No, 1=Yes) about biological characteristics of the SNPs.
The Bayesian approach corresponds to a Bayesian joint shrinkage model similar to Bayesian ridge regression and biological information informs the shrinkage applied to each of the SNPs. The model itself was implemented in the BhGLM package developed by Yi and Ma 2011 (available to donwload in the BhGLM repository in my github account). Details on how prior knowledge was incorporated in the Bayesian model thought modulation of shrinkage and the concept of the variance ratio can be found in the following paper: Pereira et al. Genetic Epidemioloy 2017.
The output of the analyses is the rank of the LD block of the best SNPs ranked according to their frequentist or Bayesian p-values, respectively. At the beggining of the script there are some processing steps assigning the SNPs to groups according to their LD blocks. This is important and the rationale of the approach is to identify small genomic regions that contain the causal SNP.


Data files used by this script: 
- fulldata.csv - phenotype and genotype data. Contains height and weight (used to calculate bmi in the script) and other variables which were adjusted for in the analyses. Also contains the dosage data for the SNPs used in this analyses.
- priorAnswers.csv - matrix of prior knowledge. Contains the binary answers to each prior knowledge question. It is used to calculate a prior knowledge score. There are three alternatives in the script corresponding to a simple sum of the score in each question (object qSum in the code) and two others where each question has a different weight (weighted Sum). These are commented out in the current script.
- snp6_24_ldblock.csv - Informaiton about the SNP (e.g. Chr and position) and allocation to LD blocks


--------------------------------------------------------------------------------------------------------------------------------------
3) ggplot example (data vizualization) - analysis of prior knowledge.R processes data of SNP prior knowledge data from the top hits of 2 different analyses of the UK Biobank data in the context of a genetic association study on lung function (outcome measures: FVC and the ratio FEV1/FVC). For each outcome, there is the standard statistical analysis (called "standard" in the script), where each SNP was analysed separately, and a bayesian joint model where prior biological knowledge about the SNPs was included in the analysis (called "bayes" in the script). Prior biological knowledge consisted of a set of 10 binary questions (0=No, 1=Yes) about biological characteristics of the SNPs. The purpose of this script was to build several plots to:
  - compare the distribution of prior biological knowledge across the two types of analyses
  - study the distribution of prior biological knowledge per question
  - check which prior knowledge questions have more replicated hits
  - compare the different performance of the model between the two outcomes: FVC and ratio FEV1/FVC

Data files used by this script: 
- fvctop500000withLDblocks.tsv - top 50,000 hits from the standard analysis, outcome: FVC
- ratiotop500000withLDblocks.tsv - top 50,000 hits from the standard analysis, outcome: ratio FEV1/FVC
- priorknowledge_fvc.csv - matrix of prior knowledge for FVC
- priorknowledge_ratio.csv - matrix of prior knowledge for the ratio FEV1/FVC
- fvc standard results with eaf.csv - top hits for the standard analysis, outcome: FVC
- fvc bayes results with eaf.csv  - top hits for the bayesian analysis, outcome: FVC
- ratio standard results with eaf.csv - top hits for the standard analysis, outcome: ratio FEV1/FVC
- ratio bayes results with eaf.csv - top hits for the bayesian analysis, outcome: ratio FEV1/FVC
