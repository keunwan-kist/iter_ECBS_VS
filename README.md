# Iteative evolutionary chemical binding similarity searches with experimental validation data 

## Description
Scripts to update, train, test, and screen ECBS (Target-ECBS)


## Directory: 
scripts/ - scripts files to train and test ECBS 

example/ - test examples for ECBS update, train, and test (example outputs in Q02750.PP_NP.tgz)

data_files/ - data files used by scripts


## Prerequisites: 
Before running ECBS script, please install the prerequisite programs and edit the path variables in the scripts. 

R - tested version 3.4.2

Perl - tested version v5.16.3 

R ranger package (0.8.0 or higher)

ChemmineOB package, 
```
source("https://bioconductor.org/biocLite.R")
biocLite("ChemmineOB") 
biocLite("ChemmineR")
```

## Test procedure for ECBS model with MEK1 (Q02750) validation data
### Unzip chemical feature files 
```
cd data_files/
bunzip2 mat.storable.data.bz2
bunzip2 BindingDB_All_2D_v2.1uM_Affinity.bz2
```

### Initial screening data for MEK1
* MEK1 (Q02750) validation data from [Durai et al. BMC Bioinformatics (2020)](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-03643-x)  
* chemical data file: data_files/new_mol.txt   
|ChemID|POC|IC50|Class|Target|
|------|---|----|-----|------|
|27846780|100|NA|0|Q02750|

```
cd example

# chem pair data generation with a PP_NP scheme 
perl ../scripts/update_PN_chempair_to_ECBS.pl ../data_files/int.drug_target_map.txt ../data_files/new_mol.txt ../data_files/BindingDB_All_2D_v2.1uM_Affinity Q02750 NA Q02750.PP_NP PP_NP >& log.Q02750.PP_NP.txt

# train/test set generation for cross validation (7-fold CV)
perl ../scripts/split_testdata_for_ensemble_model_learning_TN.pl Q02750.PP_NP/ ../data_files/new_mol.txt 7

# trained model test 
cd example/Q02750.PP_NP/
Rscript ../../scripts/splitset_train_and_test_target-wise.r trainset1.target.data testset1.target.data > pr_auc.trainset1.target.data.out
```

### Output files 
* AUC in PR curve
**pr_auc.trainset1.target.data.out**
* test score for OOB samples (Random Forest)
**oob.trainset1.target.testset1.target.score**
* test score for indepenent test set1 
**result.trainset1.target.testset1.target.score**

## Virtual Screening with Target-ECBS 
```
cd example

# Random Forest model built with chemical pair data (ERCPs by evolutionary target relationships)
Rscript ../scripts/Indv_Model_Generator.r Q02750.PP_NP/target.data target.data.rf_model

# Chemical pair data generation with known active chems for MEK1
perl ../scripts/gen_pair_FP.pl -fp1 reduced.test.data -fp2 Q02750.PP_NP/seed_target.data -out out_pair_fn.mat 

# Chemical similarity predictions
Rscript ../scripts/GeneralDrugPairPredictMachine.r target.data.rf_model out_pair_fn.mat out_pair_fn.score	# test_chem seed_chem score 
```

## Contributors
Keunwan Park (keunwan@kist.re.kr)

