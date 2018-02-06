# Example analysis scripts of "PhenotypeSeeker"

After the successful installation of the PhenotypeSeeker you can launch the example analysis scripts to aquire better intuition about the program's workflows.

## An example of "PhenotypeSeeker modeling"

To launch the example analysis script of "PhenotypeSeeker modeling", open the command-line, change to the "PhenotypeSeeker/example" directory and type in:
```
./test_PS_modeling.sh
```
The process starts with downloading the folder containing 30 *C.difficile* genomes (174 MB) originating from European Nucleotide Archive [EMBL:PRJEB11776 ((http://www.ebi.ac.uk/ena/data/view/PRJEB11776)] and inputfile for "PhenotypeSeeker modeling" containing the binary phenotypes of azithromycin resistance for these genomes, adapted from Lees et al. (2016. Sequence element enrichment analysis to determine the genetic basis of bacterial phenotypes. Nature Communications, 7, 12797. https://doi.org/10.1038/ncomms12797). The folder will be unpacked and after that the "PhenotypeSeeker modeling" with the default parameters is launched automatically:
```
phenotypeseeker modeling PS_modeling_example_files/data.pheno
```
The "PhenotypeSeeker modeling" tests every k-mer presented in genomes for association with 0 or 1 phenotype. The k-mers having chi-square test p-value lower than 0.05 (default) are saved in **k-mers_filtered_by_pvalue_Azithromycin.txt**. 1000 (default) lowest p-valued k-mers are selected as features for logistic regression phenotype prediction model generation. The model is saved in **log_reg_model_Azithromycin.pkl**. All the information about the conducted regression analysis is saved in **summary_of_log_reg_analysis_Azithromycin.txt**. As it shows, the generated models' accuracy on the 8 test samples (~25% of input samples by default) is 1.0. The k-mers used in model generation are saved in **k-mers_and_coefficients_in_log_reg_model_Azithromycin.txt** with their coefficients in model. The presence of k-mers with positive coefficients (from *ermB* gene) predict phenotype 1 and k-mers with negative coefficients predict phenotype 0.

The example analysis takes approximately 20 minutes and consumes approximately 1 GB of memory.

All the mentioned outputfiles are saved in "PhenotypeSeeker/example" directory.

## An example of "PhenotypeSeeker prediction"

Before launching the "Phenotype prediction" example analysis script, the phenotype prediction model has to be created with "PhenotypeSeeker modeling".

To launch the example analysis script of "PhenotypeSeeker prediction", open the command-line, change to the "PhenotypeSeeker/example" directory and type in:
```
./test_PS_prediction.sh
```
The process starts with downloading the folder containing 10 *C.difficile* genomes (42 MB) originating from European Nucleotide Archive [EMBL:PRJEB11776 ((http://www.ebi.ac.uk/ena/data/view/PRJEB11776)] and inputfiles for "PhenotypeSeeker prediction". The folder will be unpacked and after that the "PhenotypeSeeker prediction" with the default parameters is launched automatically:
```
phenotypeseeker prediction 
```
The "PhenotypeSeeker prediction" detects the presence or absence of the model specific k-mers in every input sample and uses this data to predict the phenotypes for samples.

The results of predictions are saved into **"PhenotypeSeeker/example/predictions_Azithromycin.txt"**. In this example the model predicts that 5 of input samples have phenotype 0 and 5 of inputsamples have phenotype 1.

The example analysis takes approximately 2 minutes and consumes approximately 2 GB of memory.
