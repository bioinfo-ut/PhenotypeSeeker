# An example of "PhenotypeSeeker modeling"
 
After the successful installation of the PhenotypeSeeker you can launch the example script to aquire better intuition about the program's workflow.

In order to do so, open the command-line, change to the "PhenotypeSeeker/example" directory and type in:
```
./test_PS_modeling.sh
```
The process starts with downloading the folder containing 30 *C.difficile* genomes (174 MB) originating from European Nucleotide Archive [EMBL:PRJEB11776 ((http://www.ebi.ac.uk/ena/data/view/PRJEB11776)] and inputfile for "PhenotypeSeeker modeling" containing the binary phenotypes of azithromycin resistance for these genomes adapted from Lees et al. (2016. Sequence element enrichment analysis to determine the genetic basis of bacterial phenotypes. Nature Communications, 7, 12797. https://doi.org/10.1038/ncomms12797). The folder will be unpacked and after that the "PhenotypeSeeker modeling" with the default parameters is launched automatically:
```
phenotypeseeker modeling PS_modeling_example_files/data.pheno
```
The PhenotypeSeeker tests every k-mer presented in genomes for association with 0 or 1 phenotype. The k-mers having chi-square test p-value lower than 0.05 (default) are saved in **k-mers_filtered_by_pvalue_Azithromycin.txt**. 1000 (default) lowest p-valued k-mers are selected as features for logistic regression phenotype prediction model generation. The model is saved in **log_reg_model_Azithromycin.pkl**. All the information about the conducted regression analysis is saved in **summary_of_log_reg_analysis_Azithromycin.txt**. As it shows, the generated models' accuracy on the 8 test samples (~25% of input samples by default) is 1.0. The k-mers used in model generation are saved in **k-mers_and_coefficients_in_log_reg_model_Azithromycin.txt** with their coefficients in model. The presence of k-mers with positive coefficients (from *ermB* gene) predict phenotype 1 and k-mers with negative coefficients predict phenotype 0.

The example analysis takes approximately 20 minutes and consumes approximately 1 GB of memory.

All the mentioned outputfiles of PS are saved in "PhenotypeSeeker/example" directory.

