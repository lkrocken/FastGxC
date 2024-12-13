# FastGxC
Computationally efficient and statistically powerful software for detecting context-specific eQTL effects in multi-context genomic studies. 

Preprint available on [BioRxiv](https://www.biorxiv.org/content/10.1101/2021.06.17.448889v1) 

Extended data with FastGxC results on GTEx, OneK1K, and CLUES cohorts can be found [here](https://zenodo.org/record/5015123#.YNJ1WpNKjOR)

# Package Installation and Dependencies
FastGxC is an R package that can be loaded and used in any R environment. 
In order for FastGxC to run properly, the following packages are required and should be installed in R prior to using any FastGxC functions.
```
library(devtools)
library(dplyr)
library(data.table)
library(MatrixEQTL)
library(mvtnorm)
library(reshape2)
library(magrittr)
library(TreeQTL)

```

Once all dependencies are installed and loaded you can install FastGxC using:
```
  devtools::install_github("lkrocken/FastGxC")
```
Once FastGxC is installed, load all functions using 
```
library(FastGxC)
```

# Simulate toy data

To run a toy example, generate simulated data by running the following code in R:
```
  data_dir_sin = "~/simulations/"
  sim_scenario = "single_context_het"
  simulate_data(data_dir = data_dir_sim, sim_scenario = sim_scenario)
```

** Note: running the code above simulates data with default parameters (300 individuals, 10,000 SNPs, 100 genes, and 50 contexts without missing data. The total list of parameters for ```simulate_data()``` can be found at the bottom of this page.)

Running the above code will generate and save the following files in the data_dir:
(1) {sim_scenario}_SNPs.txt: SNP genotype data for 10,000 SNPs and 300 individuals (individual IDs as columns and SNP IDs as rows)
(2) {sim_scenario}_snpsloc.txt: location information of the 10,000 simulated SNPs (MatrixEQTL input format)
(3) {sim_scenario}_geneloc.txt: location information of the 100 simulated genes (MatrixEQTL input format)
(4) {sun_scenario}_simulated_expression.txt: gene expression data for the 300 simulated individuals across 100 genes and 50 contexts 

# Running FastGxC

FastGxC works in two steps. In the first step, expression is decomposed into shared and context-specific components. In the second step, eQTLs are separately mapped on these components.

*Step 1 - Decomposition:* For each individual, decompose the phenotype of interest (e.g. gene expression) across C contexts (e.g. tissues or cell-types) into one context-shared and C context-specific components using the ```decomposition_step()``` function. 
This function takes as imput a file with gene expression data for all individuals, genes, and contexts (see output of ```simulate_data()``` for the right format) and outputs one file with context-shared expression (context_shared_expression.txt) and C files with expression specific to each context (CONTEXT_NAME_specific_expression.txt). 

The following code example demonstrates how to use this function with the data we just simulated above.
  
  ```
  exp_mat_filename = "~/simulations/single_context_het_simulated_expression.txt"
  data_dir_decomp = "~/example_output_single_context_het/"
  decomposition_step(exp_mat_filename, data_dir_decomp)
  ```

*Step 2 - eQTL mapping:* FastGxC estimates genetic effects on the context-shared component and each of the C context-specific components separately using simple linear models. Note: Here we use the R package [MatrixEQTL](http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/) but any other software that can perform quick linear regression can be used (e.g. [FastQTL](http://fastqtl.sourceforge.net/) or [tensorqtl](https://github.com/broadinstitute/tensorqtl)). FastGxC implements eQTL mapping using its ```eQTL_mapping_step()``` function.

This function take as input data needed to run MatrixEQTL and outputs eQTL summary statistics in the MatrixEQTL format. In the end, you should have one file with summary statistics for shared eQTL and C files with summary statistics for each context C. 

Below is a code example to map context-specific eQTLs using the decomposed simulated data from above.
```  
out_dir = "~/example_output_single_context_het/"
input_dir = "~/simulations/"
context_names = paste0("context", seq(1,50))

for context in context_names; do
    SNP_file_name = paste0(input_dir, "single_context_het_SNPs.txt")
    snps_location_file_name = paste0(input_dir, "single_context_het_snpsloc.txt")
    expression_file_name = paste0(out_dir, context, "_specific_expression.txt")
    gene_location_file_name = paste0(input_dir, "single_context_het_geneloc.txt")
    context = context
    shared_specific = "specific"
    out_dir = 

    eQTL_mapping_step(, paste0(input_dir, )) $project_directory SNPs.txt snpsloc.txt context$i\_specific\_expression.txt geneloc.txt  context$i\_specific\_eQTLs.txt specific
 done
```
Map context-shared eQTLs
```
  Rscript run_MatrixEQTL.R $project_directory SNPs.txt snpsloc.txt context_shared_expression.txt geneloc.txt  "context_shared_eQTLs.txt" shared
```

# Multiple testing adjustment

Please download the required R packages inside `run_TreeQTL.R` before running TreeQTL. 

To adjust for multiple testing across all contexts, genes, and genetic variants tested you can use the hierarchical FDR procedures implemented in the R package [TreeQTL](http://bioinformatics.org/treeqtl/). 

TreeQTL requires that you run MatrixEQTL to do eQTL mapping (see step 2 above). If you used another eQTL mapping softwares, please make sure the output is in the format required by TreeQTL. You can also replace TreeQTL with other methods, e.g. [mashR](https://github.com/stephenslab/mashr), which can also lead to a considerable increase in power. 

Map specific-eGenes, i.e., genes with at least one context-specific eQT  
```  
project_directory=your_project_directory
Rscript run_TreeQTL.R $project_directory specific
```

Map shared-eGenes, i.e., genes with at least one context-shared eQT  
```  
project_directory=your_project_directory
Rscript run_TreeQTL.R $project_directory shared
```

This script take as input data needed to run TreeQTL and outputs shared and specific eGenes (two files) and eAssociation (C+1 files) summary statistics in the TreeQTL format. 


