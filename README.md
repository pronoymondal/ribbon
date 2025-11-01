## Dataset

The dataset we consider in this example consists of single-cell RNA-seq dataset from Rhesus monkey PBMC cells. The dataset can be downloaded from this [link](https://www.10xgenomics.com/datasets/2500_Rhesus_Monkey_PBMCs_Singleplex_5p_gem-x_Universal)

## Installation

The R package can be installed using the following command:
```
install.packages("devtools")
library(devtools)
install_github('pronoymondal/ribbon',ref='HEAD')
```

## Load the dataset
```
data.dir <- paste0(getwd(),"2.5k_rhesus_monkey/")

library(Seurat)
data <- Read10X_h5(filename = paste0(data.dir,"2500_Rhesus_Monkey_PBMCs_Singleplex_5p_gem-x_Universal_2500_Rhesus_Monkey_PBMCs_Singleplex_5p_gem-x_Universal_count_sample_filtered_feature_bc_matrix.h5"))       
```
## Normalize the dataset
```
cell_sums <- colSums(data)
cpm_data <- t(t(data) / cell_sums) * 1e6
```


## Testing for bimodality
### Using Robertson's criterion:
```
bimodal_r <- check_bimodality_robertson(cpm_data)
```
### Using Holzman and Volmer's criterion:
```
bimodal_hv <- check_bimodality(cpm_data)
```

## Estimate parameters for the bimodal distribution
```
bimodal_parameters <- estimate_bimodal(cpm_data[(bimodal_hv$indd==1),],mu1[(bimodal_hv$indd==1)],mu2[(bimodal_hv$indd==1)],sig1[(bimodal_hv$indd==1)],sig2[(bimodal_hv$indd==1)],pi[(bimodal_hv$indd==1)],D[(bimodal_hv$indd==1),(bimodal_hv$indd==1)])
```
### Output
bimodal_parameters$indd = Index denoting whether the expression of a given gene is bimodal or not, 1 denoting bimodal, 0 denoting unimodal.
bimodal_parameters$mu1 = Mean parameter of the lower mode
bimodal_parameters$mu2 = Mean parameter of the larger mode
bimodal_parameters$sd1 = Standard deviation parameter of the lower mode
bimodal_parameters$sd2 = Standard deviation parameter of the larger mode
## Estimate parameters for the unimodal distribution
```
unimodal_parameters <- estimate_bimodal(cpm_data[(bimodal_hv$indd==0),])
```
## Estimate paramaters for the unimodal and the bimodal distribution manually
```
params <- ribbon_estimate(cpm_data)
```
## Calculate p-values using RIBBON bimodal
```
bimodal_stats <- ribbon_de_bimodal(t(cpm_data[,1:1280]),t(cpm_data[,1281:2560]))
```
## Calculate p-values using RIBBON unimodal
```
unimodal_stats <- ribbon_de_unimodal(t(cpm_data[,1:1280]),t(cpm_data[,1281:2560]))
```






