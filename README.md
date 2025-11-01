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
## Using Robertson's criterion:
```
bimodal_r <- check_bimodality_robertson(cpm_data)
```
### Output
bimodal_r$indd = Index denoting whether the expression of a given gene is bimodal or not, 1 denoting bimodal, 0 denoting unimodal.<br />
bimodal_r$mu1 = Mean parameter of the lower mode.<br />
bimodal_r$mu2 = Mean parameter of the larger mode.<br />
bimodal_r$sig1 = Standard deviation parameter of the lower mode.<br />
bimodal_r$sig2 = Standard deviation parameter of the larger mode.<br />
bimodal_r$pi = Mixing proportion parameter.<br />
bimodal_r$D = An array with dimension nCells x nGenes. The (i.j)-th entry denotes the probability that the i-th cell has expression in the larger mode for gene j.
## Using Holzman and Volmer's criterion:
```
bimodal_hv <- check_bimodality(cpm_data)
```
### Output
bimodal_hv = Index denoting whether the expression of a given gene is bimodal or not, 1 denoting bimodal, 0 denoting unimodal
## Estimate parameters for the bimodal distribution
```
bimodal_parameters <- estimate_bimodal(cpm_data[(bimodal_hv$indd==1),],mu1[(bimodal_hv$indd==1)],mu2[(bimodal_hv$indd==1)],sig1[(bimodal_hv$indd==1)],sig2[(bimodal_hv$indd==1)],pi[(bimodal_hv$indd==1)],D[(bimodal_hv$indd==1),(bimodal_hv$indd==1)])
```
### Output
bimodal_parameters$mu1 = Mean parameter of the lower mode.<br />
bimodal_parameters$mu2 = Mean parameter of the larger mode.<br />
bimodal_parameters$tau1 = $\tau_1^2$ parameter of RIBBON bimodal.<br />
bimodal_parameters$tau2 = $\tau_2^2$ parameter of RIBBON bimodal.<br />
bimodal_parameters$c = $c$ parameter of RIBBON bimodal.<br />
bimodal_parameters$alpha1 = $\alpha_1$ parameter of RIBBON bimodal.<br />
bimodal_parameters$alpha2 = $\alpha_2$ parameter of RIBBON bimodal.<br />
bimodal_parameters$sig1 = Standard deviation parameter of the lower mode.<br />
bimodal_parameters$sig2 = Standard deviation parameter of the larger mode.<br />
bimodal_parameters$pi = $\pi$ parameter of RIBBON bimodal.<br />
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






