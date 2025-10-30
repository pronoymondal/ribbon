##Dataset
The dataset we consider in this example consists of single-cell RNA-seq dataset from Rhesus monkey PBMC cells. The dataset can be downloaded from this [link](https://www.10xgenomics.com/datasets/2500_Rhesus_Monkey_PBMCs_Singleplex_5p_gem-x_Universal)

##Installation
The R package can be installed using the following command:
```
install.packages("devtools")
library(devtools)
install_github('pronoymondal/ribbon',ref='HEAD')
```

##Load the dataset
```
data.dir <- paste0(getwd(),"2.5k_rhesus_monkey/")

library(Seurat)
data <- Read10X_h5(filename = paste0(data.dir,"2500_Rhesus_Monkey_PBMCs_Singleplex_5p_gem-x_Universal_2500_Rhesus_Monkey_PBMCs_Singleplex_5p_gem-x_Universal_count_sample_filtered_feature_bc_matrix.h5"))       
```

##Testing for bimodality
#Using Robertson's criterion:

