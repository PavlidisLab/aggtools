## Aggregating single cell coexpression

Basic functionality

```
devtools::install_github("PavlidisLab/aggtools")
library(aggtools)

pc_df <- read.delim("/home/amorin/Data/Metadata/refseq_select_mm10.tsv")
dat <- readRDS("/space/scratch/amorin/TR_singlecell/GSE145172/GSE145172_clean_mat_and_meta_CPM.RDS")


result <- aggr_coexpr_single_dataset(mat = dat$Mat, 
                                     meta = dat$Meta,
                                     pc_df = pc_df,
                                     cor_method = "pearson",
                                     agg_method = "FZ",
                                     verbose = TRUE)

```

NOTE: Current implementation has strict assumptions about the input data:

1) pc_df is a data.frame with a "Symbol" column whose elements are unique and exactly match the input count matrix rownames
2) meta is a data.frame that has an "ID" and "Cell_type" columns (the latter is the grouping factor for aggregation, and thus could be use for things like different subjects or conditions)
3) mat is a sparse matrix of class "dgCMatrix" (see package Matrix::)
4) all column names of mat are found in the "ID" column of meta

The multi dataset function makes these further assumptions:

5) input_df is a data.frame of single cell datasets containing columns "ID", "Cell_type", and "Path"
6) The data paths specified in input_df lead to .RDS objects of lists with two named elements: "Meta" for the metadata and "Mat" for the sparse count matrix
