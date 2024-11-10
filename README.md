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

Current implementation has strict assumptions about the input data:

1)  `pc_df` is a data.frame with a "Symbol" column whose elements are unique and contain all gene/rownames of `mat`. Only the elements in `pc_df$Symbol` are kept: `mat` will be subset to these genes. An error will be thrown if `pc_df$Symbol` has genes that are not found in the rownames of `mat`.

2)  `meta` is a data.frame that has "ID" and "Cell_type" columns -- the latter is the grouping factor for aggregation to subset for IDs/cells in `mat`. Currently if you want to subset by patient/condition/some other factor other than cell type, you would have to name the corresponding column in `meta` to "Cell_type".

3)  `mat` is a gene (rows) by cell (columns) sparse matrix of class "dgCMatrix" (see package Matrix)

4)  All column names (cell IDs) of `mat` are found in the "ID" column of `meta`.

The multi dataset function makes these further assumptions:

5)  `input_df` is a data.frame of single cell datasets containing columns "ID", "Cell_type", and "Path".

6)  The data paths specified in `input_df$Path` lead to .RDS objects of lists with two named elements: "Meta" for the metadata and "Mat" for the sparse count matrix.
