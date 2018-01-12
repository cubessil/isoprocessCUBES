# isoprocessCUBES

This is a package that will provide you with the functions used in CUBES EA and Gasbench data processing, i.e. it is a replacement for the `custom.functions.R` (which you will no longer need). To install this package, run the following from your RStudio console:

```r
install.packages("devtools") # only needed once to install the devtools package
devtools::install_github("KopfLab/isoprocessCUBES")
```

To use this package in your processing pipeline, remove the `source(file.path(base.directory, "custom.functions.R"))` and replace it simply with `library(isoprocessCUBES)` (you will also NOT need the `base.directory` variable anymore). A useful line to include outside your code chunks to make sure you know which version of **isoprocessCUBES** you have been using is the following:

```
Document knitted with isoprocessCUBES version `r packageVersion("isoprocessCUBES")`.
```

#### Examples

To download the latest example files, click on the links below and copy the contents to a new RMarkdown file (or save as csv):

 - [Gasbench corrections](https://github.com/KopfLab/isoprocessCUBES/raw/master/examples/V4_Gasbench_corrections_Master_2018Jan11.Rmd)
 - [Gasbench data](https://github.com/KopfLab/isoprocessCUBES/raw/master/examples/Gasbench_data_template.csv)
 
 
 
