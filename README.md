# omicsPattern User Manual #

omicsPattern is tool for quality control, statistical analysis, and visualization of omics data. 
**`omicsPattern` is currently under development and test and we will regularly release the documentation and tutorials**


**Citation:**


**For installation and a quick demo, read the [omicsPattern User Manual](https://github.com/broadinstitute/omicsPattern/wiki).**


## Installation ##

### Install omicsPattern in RStudio ###
1. Install devtools : 
    * ``> install.packages('devtools')``
    * ``>library(devtools)``
2. Install omicsPattern (and also all dependencies from CRAN): 
    * ``install.packages(c('dplyr','pbapply', 'lme4', 'lmerTest', 'car', 'cplm', 'pscl', 'logging', 'ggrepel', 'gridExtra', 'future', 'cowplot'), repos='http://cran.r-project.org')``
    * ``> devtools::install_github('broadinstitute/omicsPattern', force = TRUE) ``

### Input files (netome format) ###
[omicsPattern demo](https://github.com/broadinstitute/m2path/tree/master/demo)

Comming soon

## Load the library ##

```
library(omicsPattern)

```
## R script neteome_process.R to use and call the common functions ##

```
# load netome format data
loaded_data <- omicsPattern::load_data()
# check the wiki for detailed parameters

# explanatory visualization
pcoa_plots <- omicsPattern::ordplots(data = loaded_data$data, metadata = loaded_data$sample_metadata, output = output_path, outputname = 'pcoa', method = 'pcoa')

```



## Functions ##
omicsPattern support statistical analyses with visualizations. Here we dicuss several common functions between all studies and for more details on these functions and other functions check out our Wiki pages.

### Load the file of your metabolite profiles  ###


```
# load the library
library(omicsPattern)

# call the load_data function
loaded_data <- omicsPattern:::load_data(input=/path-to-file/filename.xlsx, type='known', sheet = 1, name = 'Metabolite')

data <- loaded_data$data

# ensure all data are stored as numeric
data <- omicsPattern:::numeric_dataframe(data)

# sample info
sample_info <- loaded_data$sample_metadata

# feature info ( e.g. m/z and RT)
features_info <- loaded_data$feature_metadata

```

parameters:

* input: is a excel file in `format`

* type: can be `known` for charterized metabolites with a name  or `all` for all measured features including known and uncharxterized metabolites

* sheet: defualt is `1`, to read first tab of the excel sheet of the excel files but user can use different tab number. 

* ID: is a word to use identifier for features (metabolites), the options are `Meatbolite`, `HMDB_ID`, and `Compound_ID`. 

## Output files ##

Comming soon

## [User Manual! (in the wiki)](https://github.com/omicsEye/omicsPattern/wiki) ##
