# omicsArt User Manual

omicsArt is tool for quality control, statistical analysis, and visualization of omics data. **`omicsArt` is currently under development and test and we will regularly release the documentation and tutorials**

**Citation:**

Rahnavard, A. *omicsArt: omics pattern discovery by visualization*. Version 1.0.0.0, <https://github.com/omicsEye/omicsArt> (2021).

**For installation and a quick demo, read the [omicsArt User Manual](https://github.com/omicsEye/omicsArt/wiki).**

# omicsArt user manual

## Contents

-   [Features](#features)
-   [omicsArt](#omicsArt)
    -   [omicsArt approach](#omicsArt-approach)
    -   [Requirements](#requirements)
    -   [Installation](#installation)
-   [Orindation plots omicsArt](#getting-started-with-omeClust)
    -   [t-sne plot](#t-sne-plot)
    -   [PCoA plot](#t-sne-plot)
-   [Microbial community diversity](#Microbial-community-diversity)
    -   [Overall tests](#overall-test)
    -   [Diversity metadata associations and plots](#Diversity-metadata-associations-and-plots)
-   [Multi-barplot form Tweedieverse results](#Multi-barplot-form-Tweedieverse-results)

## Installation {#installation}

### Install omicsArt in RStudio

1.  Install devtools :
    -   `> install.packages('devtools')`
    -   `>library(devtools)`
2.  Install omicsArt (and also all dependencies from CRAN):
    -   `install.packages(c('dplyr','pbapply', 'lme4', 'lmerTest', 'car', 'cplm', 'pscl', 'logging', 'ggrepel', 'gridExtra', 'future', 'cowplot',      'Hmisc', 'MultiDataSet', 'TSP', 'htmlTable', 'igraph', 'insight',      'lubridate', 'mgcv', 'mvtnorm', 'optparse', 'parameters', 'pillar',      'pkgload', 'plotly', 'rlang', 'rvest', 'seriation', 'usethis', 'viridis',    'signal', 'tsne', 'openxlsx', 'readxl', 'xfun', 'yulab.utils', "labdsv", "seriation","diffusionMap"), repos='http://cran.r-project.org')`

        `if (!requireNamespace("BiocManager", quietly = TRUE))              install.packages("BiocManager")   BiocManager::install("ropls")`

    -   `devtools::install_github('omicsEye/omicsArt', force = TRUE)`

### Input files (omics format)

[omicsArt demo](https://github.com/omicsEye/omicsArt/tree/master/demo)

Comming soon

## Load the library

    library(omicsArt)

## R script neteome_process.R to use and call the common functions

    # load omics format data
    loaded_data <- omicsArt::load_data()
    # check the wiki for detailed parameters

    # explanatory visualization
    pcoa_plots <- omicsArt::ordplots(data = loaded_data$data, metadata = loaded_data$sample_metadata, output = output_path, outputname = 'pcoa', method = 'pcoa')

## Functions

omicsArt support statistical analyses with visualizations. Here we dicuss several common functions between all studies and for more details on these functions and other functions check out our Wiki pages.

### Load the file of your metabolite profiles

    # load the library
    library(omicsArt)

    # call the load_data function
    loaded_data <- omicsArt:::load_data(input=/path-to-file/filename.xlsx, type='known', sheet = 1, name = 'Metabolite')

    data <- loaded_data$data

    # ensure all data are stored as numeric
    data <- omicsArt:::numeric_dataframe(data)

    # sample info
    sample_info <- loaded_data$sample_metadata

    # feature info ( e.g. m/z and RT)
    features_info <- loaded_data$feature_metadata

parameters:

-   input: is a excel file in `format`

-   type: can be `known` for characterized metabolites with a name or `all` for all measured features including known and uncharxterized metabolites

-   sheet: defualt is `1`, to read first tab of the excel sheet of the excel files but user can use different tab number.

-   ID: is a word to use identifier for features (metabolites), the options are `Meatbolite`, `HMDB_ID`, and `Compound_ID`.

## Diversity test and visualization

-   data: includes taxa profiles n\*m where n is number of observations (samples) and m are number features (e.g. microbila species or taxa).

-   metadata: includes n\*p where n is number of observations (samples) and p are number of metadat or information about samples (e.g. age, sex, and health status).

<!-- -->

    infant_results <- omicsArt::diversity_test(data, metadat)
    infant_alpha_diversity_data <- infant_results$alpha_diversity_data
    infant_alpha_diversity_test <- infant_results$alpha_diversity_test
    infant_alpha_diversity_plots <- infant_results$diversity_test_plots
    infant_overall_diversity_barplot <- infant_results$overall_diversity_barplot


    pdf(
      paste('analysis', '/alpaha_diversity.pdf', sep = ''),
      width = 2.4,
      height = 2.25,
      onefile = TRUE
    )

    for (meta in unique(colnames(metadata))){
      tryCatch({
        stdout <-
          capture.output(print(alpha_diversity_plots[[meta]]), type = "message")
      }, error = function(e) {
        print(meta)
        print(paste('error:', e))
      })
    }

## Multi-barplot form Tweedieverse results {#multi-barplot-form-tweedieverse-results}

    library(tidyr)
    library(dplyr)
    library(reshape2)
    library(deepath)
    #setting the working directory
    setwd("~/Projects/")

    number_of_sig_to_keep <- 20
    sig_threshold <- 0.05

    ## read metabolites
    metabolites_Tweedieverse <- read.delim(
      "analysis/meatbolites_Tweedieverse/all_results.tsv",
      sep = '\t',
      header = T,
      fill = F,
      comment.char = "" ,
      check.names = F,
      #row.names = NA
    )
    metabolites_score_data_severe <- metabolites_Tweedieverse[metabolites_Tweedieverse$metadata=="Group" & metabolites_Tweedieverse$value=="Severe" ,]
    rownames(metabolites_score_data_severe) <- metabolites_score_data_severe$feature

    metabolites_score_data_non_severe <- metabolites_Tweedieverse[metabolites_Tweedieverse$metadata=="Group" & metabolites_Tweedieverse$value=="non-Severe" ,]
    rownames(metabolites_score_data_non_severe) <- metabolites_score_data_non_severe$feature

    metabolites_score_data_non_covid <- metabolites_Tweedieverse[metabolites_Tweedieverse$metadata=="Group" & metabolites_Tweedieverse$value=="non-COVID-19" ,]
    rownames(metabolites_score_data_non_covid) <- metabolites_score_data_non_covid$feature



    # use score_data_severe is reference
    order_sig <- rownames(metabolites_score_data_severe)[1:number_of_sig_to_keep]
    metabolites_score_data_severe <- metabolites_score_data_severe[order_sig,]
    metabolites_score_data_severe<- metabolites_score_data_severe[order(metabolites_score_data_severe$coef),]
    order_sig <- rownames(metabolites_score_data_severe)
    metabolites_score_data_severe <- within(metabolites_score_data_severe,
                                            feature <- factor(feature,
                                                              levels=order_sig))
    metabolites_score_data_non_severe <- metabolites_score_data_non_severe[rownames(metabolites_score_data_severe),]
    metabolites_score_data_non_severe <- within(metabolites_score_data_non_severe,
                                                feature <- factor(feature,
                                                                  levels=order_sig))


    metabolites_score_data_non_covid <- metabolites_score_data_non_covid[rownames(metabolites_score_data_severe),]
    metabolites_score_data_non_covid <- within(metabolites_score_data_non_covid,
                                               feature <- factor(feature,
                                                                 levels=order_sig))

    metabolites_severe_temp_diff_bar <- diff_bar_plot(metabolites_score_data_severe, threshold = sig_threshold, pvalue_col = "pval",  method = "none",
                                                      fdr ="qval", orderby = NA, x_label = 'Coefficient', y_label = '')
    metabolites_non_severe_temp_diff_bar <- diff_bar_plot(metabolites_score_data_non_severe, threshold = sig_threshold, pvalue_col = "pval",  method = "none",
                                                          fdr ="qval", orderby = NA, x_label = 'Coefficient', y_label = '')
    metabolites_non_covid_temp_diff_bar <- diff_bar_plot(metabolites_score_data_non_covid, threshold = sig_threshold, pvalue_col = "pval",  method = "none",
                                                         fdr ="qval", orderby = NA, x_label = 'Coefficient', y_label = '')



    ## read association
    box_association <- readRDS("/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/COVID-Omics/analysis/meatbolites_Tweedieverse/figures/Group_gg_associations.RDS")
    ## do plots

    fig2_metabolites <- ggdraw() +
      draw_plot(metabolites_severe_temp_diff_bar,
                x = 0, y = .47, width = .55, height = .53) +
      draw_plot(metabolites_non_severe_temp_diff_bar + theme(axis.title.y = element_blank(),
                                                             axis.text.y = element_blank(),
                                                             axis.ticks.y = element_blank(),
                                                             axis.line.y = element_blank()),
                x = .55, y = .47, width = .225, height = .53) +
      draw_plot(metabolites_non_covid_temp_diff_bar + theme(axis.title.y = element_blank(),
                                                            axis.text.y = element_blank(),
                                                            axis.ticks.y = element_blank(),
                                                            axis.line.y = element_blank()),
                x = .775, y = .47, width = .225, height = .53) +
      draw_plot(box_association[[11]] + theme(
        axis.title.x = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.y = element_text(size = 5)), x = 0, y = 0, width = .25, height = .45) +
      draw_plot(box_association[[52]] + theme(
        axis.title.x = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.y = element_text(size = 5)), x = .25, y = 0, width = .25, height = .45) +
      draw_plot(box_association[[139]] + theme(
        axis.title.x = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.y = element_text(size = 5)), x = .5, y = 0, width = .25, height = .45) +
      draw_plot(box_association[[2]] + theme(
        axis.title.x = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.y = element_text(size = 5)), x = .75, y = 0, width = .25, height = .45) +

      draw_plot_label((label = c("a",  "Severe", "non-Severe", "non-COVID", "b", "c", "d", "e")),
                      size = 7,x = c(0, .28, .53, .76, 0, .25, .5, .75), y = c(1, 1, 1, 1, 0.47, 0.47, 0.47, 0.47))
    fig3_metabolites

    ggsave(filename = 'figures/fig3/fig#_barplot.pdf', plot=fig2_metabolites, width = 183, height = 110, units = "mm", dpi = 350)
    ggsave(filename = 'figures/fig3/fig#_barplot.pdf', plot=fig2_metabolites, width = 183, height = 110, units = "mm", dpi = 350)

## [User Manual! (in the wiki)](https://github.com/omicsEye/omicsArt/wiki)
