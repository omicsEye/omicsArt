
# Stand-alone utility functions

#' PERMANOVA with repeat measure-aware permutations. Block sizes are allowed to
#' differ.
#'
#' @param D An N-by-N distance matrix (must be a \code{dist} object).
#' @param permute_within Data frame with N rows containing metadata to test per
#' sample
#' @param blocks A length-N vector containing the block structure.
#' @param block_data Data frame with per-block metadata. If \code{blocks} is
#' numeric, \code{block_data}'s rows must match those indices. If \code{blocks}
#' is a \code{factor}, then the row ordering must match the factor levels. If
#' \code{blocks} is a character vector, \code{block_data} must have rownames
#' matching the contents of \code{blocks}.
#' @param permutations Number of permutations to test
#' @param metadata_order Order of the metadata in the model. If not given, this
#' is assumed to be within-block metadata first, followed by block metadata, in
#' the order given in \code{permute_within} and \code{block_data}.
#' @return Same structure as \code{adonis}, with p-values recalculated based on
#' permutations that are aware of the block structure of the data.
#' Metadata in \code{permute_within} are permuted within blocks, whereas
#' metadata in \code{block_data} are first permuted across blocks, and then
#' assigned to samples according to the block structure.
PERMANOVA_repeat_measures <- function(
    D, permute_within, blocks = NULL, block_data, permutations=999,
    metadata_order = c(names(permute_within), names(block_data)),
    na.rm=F) {

    # Make sure D is a dist object
    if (class(D) != "dist") {
        stop("D must be a dist object")
    }

    # Default to free permutations if blocks is not given
    if (!missing(block_data) && is.null(blocks)) {
        stop("blocks must be given if block_data is present")
    } else if (is.null(blocks)) {
        blocks <- rep(1, nrow(permute_within))
        block_data <- as.data.frame(matrix(0, nrow=1, ncol=0))
    } else if (length(unique(blocks)) == 1) {
        warning("blocks only contains one unique value")
    }

    # Ensure no metadata overlap between permute_within and block_data
    if (length(intersect(names(permute_within), names(block_data))) > 0) {
        stop("metadata is repeated across permute_within and block_data")
    }

    # Ensure that metadata_order only contains stuff in permute_within and block_data
    if(length(setdiff(metadata_order, union(names(permute_within), names(block_data)))) > 0) {
        stop("metadata_order contains metadata not in permute_within and block_data")
    }

    # Ensure that the data in permute_within matches that in dist
    ord <- rownames(as.matrix(D))
    if (length(ord) != nrow(permute_within) || length(blocks) != length(ord)) {
        stop("blocks, permute_within, and D are not the same size")
    }
    if (is.null(rownames(permute_within))) {
        warning("permute_within has no rownames - can't verify sample orders")
    } else if (!all(ord == rownames(permute_within))) {
        stop("rownames do not match between permute_within and D")
    }

    # Ensure matching between blocks and block_data
    if (any(is.na(blocks))) {
        stop("NAs are not allowed in blocks")
    }
    if (is.factor(blocks)) {
        if (length(levels(blocks)) != nrow(block_data)) {
            stop("block_data does not have as many rows as blocks has levels")
        }
        if (!is.null(rownames(block_data)) && any(rownames(block_data) != levels(blocks))) {
            stop("block_data rownames does not match the levels of blocks")
        }
        # Discard level information
        blocks <- as.numeric(blocks)
    } else if (is.numeric(blocks)) {
        if (blocks < 1 || max(blocks) > nrow(block_data)) {
            stop("Numeric blocks has indices out of range")
        }
    } else if (is.character(blocks)) {
        if (is.null(rownames(block_data)) || !all(blocks %in% rownames(block_data))) {
            stop("blocks does not match the rownames of block_data")
        }
        # Transform to numeric
        blocks <- match(blocks, rownames(block_data))
    } else {
        stop("blocks must be a numeric, factor, or character vector")
    }

    # Error out on NA metadata rather than allowing adonis to error out with
    # a totally nonsensical error message
    if (any(is.na(permute_within)) || any(is.na(block_data))) {
        if (na.rm) {
            n_prerm <- length(blocks)

            # Remove NAs in block_data
            hasna <- (rowSums(is.na(block_data)) > 0) | (sapply(split(rowSums(is.na(permute_within)) > 0, blocks), mean) == 1)
            block_data <- block_data[!hasna,, drop=F]
            keep <- !hasna[blocks]
            blocks <- cumsum(!hasna)[blocks]

            blocks <- blocks[keep]
            permute_within <- permute_within[keep,, drop=F]
            D <- as.matrix(D)[keep, keep]
            # block_data is not subset, as the rows with NAs are no longer referenced in blocks

            # Remove NAs in permute_within
            keep <- rowSums(is.na(permute_within)) == 0
            blocks <- blocks[keep]
            permute_within <- permute_within[keep,, drop=F]
            D <- as.dist(D[keep, keep])

            if (length(blocks) < ncol(permute_within) + ncol(block_data)) {
                stop(sprintf("After omitting samples with NAs, the number of samples (%d) is less than the number of metadata (%d)",
                             length(blocks), ncol(permute_within) + ncol(block_data)))
            } else if (length(blocks) < n_prerm * 0.5) {
                warning(sprintf("Removed %d samples with NA metadata", n_prerm - length(blocks)))
            }
        } else {
            stop("Some metadata is NA! adonis does not support any NA in the metadata")
        }
    }

    # Warn on some suspicious input
    persample <- apply(permute_within, 1, function(x)is.factor(x) && !any(duplicated(x)))
    if (any(persample)) {
        warning(sprintf("%s in permute_within has one DOF per sample.", colnames(permute_within)[which(persample)[1]]))
    }
    if (length(unique(blocks)) < nrow(block_data)) {
        warning("Not all blocks have a sample associated with them. Block permutations will still be performed over the full set of blocks - if this is not desired, subset block_data to only the blocks which appear in the data.")
    }
    if (!any(duplicated(blocks))) {
        warning("blocks contains no duplicated elements")
    }

    library(vegan)
    library(permute)

    # Test statistic from non-permuted data
    mtdat <- cbind(permute_within, block_data[blocks,,drop=F])
    ad <- adonis(D ~ ., permutations=0, data=mtdat[, metadata_order, drop=F])
    R2 <- ad$aov.tab$R2
    names(R2) <- rownames(ad$aov.tab)

    # Permutations
    nullsamples <- matrix(NA, nrow=length(R2), ncol=permutations)
    for (i in seq_len(permutations)) {
        within.i <- shuffle(nrow(permute_within), control=how(blocks=blocks))
        block.i <- sample(seq_len(nrow(block_data)))
        mtdat <- cbind(
            permute_within[within.i,,drop=F],
            block_data[block.i,,drop=F][blocks,,drop=F])
        perm.ad <- adonis(D ~ ., permutations=0, data=mtdat[, metadata_order, drop=F])

        nullsamples[,i] <- perm.ad$aov.tab$R2
    }

    # For residuals, test the other direction (i.e. p-value of all covariates)
    n <- length(R2)
    R2[n-2] <- 1 - R2[n-2]
    nullsamples[n-2,] <- 1 - nullsamples[n-2,]

    # P value calculation similar to adonis's
    exceedances <- rowSums(nullsamples > R2)
    P <- (exceedances + 1) / (permutations + 1)

    P[n] <- NA    # No p-values for "Total"
    ad$aov.tab$`Pr(>F)` <- P

    return (ad)
}

source("./common/theme_nature.r")

#' PERMANOVA results "heatmap"
#'
#' @param R2 A matrix of R2 values from \code{adonis}.
#' @param P A matrix of P-values from \code{adonis}.
#' @param fontsize The desired font size of the % R2 values in the heatmap.
#' @param FDR If \code{T}, P-values are first BH-adjusted, and significance
#' stars are shown as *** < 0.001, ** < 0.01, * < 0.05.
#' @return ggplot object.
PERMANOVA_heatmap <- function(R2, P, fontsize=6, FDR=T) {
    library(reshape2)
    library(ggplot2)
    library(RColorBrewer)

    df <- melt(R2)
    colnames(df) <- c("Feature", "Dataset", "R2")
    df$P <- melt(P)[,3]
    if (FDR) {
        df$P <- p.adjust(df$P, method="fdr")
    }
    df$varExpPct <- sprintf("%.1f%%", 100*df$R2)
    df$varExpPct[is.na(df$R2)] <- ""
    df$NAtext <- ""
    df$NAtext[is.na(df$R2)] <- "N/A"
    df$stars <- cut(df$P, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), c("***", "**", "*", ""))
    df$stars[is.na(df$R2)] <- ""
    df$R2[is.na(df$R2)] <- 0

    # Try to make a reasonable color scheme that has contrast where needed
    colors <- brewer.pal(name="Oranges", n=9)
    R2Q <- quantile(R2, c(0.25, 0.75), na.rm=T)
    labhat <- optim(par=c(0, 0), method="Nelder-Mead",
                    fn=function(lab) sum((pbeta(c(0.25, 0.75), exp(lab[1]), exp(lab[2])) - R2Q)^2))
    abhat <- exp(labhat$par)
    colorvalues <- pbeta(seq(0, 1, length=length(colors)), abhat[1], abhat[2])

    # Label colors flip to white when the color is too dark
    df$lblcolor <- ifelse(qbeta(df$R2, abhat[1], abhat[2]) < 0.85, "black", "white")

    ggp <- ggplot(data=df, aes(x=Dataset, y=Feature)) +
        geom_tile(aes(fill=R2)) +
        geom_text(aes(label=varExpPct, color=lblcolor), size=fontsize/(14/5), nudge_y=-0.15) +
        geom_text(aes(label=NAtext), color="grey", size=fontsize/(14/5), nudge_y=-0.15) +
        geom_text(aes(label=stars, color=lblcolor), fontface="bold", size=1.25*fontsize/(14/5), nudge_y=0.12) +
        scale_fill_gradientn(colors=colors, values=colorvalues, limits=c(0, 1)) +
        scale_color_manual(values=c(black="black", white="white")) +
        scale_x_discrete(expand=c(0,0)) + xlab(NULL) +
        scale_y_discrete(expand=c(0,0), position = "right", limits = rev(levels(df$Feature))) + ylab(NULL) +
        guides(color="none",
               fill=guide_colourbar(title=NULL, barheight=unit(40,"mm"), label.position = "left")) +
        theme_nature() +
        theme(axis.text.x = element_text(angle=-17, hjust=0),
              #panel.border=element_rect(),
              legend.position = "left", axis.ticks.y = element_blank())

    return (ggp)
}


# -----------------------------------------------------------------------------
#   HMP2-specific

# Generates the source of Fig. 1d as ./overview/variance_explained_max.pdf
# P-values are dumped to ./overview/variance_explained_max_p.txt
# Significance is to be merged in separately

source("./common/merge_metadata.r")
source("./common/merge_BMI.r")
source("./common/disease_activity.r")
source("./common/uc_disease_locations.r")
source("./common/biopsy_metadata.r")


# What covariates to include
covariates <- list(
    # "Name in heatmap" = c(list of metadata names to include as covariates)
    Age = "consent_age",
    Sex = "sex",
    BMI = "BMI_imputed",
    Race = "race_simplified",
    Antibiotics = "abx",
    Immunosuppressants = "immsup",
    "Bowel Surgery" = "bowelsurgery",
    "Recruitment Site" = "site_name",
    #Pilot = "pilot", # whether the sample was from the pilot
    "Biopsy Location" = "biopsy_location_ilecolrec",
    Inflammation = "is_inflamed",
    Disease = "diagnosis",
    "Disease Location (UC)" = "uc_baseline_location_nona",
    "Disease Location (CD)" = c("montreal_l123", "montreal_l4"),
    "Dysbiotic (CD)" = "active_cd",
    "Dysbiotic (UC)" = "active_uc",
    "Dysbiotic (non-IBD)" = "active_nonibd",
    "Dysbiotic" = "active",
    Subject = "subject"
)

# These covariates are treated as per-individual (i.e. labels are permuted across
# individuals rather than longitudinally)
perindividual_covariates <- c(
    "consent_age", "sex", "BMI_imputed", "race_simplified", "site_name", "diagnosis",
    "montreal_l123", "montreal_l4", "subject", "uc_baseline_location_nona"
)

# Some covariates have their R2 calculated from only a subset of samples
# (e.g. disease location only matters within the respective disease)
subsetted_variables <- list(
    # c(variable, subset/exclude, value)
    "Disease Location (CD)" = c("diagnosis", T, "CD", "montreal_l123", F, "None"),
    "Disease Location (UC)" = c("diagnosis", T, "UC", "uc_baseline_location_nona", F, "None"),
    "Dysbiotic (CD)" = c("diagnosis", T, "CD"),
    "Dysbiotic (UC)" = c("diagnosis", T, "UC"),
    "Dysbiotic (non-IBD)" = c("diagnosis", T, "nonIBD")
)

# Covariates and order to include in the final heatmap
put_in_heatmap <- c("Age", "Sex", "BMI", "Race", "Antibiotics", "Immunosuppressants",
    "Recruitment Site", "Biopsy Location", "Inflammation", "Disease",
    "Disease Location (CD)", "Disease Location (UC)", "Dysbiotic (CD)",
    "Dysbiotic (UC)", "Subject", "All")

# Initialize the heatmap contents
var_explained.separate <- matrix(NA, nrow=length(covariates)+1, ncol=0)
rownames(var_explained.separate) <- c(names(covariates), "All")
var_explained.separate <- as.data.frame(var_explained.separate)
var_explained.separate.p <- var_explained.separate


merge_metadata_for_omnibus <- function(pcl, act_lenience=0, week_offset=0) {
    # Merge in/create the expected metadata
    pcl <- pcl %>% merge_metadata(c(
            abx = "Antibiotics",
            immsup = "Immunosuppressants..e.g..oral.corticosteroids.",
            bowelsurgery = "X6..Have.you.ever.had.bowel.surgery."),
        week_offset=week_offset) %>%
        merge_biopsy_metadata
    pcl$meta$ped_disease <- factor(sprintf("%s%s", ifelse(pcl$meta$diagnosis=="nonIBD", "", ifelse(pcl$meta$adult, "Adult", "Pediatric")), pcl$meta$diagnosis))

    if (!("is_inflamed" %in% colnames(pcl$meta))) {
        pcl$meta$is_inflamed <- F
    }
    if (!("biopsy_location_ilecolrec" %in% colnames(pcl$meta))) {
        pcl$meta$biopsy_location_ilecolrec <- "None"
    }

    pcl <- merge_BMI(pcl)
    pcl$meta$BMI_imputed <- pcl$meta$BMI
    pcl$meta$BMI_imputed[is.na(pcl$meta$BMI_imputed)] <- mean(pcl$meta$BMI_imputed[!duplicated(pcl$meta$subject)], na.rm=T)

    pcl$meta$montreal_l123 <- rep(factor("None", levels=c("None", "L1", "L2", "L3")))
    pcl$meta$montreal_l123[grepl("L1", pcl$meta$baseline_montreal_location)] <- "L1"
    pcl$meta$montreal_l123[grepl("L2", pcl$meta$baseline_montreal_location)] <- "L2"
    pcl$meta$montreal_l123[grepl("L3", pcl$meta$baseline_montreal_location)] <- "L3"
    pcl$meta$montreal_l4 <- factor(ifelse(grepl("L4", pcl$meta$baseline_montreal_location), "None", "L4"), levels=c("None", "L4"))

    pcl <- merge_uc_disease_locations(pcl)
    pcl$meta$uc_baseline_location_nona <- factor(pcl$meta$uc_baseline_location, levels=c("None", levels(pcl$meta$uc_baseline_location)))
    pcl$meta$uc_baseline_location_nona[is.na(pcl$meta$uc_baseline_location_nona)] <- "None"
    #pcl$meta$baseline_uc_extent_nona <- factor(pcl$meta$baseline_uc_extent, levels=c("None", levels(pcl$meta$baseline_uc_extent)))
    #pcl$meta$baseline_uc_extent_nona[is.na(pcl$meta$baseline_uc_extent_nona)] <- "None"

    pcl <- merge_disease_activity(pcl, lenience=act_lenience, week_offset=week_offset)
    pcl$meta$active_cd <- (pcl$meta$diagnosis == "CD") & pcl$meta$active
    pcl$meta$active_cd[is.na(pcl$meta$active_cd)] <- F
    pcl$meta$active_uc <- (pcl$meta$diagnosis == "UC") & pcl$meta$active
    pcl$meta$active_uc[is.na(pcl$meta$active_uc)] <- F
    pcl$meta$active_nonibd <- (pcl$meta$diagnosis == "nonIBD") & pcl$meta$active
    pcl$meta$active_nonibd[is.na(pcl$meta$active_nonibd)] <- F

    pcl$meta$race_simplified <- pcl$meta$race
    pcl$meta$race_simplified[pcl$meta$race_simplified=="American Indian or Alaska Native"] <- "Other"

    pcl$meta$bowelsurgery[pcl$meta$bowelsurgery==""] <- NA

    return (pcl)
}

bc_omnibus_tests <- function(pcl, week_offset=0, method="bray", act_lenience=2, Nperms=4999) {
    # Main function for performing overall and per-covariate PERMANOVAs for a given dataset

    # Merge in relevant metadata
    pcl <- merge_metadata_for_omnibus(pcl, act_lenience=act_lenience, week_offset=week_offset)

    # Subset to samples with no NA's (or adonis will puke with a bizarre error message)
    covars.unl <- unlist(covariates)
    ok_covars <- colMeans(is.na(pcl$meta[,covars.unl])) < 0.5
    # NA removal is now done in PERMANOVA_repeat_measures
    # keep <- rowSums(is.na(pcl$meta[,covars.unl[ok_covars]])) == 0
    # if (!all(keep)) {
    #     # Warn in case too many get dropped.. this is usually <1%, but
    #     # silent upstream errors can cause this to jump
    #     warning(sprintf("Dropping %d samples due to NAs in metadata..", sum(!keep)))
    #     pcl <- pcl %>% pcl.filter.s(keep=keep) %>% pcl.filter.s(any(x>0))
    # }
    keep <- !is.na(pcl$meta$subject) & (pcl.apply.s(pcl, any(x>0)))
    if (!all(keep)) {
        warning(sprintf("Dropping %d samples due to NA subject or zero profiles.", sum(!keep)))
        pcl <- pcl %>% pcl.filter.s(keep=keep)
    }

    # Get the distance matrix
    D <- vegdist(pcl$x, method, binary=method=="jaccard")

    # Find relevant covariates
    ok_covars[ok_covars] <- apply(pcl$meta[,covars.unl[ok_covars]], 2, function(x)!all(x==na.omit(x)[1], na.rm=T))
    covars <- covars.unl[ok_covars]

    # Make a per-person covariate table
    pcl$meta$subject <- factor(pcl$meta$subject)
    subj_meta <- pcl$meta[order(pcl$meta$subject),,drop=F]
    subj_meta <- subj_meta[!duplicated(subj_meta$subject), covars[covars %in% perindividual_covariates], drop=F]
    rownames(subj_meta) <- subj_meta$subject

    # Do the overall PERMANOVA
    overall.ad <- PERMANOVA_repeat_measures(
        D, permutations=Nperms, na.rm=T,
        permute_within=pcl$meta[,covars, drop=F])

    # Modify the main PERMANOVA to match the expected covariate structure
    ad <- overall.ad
    ad$aov.tab <- ad$aov.tab[seq_len(length(covariates)+1),]
    ad$aov.tab[,] <- NA
    rownames(ad$aov.tab) <- c(names(covariates), "All")
    ad$aov.tab["All",] <- colSums(overall.ad$aov.tab[seq_len(nrow(overall.ad$aov.tab)-2),])
    ad$aov.tab["All","Pr(>F)"] <- overall.ad$aov.tab$`Pr(>F)`[nrow(overall.ad$aov.tab)-1]

    ad$N <- pcl$ns

    # Do per-covariate tests
    for (var in names(covariates)) {
        varsub <- covariates[[var]]
        if (any(ok_covars[varsub])) {
            if (any(!ok_covars[varsub])) {
                varsub <- varsub[ok_covars[varsub]]
            }

            # Subset the data if necessary for this covariate
            mt <- pcl$meta
            subj_meta_sub <- subj_meta
            Dsub <- D
            if (var %in% names(subsetted_variables)) {
                ss <- subsetted_variables[[var]]
                for (i in 1:(length(ss)/3)) {
                    j <- (i-1) * 3 + 1
                    keep <- mt[,ss[j]] %in% ss[j+2]
                    keep <- if (ss[j+1]) {keep} else {!keep}
                    mt <- mt[keep,, drop=F]
                    Dsub <- as.dist(as.matrix(Dsub)[keep, keep])
                }
                mt$subject <- factor(mt$subject)
                subj_meta_sub <- subj_meta[match(levels(mt$subject), subj_meta$subject),]
            }

            if ("subject" %in% varsub) {
                # Subject is among the covariates, so the permutations must be free
                sep.ad <- PERMANOVA_repeat_measures(
                    D, permutations=Nperms, na.rm=T,
                    permute_within=mt[,varsub, drop=F])
            } else {
                # Per-covariate test with subject-aware permutations
                sep.ad <- PERMANOVA_repeat_measures(
                    D=Dsub, permutations=Nperms, na.rm=T,
                    permute_within=mt[,varsub[!(varsub %in% perindividual_covariates)], drop=F],
                    blocks=mt$subject,
                    block_data=subj_meta_sub[,varsub[varsub %in% perindividual_covariates], drop=F])
            }

            # Store results
            ad$aov.tab[var,] <- colSums(sep.ad$aov.tab[seq_len(nrow(sep.ad$aov.tab)-2),])
            ad$aov.tab[var,"Pr(>F)"] <- sep.ad$aov.tab$`Pr(>F)`[nrow(sep.ad$aov.tab)-1]
        }
    }

    return (ad)
}




# PERMANOVA tests

source("./common/load_bugs.r")
source("./common/load_kos.r")
source("./common/load_metabolites.r")
source("./common/load_proteins.r")
source("./common/load_biopsies.r")

source("./common/load_diet.r")
hmp2_diet_num_nona.pcl <- hmp2_diet_num.pcl %>%
    pcl.filter.s(!any(is.na(x)))

hmp2_omnibus_tests <- function(week_offset=0) {
    sink(sprintf("./overview/omnibus_tests_dt%d.txt", week_offset))
    cat("\n  Omnibus test results")
    cat("\n  --------------------")
    cat("\n\nThis file is generated by omnibus_tests.r\n\n")

    R2 <- matrix(0, nrow=length(covariates) + 1, ncol=0)
    rownames(R2) <- c(names(covariates), "Overall")
    dofs <- R2
    P <- R2
    ds_dofs <- c()
    add_tests <- function(ad, name) {
        print(ad)
        R2 <<- cbind(R2, ad$aov.tab$R2)
        colnames(R2)[ncol(R2)] <<- name
        dofs <<- cbind(dofs, ad$aov.tab$Df)
        colnames(dofs)[ncol(dofs)] <<- name
        P <<- cbind(P, ad$aov.tab$`Pr(>F)`)
        colnames(P)[ncol(P)] <<- name
        ds_dofs <<- c(ds_dofs, ad$N)
    }

    cat("\n ------------------------------------------------------------------------------")
    cat("\n    Bugs: Species-level Bray-Curtis\n")
    add_tests(bc_omnibus_tests(bugs.pcl %>% pcl.only(rank="s") %>% pcl.nicenames, week_offset), "Taxonomy")

    cat("\n ------------------------------------------------------------------------------")
    cat("\n    KOs (DNA): Unstratified KO DNA Bray-Curtis\n")
    add_tests(bc_omnibus_tests(ko.dna.unstrat.pcl %>% pcl.normalize, week_offset), "KOs-DNA")

    cat("\n ------------------------------------------------------------------------------")
    cat("\n    KOs (RNA): Unstratified KO RNA Bray-Curtis\n")
    add_tests(bc_omnibus_tests(ko.rna.unstrat.pcl %>% pcl.normalize, week_offset), "KOs-RNA")

    # cat("\n ------------------------------------------------------------------------------")
    # cat("\n    ECs (DNA): Unstratified EC DNA Bray-Curtis\n")
    # add_tests(bc_omnibus_tests(ec.dna.unstrat.pcl, week_offset), "ECs-DNA")
    #
    # cat("\n ------------------------------------------------------------------------------")
    # cat("\n    Pathways (DNA): Unstratified pathway DNA Bray-Curtis\n")
    # add_tests(bc_omnibus_tests(pwy.dna.unstrat.pcl, week_offset), "Pathways-DNA")
    #
    # cat("\n ------------------------------------------------------------------------------")
    # cat("\n    ECs (RNA): Unstratified EC RNA Bray-Curtis\n")
    # add_tests(bc_omnibus_tests(ec.rna.unstrat.pcl, week_offset), "ECs-RNA")
    #
    # cat("\n ------------------------------------------------------------------------------")
    # cat("\n    Pathways (RNA): Unstratified pathway RNA Bray-Curtis\n")
    # add_tests(bc_omnibus_tests(pwy.rna.unstrat.pcl, week_offset), "Pathways-RNA")

    #cat("\n ------------------------------------------------------------------------------")
    #cat("\n    Bug Expression (RNA/DNA): Euclidean distance of logratio of summed stratified ECs\n")
    #print(bc_omnibus_tests(bug.rnadna.pcl, week_offset), "EC RNA/DNA")

    cat("\n ------------------------------------------------------------------------------")
    cat("\n    Proteins (KO rollup): Bray-Curtis\n")
    add_tests(bc_omnibus_tests(proteins.kos.pcl %>% pcl.normalize, week_offset), "Proteins")

    cat("\n ------------------------------------------------------------------------------")
    cat("\n    Metabolites: Bray-Curtis, with normalization within Method\n")
    add_tests(bc_omnibus_tests(metabolites.pcl.nrm, week_offset), "Metabolites")

    #cat("\n ------------------------------------------------------------------------------")
    #cat("\n    Viruses: Jaccard\n")
    #add_tests(bc_omnibus_tests(metabolites.pcl.nrm, week_offset, method="jaccard"), "Viruses")

    if (week_offset == 0) {
        cat("\n ------------------------------------------------------------------------------")
        cat("\n    Biopsy(16S): Bray-Curtis\n")
        add_tests(bc_omnibus_tests(biopsy_16s.pcl %>% pcl.normalize, week_offset, act_lenience=4), "Biopsy (16S)")

        cat("\n ------------------------------------------------------------------------------")
        cat("\n    Biopsy(HTX): Bray-Curtis\n")
        add_tests(bc_omnibus_tests(biopsy_htx.counts.pcl %>% pcl.normalize, week_offset, act_lenience=4), "Biopsy (HTX)")
    }

    cat("\n ------------------------------------------------------------------------------")
    cat("\n    Diet: Bray-Curtis\n")
    add_tests(bc_omnibus_tests(hmp2_diet_num_nona.pcl, week_offset, act_lenience=4, method="manhattan"), "Diet")

    cat("\n\n")
    sink()

    res <- list(R2=R2, P=P, dofs=dofs, ds_dofs=ds_dofs)
    save(res, file=sprintf("./overview/variance_explained_results_dt%d.RData", week_offset))

    return (res)
}

vis_omnibus_test <- function(res, week_offset) {
    R2 <- res$R2
    dofs <- res$dofs
    ds_dofs <- res$ds_dofs
    P <- res$P

    # Add DOF to covariate names
    dofnames <- rownames(R2)
    min_dof <- apply(dofs, 1, min, na.rm=T)
    max_dof <- apply(dofs, 1, max, na.rm=T)
    for (i in seq_along(covariates)) {
        dofnames[i] <- if (min_dof[i] == max_dof[i]) {
            sprintf("%s [%.0f]", dofnames[i], min_dof[i])
        } else {
            sprintf("%s [%.0f-%.0f]", dofnames[i], min_dof[i], max_dof[i])
        }
    }
    rownames(R2) <- dofnames
    ds_dofnames <- colnames(R2)
    for (i in seq_len(ncol(R2))) {
        ds_dofnames[i] <- sprintf("%s [%.0f]", ds_dofnames[i], ds_dofs[i])
    }
    colnames(R2) <- ds_dofnames


    # Main heatmap
    ggp <- PERMANOVA_heatmap(R2, P, fontsize=5)

    pdf(sprintf("./overview/variance_explained_max_dt%d.pdf", week_offset), 3.4, 2.8, onefile=F)
    print(ggp)
    dev.off()
}

omn_t0 <- hmp2_omnibus_tests()
vis_omnibus_test(omn_t0, 0)

omn_t2 <- hmp2_omnibus_tests(2)
vis_omnibus_test(omn_t2, 2)
omn_tm2 <- hmp2_omnibus_tests(-2)
vis_omnibus_test(omn_tm2, -2)

omn_t4 <- hmp2_omnibus_tests(4)
vis_omnibus_test(omn_t4, 4)
omn_tm4 <- hmp2_omnibus_tests(-4)
vis_omnibus_test(omn_tm4, -4)
