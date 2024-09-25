# This script contains functions to manipulate metbolite profiling data
#' @export
numeric_dataframe <- function(input) {
  input[, c(1:ncol(input))] <-
    sapply(sapply(input[, c(1:ncol(input))], as.character), as.numeric)
  return (input)

}
# convert2mxpformat <- function(data, sample_metadata=NA, feature_metdata=NA){
#
#   if(!is.na(feature_metdata)){
#     combined <- as.matrix(cbind(t(feature_metadata), data))
#   }
#   if(!is.na(sample_metadata)){
#     empty_corner <- matrix("", nrow = dim(sample_metadata)[1], ncol = dim(feature_metadata)[2])
#     ext_sample_metadata <- as.matrix(cbind(empty_corner,sample_metadata))
#     ext_sample_metadata[,dim(feature_metadata)[2]] <- rownames(sample_metadata)
#     combined2 <- as.matrix(rbind(ext_sample_metadata,combined))
#   }
#
# }

find_indx <- function(df, word = 'Metabolite')
{
  colmn <- which(df == word) %/% nrow(df) + 1
  row <- which(df == word) %% nrow(df)
  return(c(row, colmn))
}
#' @export
load_data <-
  function(input,
           type = 'known',
           sheet = 1,
           ID = 'Metabolite') {
    if (is.character(input)) {
      #df <- read.xls (input, sheet = sheet, header = FALSE)
      df <-
        readxl::read_excel(input, sheet = sheet, col_names = FALSE) #much faster
    } else {
      df <- input
    }
    df <- as.data.frame(df)
    df[df == "n/a"] <- NA
    df[df == "NA"] <- NA
    df[df == ""] <- NA
    # remove columns with all NAs
    df <- Filter(function(x)
      !all(is.na(x)), df)

    # find the occutance of Metabolite word in the file
    ord <- find_indx(df, word = 'Metabolite')
    row <- ord[1]
    col <- ord[2]

    Compound_ord <- find_indx(df, word = 'Compound_ID')
    HMDB_ord <- find_indx(df, word = 'HMDB_ID')
    MZ_ord <- find_indx(df, word = 'MZ')
    RT_ord <- find_indx(df, word = 'RT')
    # seprate sample metadata
    col_names <- NA
    if (row > 1) {
      sample_metadata <- df[1:(row - 1), col:dim(df)[2]]
      meta_rownames <- sample_metadata[[1]]
      sample_metadata <-
        sample_metadata[1:dim(sample_metadata)[1], 2:dim(sample_metadata)[2]]
      col_names <- sapply(df[row, (col + 1):dim(df)[2]], as.factor)
      colnames(sample_metadata) <- col_names
      rownames(sample_metadata) <- meta_rownames
    } else{
      col_names <- sapply(df[row, (col + 1):dim(df)[2]], as.factor)
      sample_metadata <- NA
    }
    # seprate feature metadata
    f_col_names <- NA
    if (col > 1) {
      feature_metadata <- df[row:dim(df)[1], 1:col]
      f_meta_rownames <- feature_metadata[[Compound_ord[2]]]
      feature_metadata <-
        feature_metadata[2:dim(feature_metadata)[1], 1:dim(feature_metadata)[2]]
      f_col_names <- sapply(df[row, 1:col], as.factor)
      colnames(feature_metadata) <- f_col_names
      rownames(feature_metadata) <-
        f_meta_rownames[-1] # ro remove "Compound" from the names
    } else{
      feature_metadata <- NA
    }
    df <- as.data.frame(df)


    #seprate known data or generate a ID for unknows in the case of using all
    if (type == 'known') {
      if (ID == 'Metabolite') {
        df <- df[!is.na(df[, col]), ]
      } else if (ID == 'HMDB_ID') {
        df <- df[!is.na(df[, HMDB_ord[2]]), ]
        df <- df[!grepl('\\*', df[, HMDB_ord[2]]),]
        df <- df[df[, HMDB_ord[2]] != "", ]
      }
    } else if (type == 'all') {
      if (ID == 'Compound_ID') {
        name_col <- df[, Compound_ord[2]] #df[,col-1]
      } else{
        if (ID == 'Metabolite') {
          name_col = col
        } else if (ID == 'HMDB_ID') {
          name_col = HMDB_ord[2]
        }
        df[, col] <- ifelse(is.na(df[, name_col]),
                            paste(
                              df[, Compound_ord[2]],
                              'MZ',
                              round(as.numeric(df[, MZ_ord[2]]), digits = 4),
                              'RT',
                              round(as.numeric(df[, RT_ord[2]]), digits = 2),
                              sep = '_'
                            ),
                            as.character(df[, name_col]))
      }
    }

    #remove redundant ion
    df[is.na(df)] <- "NA"
    if (!is.na(HMDB_ord[2])) {
      df <- df[which(df[, HMDB_ord[2]] != "redundant ion"),]
      df <- df[df[, HMDB_ord[2]] != "Internal Standard", ]
    }
    df <- df[df[, col] != "Metabolite", ]
    # remove feature start with NH4_
    df <- df[!startsWith(sapply(df[, col], as.character), "NH4_"),]
    df[df == "NA"] <- NA
    df <- df[!is.na(df[, col]), ]

    if (ID == 'HMDB_ID') {
      row_names <-
        gsub('\\*', '', as.matrix(df[row:dim(df)[1], HMDB_ord[2]]))
    } else if (ID == 'Metabolite') {
      row_names <- unlist(sapply(df[row:dim(df)[1], col], as.factor))
    } else if (ID == 'Compound_ID') {
      row_names <-
        unlist(sapply(df[row:dim(df)[1], Compound_ord[2]], as.factor))
    }
    data <- df[row:dim(df)[1], col:dim(df)[2]]
    data <- data[,-1]
    # check if row_names has unique values to be used for row names
    dup_names <- row_names[duplicated(row_names)]
    if (length(dup_names) > 0) {
      print(sprintf("There are duplicated values: %s", dup_names))
    }
    rownames(data) <- row_names
    colnames(data) <- col_names
    data <- as.data.frame(data)
    result <- list()
    # put data in the right ordination rows(observation or samples) columns features
    result$sample_metadata <- as.data.frame(t(sample_metadata))
    result$feature_metadata <- feature_metadata
    data <- as.data.frame(t(data))

    # make sure data frame is numeric
    result$data <- numeric_dataframe(data)
    return (result)
  }
stats_tabletest_2d <- function(df1, df2) {
  results <-
    setNames(data.frame(matrix(
      ncol = dim(df1)[2], nrow = dim(df2)[2]
    )), colnames(df2))
  rownames(results) <- colnames(df1)
  for (i in 1:dim(df1)[2]) {
    for (j in 1:dim(df2)[2]) {
      if (i == j)
        results[i, j] <- t.test(df1[, i], df2[, j])
    }
  }
  return(results)
}
cor_2d <-
  function(df1,
           df2,
           method = 'spearman',
           use = "na.or.complete") {
    results <-
      setNames(data.frame(matrix(
        ncol = dim(df1)[2], nrow = dim(df2)[2]
      )), colnames(df2))
    rownames(results) <- colnames(df1)
    for (i in 1:dim(df1)[2]) {
      for (j in 1:dim(df2)[2]) {
        results[i, j] <-
          cor(df1[, i], df2[, j], method = method , use = use)
      }
    }
    return(results)
  }

cor_within_and_between_2df <- function(df1,
                                       df2,
                                       output_path,
                                       metadata = NA,
                                       fillby = NA,
                                       method = 'spearman',
                                       use = "pairwise.complete.obs",
                                       all = F,
                                       plot_cor_threshold = .65,
                                       plot_min_points_threshold = 10) {
  commom_rows <- intersect(rownames(df1), rownames(df2))
  df1 <- df1[commom_rows,]
  df2 <- df2[commom_rows,]

  n_row = dim(df1)[2]
  if (all) {
    n_row = dim(df1)[2] * (dim(df1)[2] - 1)  + dim(df2)[2] * (dim(df2)[2] -
                                                                1) + dim(df1)[2]
  }
  stats_table <-
    setNames(
      data.frame(matrix(ncol = 7, nrow = n_row)),
      c(
        "Feature1",
        "Feature2",
        "Correlation",
        "P.Value",
        "Between",
        "Dataset",
        "Characterization"
      )
    )
  N <- 0


  # within dataset 1
  if (all) {
    for (i in 1:dim(df1)[2]) {
      for (j in 1:dim(df1)[2]) {
        if (i == j)
          next
        N <- N + 1
        stats_table [N, ] <- c(
          colnames(df1)[i],
          colnames(df1)[j],
          cor(df1[, i], df1[, j], method = method , use =
                use),
          NA,
          "Within",
          '1'
        )
      }
    }

    # within dataset 2
    for (i in 1:dim(df2)[2]) {
      for (j in 1:dim(df2)[2]) {
        if (i == j)
          next
        N <- N + 1
        stats_table [N, ] <- c(
          colnames(df2)[i],
          colnames(df2)[j],
          cor(df2[, i], df2[, j], method = method , use =
                use),
          NA,
          "Within",
          '2'
        )
      }
    }
  }
  # between two datasets
  for (i in 1:dim(df1)[2]) {
    bothNArows <- which(is.na(df1[, i]) | is.na(df2[, i]))
    n_not_NA <- dim(df1)[1] - length(bothNArows)
    if (n_not_NA < plot_min_points_threshold)
      next
    N <- N + 1
    Characterization <-
      ifelse(startsWith(colnames(df1)[i], 'QI'), 'Unknown', 'Known')
    cor_test <-
      cor.test(df1[, i], df2[, i], method = method , use = use)
    #cor(df1[,i], df2[,i], method = method , use=use)
    stats_table [N, ] <- c(
      colnames(df1)[i],
      colnames(df2)[i],
      cor_test$estimate,
      cor_test$p.value,
      "Between",
      3,
      Characterization
    )
  }
  stats_table <- stats_table[!is.na(stats_table$Correlation),]
  stats_table <-
    stats_table[order(-abs(as.numeric(stats_table$Correlation))),]
  # write the stats table
  write.table(
    stats_table,
    paste(output_path, '/', 'correlation_between_features.txt', sep = ''),
    sep = "\t",
    eol = "\n",
    quote = F,
    col.names = NA,
    row.names = T
  )

  # distribution of corrlations
  plot_file <-
    paste(output_path, "/correlation_plots", fillby, ".pdf", sep = "")
  pdf(plot_file,
      width = 2.65,
      height = 2.5,
      onefile = TRUE)
  steps <- 0.05
  ggp <-
    ggplot(stats_table,
           aes(
             x = as.numeric(Correlation),
             #fill = Characterization,
             na.rm = TRUE
           )) +
    #geom_density(alpha = 0.35, size = 0.1) +
    geom_density(
      aes(y = ..count..*steps, fill = Characterization, colour = Characterization),
      alpha = 0.25,
      size = 0.3
    ) +
    #facet_wrap(~Characterization)+
    geom_histogram(
      aes(y = ..count..),
      breaks = seq(-1, 1.0, by = steps),
      size = 0.1,
      col = "black",
      fill = "gold",
      alpha = .25
    ) +
    guides(fill = guide_legend(title = NULL)) + scale_colour_discrete(guide = 'none') +
    scale_x_continuous(limits = c(min(as.numeric(
      stats_table$Correlation
    )), max(as.numeric(
      stats_table$Correlation
    )))) +
    xlab("Spearman Correlation") + ylab("Count") +
    theme_omicsEye() +
    theme(legend.justification = c(1, 1),
          legend.position = c(1, 1))
  stdout <- capture.output(print(ggp), type = "message")
  if (length(stdout) > 0)
    logging::logdebug(stdout)

  stats_table_plots <-
    stats_table[abs(as.numeric(stats_table$Correlation)) > .7,]
  scatter_plots <-
    vector(mode = "list", length = dim(stats_table_plots)[1])
  names(scatter_plots) <- stats_table_plots$Feature1
  logging::loginfo("Plotting data ...")


  for (i in 1:dim(stats_table_plots)[1]) {
    #i = 1
    x_label <- as.character(stats_table_plots[i, 'Feature1'])
    y_label <- as.character(stats_table_plots[i, 'Feature2'])
    pval <- as.numeric(stats_table_plots[i, 'P.Value'])
    coef_val <- as.numeric(stats_table_plots[i, 'Correlation'])
    input_df <- cbind(df1[x_label], df2[y_label])
    n_not_NA <- 0
    colnames(input_df) <- c("x", "y")
    fillby_one_color <- NA
    if (!is.na(metadata) & !is.na(fillby)){
      input_df[, fillby] <- metadata[rownames(input_df), fillby]
    }else {
      input_df[, fillby] <- NA
      fillby_one_color <- 'darkolivegreen4'
    }
    bothNArows <- which(is.na(input_df$x) | is.na(input_df$y))
    n_not_NA <- dim(input_df)[1] - length(bothNArows)
    print(n_not_NA)
    # if Metadata is continuous generate a scatter plot
    # Continuous is defined as numerical with more than 2 values (to exclude binary data)
    if (plot_min_points_threshold > 0)
    {
      if (n_not_NA < plot_min_points_threshold)
        # for positive correlation
        next
    }else if (n_not_NA > plot_min_points_threshold) {
      # for negative correlations
      next
    }
    temp_plot <- NULL
    logging::loginfo("Creating scatter plot for continuous data, %s vs %s",
                     x_label,
                     y_label)
    temp_plot <- ggplot2::ggplot(data = input_df,
                                 ggplot2::aes(as.numeric(as.character(x)), as.numeric(as.character(y)))) +
      ggplot2::geom_point(
        aes(fill = get(fillby)),#'darkolivegreen4',
        #fill = fillby_one_color,#'darkolivegreen4',
        color = 'black',
        alpha = .5,
        shape = 21,
        size = 1,
        stroke = 0.15
      ) +
      ggplot2::scale_x_continuous(limits = c(min(input_df['x']), max(input_df['x']))) +
      ggplot2::scale_y_continuous(limits = c(min(input_df['y']), max(input_df['y']))) +
      ggplot2::stat_smooth(
        method = "glm",
        color = 'blue',
        size = .5,
        na.rm = T
      ) +
      ggplot2::guides(alpha = 'none') + ggplot2::labs("") +
      ggplot2::xlab(x_label) +  ggplot2::ylab(y_label) + theme_omicsEye() +
      guides(fill = guide_legend(title = fillby))
      ggplot2::annotate(
        geom = "text",
        x = Inf,
        y = Inf,
        hjust = 1,
        vjust = 1,
        label = sprintf(
          "P-Value: %s\nCoefficient: %s",
          formatC(pval, format = "e", digits = 3),
          formatC(coef_val, format = "e", digits = 2)
        ) ,
        color = "black",
        size = 2,
        fontface = "italic"
      )
    if (!is.na(fillby_one_color)){
      temp_plot <- temp_plot + ggplot2::geom_point(fill = fillby_one_color)
    }
    scatter_plots[[x_label]] <- temp_plot
    stdout <- capture.output(print(temp_plot), type = "message")
    if (length(stdout) > 0)
      logging::logdebug(stdout)
  }

  dev.off()
  results <- list()
  results$stats_table  <- stats_table
  results$dist_plot <- ggp
  results$scatter_plots <- scatter_plots
  return(results)
}

# transformation function for reverse log1p axis
revlog_trans <- function(base = exp(1)) {
  trans <- function(x)
    - log1p(x)
  inv <- function(x)
    expm1(-x)
  scales::trans_new("revlog1p", trans, inv, domain = c(0, Inf))
}

#' @export
volcano_plot <-
  function(stats_table,
           threshold = 0.05,
           method = 'nominal',
           pvalue_col = "pval",
           fdr_col = "qval",
           orderby = 'coef',
           y_label = "P-Value (-log10)",
           x_label = 'Coefficient') {
    colnames(stats_table) <- gsub(orderby, "coef", colnames(stats_table))
    colnames(stats_table) <- gsub(pvalue_col, "P.Value", colnames(stats_table))
    colnames(stats_table) <- gsub(fdr_col, "fdr", colnames(stats_table))
    if (method == 'nominal')
      stats_table$fdr <- stats_table$P.Value
    stats_table$feature[stats_table$fdr  >= threshold] <- NA
    p <-
      ggplot(stats_table, aes(
        x = coef,
        y = -log10(P.Value),
        label = feature
      )) +
      scale_fill_gradient(low = "lightgray", high = "navy") +
      scale_color_gradient(low = "lightgray", high = "navy") +
      #scale_y_continuous(trans = revlog_trans(), expand = c(0.005, 0.005)) +
      #scale_y_log10()+
      #expand_limits(y = c(0, 1)) +
      stat_density_2d(aes(fill = ..level..),
                      geom = "polygon",
                      show.legend = FALSE) +
      geom_point(
        data = subset(stats_table, fdr <= threshold & coef > 0.0),
        fill = "green",
        color = 'black',
        alpha = .5,
        shape = 21,
        size = 1,
        stroke = 0.05
      ) +
      geom_point(
        data = subset(stats_table, fdr <= threshold & coef < 0.0),
        fill = 'red',
        color = 'black',
        alpha = .5,
        shape = 21,
        size = 1,
        stroke = 0.05
      ) +
      geom_point(
        data = subset(stats_table, fdr > threshold),
        fill =  'gray',
        color = "black",
        alpha = 0.2,
        shape = 21,
        size = 1,
        stroke = 0.05
      ) +
      geom_vline(xintercept = 0, size = 0.1) +
      geom_hline(yintercept = 0, size = 0.1) +
      geom_hline(
        yintercept = -log10(threshold),
        linetype = "dashed",
        size = 0.1,
        color = "red"
      ) +
      geom_vline(
        xintercept = c(-1, 1),
        linetype = "dashed",
        size = 0.1,
        color = "red"
      ) +
      theme_linedraw() +
      theme(panel.grid = element_blank()) +
      xlab(x_label) +
      ylab(y_label) +
      annotate(
        "text",
        x = min(stats_table$coef) + .1,
        y = min(-log2(stats_table$P.Value)) + .75,
        label = "Significant up",
        size = 3,
        color = "black",
        hjust = 0
      ) +
      annotate(
        "point",
        x = min(stats_table$coef),
        y = min(-log2(stats_table$P.Value)) + .75,
        color = "green"
      ) +

      annotate(
        "text",
        x = min(stats_table$coef) + .1,
        y = min(-log2(stats_table$P.Value)) + .5,
        label = "Significant down",
        size = 3,
        color = "black",
        hjust = 0
      ) +
      annotate(
        "point",
        x = min(stats_table$coef),
        y = min(-log2(stats_table$P.Value)) + .5,
        color = "red"
      ) +

      annotate(
        "text",
        x = min(stats_table$coef),
        y = min(-log2(stats_table$P.Value)) + 1.0,
        label = paste(method, " threshold: ", threshold, sep = ""),
        size = 3,
        color = "black",
        hjust = 0
      ) +
      geom_text_repel(
        size = 2,
        force = 1,
        fontface =  "italic")#,
        #box.padding = unit(.01, "lines"),
        #point.padding = unit(0.01, "lines")
      #)
    return (p)
  }
mean_test <-
  function(data,
           metadata,
           meta ,
           case = NA,
           control = NA,
           method = 'wilcox.test',
           paired = F) {
    stats_table <-
      setNames(
        data.frame(matrix(
          ncol = 8, nrow = dim(data)[2]
        )),
        c(
          "logFC",
          "Carbon length",
          "Double bond",
          "statistic",
          "P.Value",
          "adj.P.Val",
          "B",
          "fdr"
        )
      )
    rownames(stats_table) <- colnames(data)

    commom_rows <- intersect(rownames(data), rownames(metadata))
    # response variables must be numeric
    data <- numeric_dataframe(data[commom_rows,])
    metadata <- metadata[commom_rows, ]
    if (!is.na(case)) {
      case <- data[which(metadata[meta,] == case)]
    }
    if (!is.na(control)) {
      control <- data[which(metadata[meta,] == control)]
    }
    if (paired) {
      commom_rows2 <- intersect(rownames(case), rownames(control))
      case <- case[commom_rows2,]
      control <- control[commom_rows2,]
    }
    for (i in 1:dim(stats_table)[1]) {
      if (all(is.na(case[, i])) || all(is.na(control[, i]))) {
        print(i)
        next
      }
      stats_table[i, 'logFC'] <-
        log2(mean(case[, i], na.rm = TRUE) / mean(control[, i], na.rm = TRUE))
      tryCatch({
        if (method == 't.test') {
          stats_table[i, 'P.Value'] <-
            t.test(case[, i], control[, i], paired = paired)$p.value
          stats_table[i, 'statistic'] <-
            t.test(case[, i], control[, i], paired = paired)$statistic
        } else{
          stats_table[i, 'P.Value'] <-
            wilcox.test(case[, i], control[, i], paired = paired)$p.value
          stats_table[i, 'statistic'] <-
            wilcox.test(case[, i], control[, i], paired = paired)$statistic
        }
      }, error = function(e) {
        stats_table[i, 'P.Value'] <- NA
        stats_table[i, 'statistic'] <- NA
      })
      stats_table[i, "Carbon length"] <-
        as.numeric(unlist(strsplit(unlist(
          strsplit(rownames(stats_table)[i], ":")
        )[1], "C"))[2])
      stats_table[i, "Double bond"] <-
        as.numeric(unlist(strsplit(unlist(
          strsplit(rownames(stats_table)[i], ":")
        )[2], " "))[1])
      stats_table[i, "Family"] <-
        gsub("^.* ", "", rownames(stats_table)[i])
    }

    stats_table$fdr <-
      p.adjust(
        as.vector(stats_table$P.Value),
        method = 'BH',
        n = length(stats_table$P.Value)
      )
    stats_table$feature <- rownames(stats_table)
    return (stats_table)
  }

#' @export
stats_2groups <-
  function(case,
           control,
           test_type = 'wilcox.test',
           paired = F) {
    stats_table <-
      setNames(
        data.frame(matrix(
          ncol = 5, nrow = dim(case)[2]
        )),
        c(
          "logFC",
          "statistic",
          "P.Value",
          "adj.P.Val",
          "fdr"
        )
      )
    rownames(stats_table) <- colnames(case)
    #i <- 100
    for (i in 1:dim(stats_table)[1]) {
      if (all(is.na(case[, i])) || all(is.na(control[, i]))) {
        #print(i)
        next
      }
      stats_table[i, 'logFC'] <-
        log2(mean(case[, i], na.rm = TRUE)) - log2(mean(control[, i], na.rm = TRUE))
      tryCatch({
        if (test_type == 't.test') {
          temp_test_result <- NA
          temp_test_result <- t.test(case[, i], control[, i], paired = paired)
          stats_table[i, 'P.Value'] <-temp_test_result$p.value
          stats_table[i, 'statistic'] <-temp_test_result$statistic
        } else{
            temp_test_result <- NA
            temp_test_result <-wilcox.test(case[, i], control[, i], paired = paired)
            stats_table[i, 'P.Value'] <- temp_test_result$p.value
          stats_table[i, 'statistic'] <-temp_test_result$statistic
        }
      }, error = function(e) {
        stats_table[i, 'P.Value'] <- NA
        stats_table[i, 'statistic'] <- NA
      })
    }

    stats_table$fdr <-
      p.adjust(
        as.vector(stats_table$P.Value),
        method = 'BH',
        n = length(stats_table$P.Value)
      )
    stats_table$feature <- rownames(stats_table)
    return (stats_table)
  }

#' @export
diff_bar_plot <-
  function(stats_table,
           threshold = 0.05,
           method = 'nominal',
           coef = "coef",
           pvalue_col = "P.Value",
           fdr_col = "fdr",
           orderby = "coef",
           y_label = "Y label",
           x_label = 'X label',
           feature_col = "feature") {
    if(coef %in% colnames(stats_table) ){
      colnames(stats_table) <- gsub(coef, "coef", colnames(stats_table))
      if (!is.na(orderby) && orderby == coef)
        orderby= "coef"
    }
    if(pvalue_col %in% colnames(stats_table)){
      colnames(stats_table) <- gsub(pvalue_col, "P.Value", colnames(stats_table))
      pvalue_col = "P.Value"
      if (!is.na(orderby) && orderby == pvalue_col)
        orderby= "P.Value"
    }
    if(fdr_col %in% colnames(stats_table)){
      colnames(stats_table) <- gsub(fdr_col, "fdr", colnames(stats_table))
      fdr_col = "fdr"
      if (!is.na(orderby) && orderby == fdr_col)
        orderby= "fdr"
    }
    if(feature_col %in% colnames(stats_table)){
      colnames(stats_table) <- gsub(feature_col, "feature", colnames(stats_table))
      if (!is.na(orderby) && orderby == feature_col)
        orderby= "feature"

    }
    if (method == 'nominal') {
      stats_table <-
        stats_table[which(stats_table[,pvalue_col] <= threshold),]
    }else if (method == 'fdr') {
      stats_table <-
        stats_table[which(stats_table[,fdr_col] <= threshold),]
    }
    #stats_table2 <- stats_table
    #stats_table <- stats_table2
    if (!is.na(orderby)){
      stats_table <- stats_table[order(stats_table[orderby]),]
      order_sig <- rownames(stats_table)
      stats_table <- within(stats_table,
                            feature <- factor(feature,
                                              levels=order_sig))
      #stats_table <- transform( stats_table,
      #                          feature = ordered(feature, levels = names( sort(-table(feature)))))
      #stats_table <- within(stats_table,
      #                      feature <- factor(feature,
      #                                        levels=names(sort(table(feature),
      #                                                          decreasing=TRUE))))
      #stats_table2 <- stats_table[order(stats_table[orderby]),]
    }
    p <- ggplot(stats_table) +
      geom_bar(
        aes(
          x= feature, #reorder(rownames(stats_table), stats_table$logFC),
          y = coef,
          fill = coef
        ),
        show.legend =  F,
        stat = "identity",
        position = "identity"
      ) +
      #scale_alpha_manual(values=-log10(stats_table$qval),guide=F)+
      geom_text(aes(x= feature,
                    y = coef,
                    label=ifelse(stats_table[fdr_col] < threshold, "*", "") ##sprintf("%s",formatC(qval, format = "e", digits = 3)
                    ),
                size = 2
                )+
      xlab(y_label) + ylab(x_label) +
      #guides(colour = guide_legend(show = FALSE) )+
      #guides(fill = guide_legend(title = "")) + theme(legend.justification =
       #                                                 c(0, 0),
       #                                               legend.position = c(.3, .2)) +
      #theme(legend.position="top", legend.direction="vertical")+
      theme(plot.title = element_text(face = "bold", hjust = .9)) +
      #labs(title = "")+
      scale_fill_gradient2(
        low = 'red',
        mid = 'snow3',
        high = 'darkgreen',
        space = 'Lab'
      ) +
      theme_omicsEye() +

      #theme(text = element_text(),axis.text.x = element_text(angle = 0, hjust = 1), axis.text.y = element_text(angle = 0, hjust = 1))+
      #geom_text(aes(label=number_of_genome_ref), position=position_dodge(width=0.9), vjust=-0.25)+
      #scale_x_discrete()+
      #scale_y_continuous(expand = c(0)) +
      coord_flip()
    return (p)
  }

check_anotation <- function(study, ref_db = NA) {
  if (is.na(ref_db)) {
    ref_db <-
      read.table(
        "/Users/rah/Documents/Metabolomics/example_data/HMDB/mastermapping.txt",
        header = TRUE,
        sep = "\t",
        fill = FALSE,
        comment.char = "" ,
        check.names = FALSE
      )
  }
  unkcharterized_hmdb <-
    rownames(data)[!rownames(data) %in% ref_db$HMDBID]
  unkcharterized_hmdb <-
    as.data.frame(unkcharterized_hmdb, drop = F)
  return (unkcharterized_hmdb)
}

#' @export
variance.test <-
  function(data, metadata, meta , method = 'levene') {
    # find common samples (rows) between data and metdata if the are not the same
    if (!(meta %in% colnames(data)) &&
        !(meta %in% colnames(metadata))) {
      print('meta does not exist in your data')
      return (NA)
    }
    commom_rows <- intersect(rownames(data), rownames(metadata))
    # response variables must be numeric
    data <- numeric_dataframe(data[commom_rows,])
    metadata2 <- as.data.frame(metadata[commom_rows, meta], drop = F)

    # combine data nad metadata for test function
    test_data <- cbind(metadata2, data)

    # group factor defining groups.
    test_data$meta <- as.factor(test_data$meta)

    # set up results table
    stats_table <-
      setNames(data.frame(matrix(
        ncol = 3, nrow = dim(data)[2]
      )),
      c("Df", "F value", "Pr(>F)"))

    #
    rownames(stats_table) <- colnames(data)
    for (i in 1:dim(stats_table)[1]) {
      if (all(is.na(data[, i]))) {
        print(i)
        next
      }
      tryCatch({
        if (method == 'levene') {
          # feature to be tested
          feature <- rownames(stats_table)[i]
          test_data[feature] <- lapply(test_data[feature], as.numeric)
          # an alternative test
          # #bf.test(cepgfr5t ~ race1c, data = data)
          test_results <-
            with(test_data, leveneTest(get(feature), meta))
          stats_table[i, 'Df'] <- test_results$`Df`[1]
          stats_table[i, 'F value'] <- test_results$`F value`[1]
          stats_table[i, 'Pr(>F)'] <- test_results$`Pr(>F)`[1]
        } else{
          return(NA)
        }
      }, error = function(e) {
        stats_table[i, 'Df'] <- NA
        stats_table[i, 'F value'] <- NA
        stats_table[i, 'Pr(>F)'] <- NA
      })
    }
    desnity_pvalues_plot <-
      ggplot(data = stats_table, aes(get('Pr(>F)'))) +
      geom_histogram(
        aes(y = ..density..),
        breaks = seq(0.0, 1.0, by = 0.05),
        col = "black",
        fill = "gold",
        alpha = .5
      ) +
      #geom_density(col=2) +
      labs(
        title = sprintf('p-values of variantion in %s groups', meta),
        x = sprintf('P-values using %s test', method),
        y = 'Density'
      )
    result <- list()
    result$stats_table <- stats_table
    result$desnity_pvalues_plot <- desnity_pvalues_plot
    return (result)
  }

regrese_out <- function(features,
                        meta,
                        random = NA,
                        fixed = NA) {
  random_effect <- factor(meta[, random])
  fixed_effect <- meta[, fixed]

  residuals_df <- NULL
  for (i in 1:dim(features)[2]) {
    #i<- 1
    lmer_results <- lmer(features[, i] ~  fixed_effect +
                           (1 |
                              random_effect), na.action = na.exclude)
    temp_residual <- resid(lmer_results)
    residuals_df <- rbind(residuals_df, temp_residual)
  }
  row.names(residuals_df) <- colnames(features)
  colnames(residuals_df) <- rownames(features)
  return (as.data.frame(t(residuals_df)))
}

#' @export
merge_cvs <- function(input_dir_path, outputfile) {
  # Download from http://smpdb.ca/downloads
  # extract the zip file to smpdb_metabolites
  #input_dir_path = "~/Downloads/smpdb_metabolites"
  #outputfile = '/Users/rah/Documents/Metabolomics/example_data/HMDB/merged_ref_db.csv'
  library(dplyr)
  library(readr)
  df <- list.files(path = input_dir_path, full.names = TRUE) %>%
    lapply(read_csv) %>%
    bind_rows
  write_csv(df, outputfile)
}
filter_by_mean <- function(data, top = .95, axis = 1) {
  if (axis == 2) {
    data <- as.data.frame(t(data))
  }
  step1.dat <-
    data[rowMeans(data, na.rm = TRUE) >= quantile(rowMeans(data, na.rm = TRUE), 1.0 -
                                                    top), ]
  if (axis == 2)
    return (as.data.frame(t(step1.dat)))
  else
    return (as.data.frame(step1.dat))
}
filter_by_mean <- function(data, top = .95, axis = 1) {
  if (axis == 2) {
    data <- as.data.frame(t(data))
  }
  step1.dat <-
    data[rowMeans(data, na.rm = TRUE) >= quantile(rowMeans(data, na.rm = TRUE), 1.0 -
                                                    top), ]
  if (axis == 2)
    return (as.data.frame(t(step1.dat)))
  else
    return (as.data.frame(step1.dat))
}
filter_by_prevelance <-
  function(data,
           prev = 2,
           type = "number",
           axis = 1) {
    if (axis == 2) {
      data <- as.data.frame(t(data))
    }
    data2 <- data
    data2[data2 > 0] <- 1
    #data[is.na(data)] <- 0
    if (type == "number")
      step1.dat <- data[rowSums(data2, na.rm = TRUE) >= prev, ]
    else if (type == "perc")
      step1.dat <-
      data[rowMeans(data, na.rm = TRUE) >= quantile(rowMeans(data, na.rm = TRUE), 1.0 -
                                                      top), ]
    else
      # if (type == "percentile")
      step1.dat <-
      data[rowMeans(data, na.rm = TRUE) >= quantile(rowMeans(data, na.rm = TRUE), 1.0 -
                                                      top), ]
    if (axis == 2)
      return (as.data.frame(t(step1.dat)))
    else
      return (as.data.frame(step1.dat))
  }

#' @export
filter_by_variance <- function(data, top = .95, axis = 1) {
  if (axis == 2) {
    data <- as.data.frame(t(data))
  }
  step2.dat <- apply(data, 1, var, na.rm = TRUE)
  step2.dat <-
    data[step2.dat >= quantile(step2.dat, 1.0 - top, na.rm = TRUE), ]
  if (axis == 2)
    return (as.data.frame(t(step2.dat)))
  else
    return (as.data.frame(step2.dat))
}

#' @export
ordplots <-
  function(data,
           metadata,
           output,
           outputname = NA,
           method = 'pcoa',
           index = 'bray/curtis') {
    dir.create(file.path(output), showWarnings = FALSE)
    commom_rows <- intersect(rownames(data), rownames(metadata))
    metadata <- metadata[commom_rows,]
    data <- data[commom_rows,]
    fakepcl <- list(
      meta = as.data.frame(metadata),
      x = as.matrix(data),
      ns = dim(data)[1],
      nf = dim(data)[2]
    )
    if (method == 'tsne') {
      ord <- pcl.tsne(fakepcl)
    } else if (method == 'pcoa') {
      ord <- pcl.pcoa(fakepcl, k = 3, index = 'bray/curtis')#, index = 'euclidean'
    } else if (method == 'pca') {
      ord <- pcl.pcoa(fakepcl, k = 3, index = 'euclidean')#, index = 'euclidean'
    }else{
      print("The ordination is not implemeted in omicsArt!")
    }
    if (is.na(outputname)) {
      outputname <- method
    }
    ord_plots <- vector(mode = "list", length = dim(metadata)[2])
    names(ord_plots) <- colnames(metadata)
    ord_plot <- pcl.ordplot(fakepcl, ord,  pcos = c(1, 2))
    #logging::loginfo("Plotting data using %s ordination", method)
    pdf(
      paste(output, '/', outputname, '.pdf', sep = ''),
      width = 4.5,
      height = 3.25,
      onefile = TRUE
    )
    for (temp_metadat in colnames(metadata)) {
      #print(temp_metadat)
      #print(is.numeric(fakepcl$meta[1, temp_metadat]))
      #print(length(fakepcl$meta[[temp_metadat]]))
      #print(fakepcl$meta[1, temp_metadat])
      # get not NA values
      not_na <- fakepcl$meta[temp_metadat]
      not_na <- not_na[!is.na(not_na)]
      l_unique <- length(unique(not_na))
      if (check.numeric(not_na[1]) & l_unique > 2)
        #tryCatch({
        fakepcl$meta[temp_metadat] <-
        sapply(fakepcl$meta[temp_metadat], as.numeric)
      #}, error = function(e){
      #  next#fakepcl$meta[temp_metadat] <- lapply(fakepcl$meta[temp_metadat], as.factor)
      #})
      else if (l_unique > 100 || l_unique < 2) {
        #print(temp_metadat)
        #print(l_unique)
        #print(not_na[1])
        #print(check.numeric(not_na[1]))
        next
      }

      ord_plot <-
        pcl.ordplot(
          fakepcl,
          ord,
          colour_title = temp_metadat,
          colour = temp_metadat,
          pcos = c(1, 2)
        )
      ord_plots[[temp_metadat]] <- ord_plot
      stdout <- capture.output(print(ord_plot), type = "message")
      if (length(stdout) > 0)
        logging::logdebug(stdout)

      #ggsave(filename=paste(output,'/', temp_metadat,'.pdf', sep=''), plot=ord_plot,
      #       width = 5, height = 4, units = "in", dpi = 300)
    }
    dev.off()
    return(ord_plots)
  }

#' @export
check.numeric <- function(N) {
  return(
    grepl(
      "[-]?^[[:digit:]]+[.]?^[[:digit:]]*|[-]?^[[:digit:]]+[L]?|[-]?^[[:digit:]]+[.]?^[[:digit:]]*[eE][0-9]+",
      N
    )
  )
}

#' @export
convert_maaslin_output2matrix <-
  function(df,
           simlarity_threshold = 0.0005,
           q_threshold = 1.0,
           first_n = NA,
           cell_value = 'pval') {
    metadata <- df$metadata
    feature <- df$feature
    value <- NA
    if (!is.na(first_n) & first_n > 0 & first_n < dim(df)[1]) {
      if (cell_value == 'coef') {
        df <- df[order(-abs(df[cell_value])) ,]
      } else{
        df <- df[order(df[cell_value]),]
      }
      df <- df[1:first_n, ]
    }
    df$coef[abs(df$qval) > q_threshold |
              abs(df$coef) < simlarity_threshold] <- 0
    value <- -log(df$qval) * sign(df$coef)
    value <- pmax(-20, pmin(20, value))

    n <- length(unique(feature))
    m <- length(unique(metadata))
    if (n < 2) {
      print(
        'There is no enough features in associations to make a heatmap plot of associations.
        Please look at association in text output file!'
      )
      return (NA)
    }
    if (m < 2) {
      print(
        'There is no enough metadata in associations to make a heatmap plot of associations.
        Please look at association in text output file!'
      )
      return (NA)
    }
    a = matrix(0, nrow = n, ncol = m)
    a <- as.data.frame(a)
    rownames(a) <- unique(feature)
    colnames(a) <- unique(metadata)
    for (i in 1:dim(df)[1]) {
      if (abs(a[as.character(feature[i]), as.character(metadata[i])]) > abs(value[i]))
        next
      a[as.character(feature[i]), as.character(metadata[i])] <-
        value[i]
    }
    #if (data_type == "pathway"){
    #  colnames(a) <-  sapply(colnames(a), shorter_pwy_names )
    #}
    #else{
    #  colnames(a) <-  sapply(colnames(a), pcl.sub )
    #}
    #rownames(b) <- rownames(a)#all_metadata
    #colnames(b) <- colnames(a)
    b <- a
    #for (i in 1:n)
    #  for (j in 1:m){
    #    b[rownames(a)[i], colnames(a)[j] ] <- a[rownames(a)[i], colnames(a)[j] ]
    #  }
    #b <-b[needed_metadata,]
    #rownames(b) <- needed_metadata_names
    b <- b[, colSums(b != 0, na.rm = TRUE) > 0, drop = F]
    b <- t(b)
    b <- b[, colSums(b != 0, na.rm = TRUE) > 0, drop = F]
    #b <- t(b)
    if (length(colnames(b)) > 2) {
      hc <- hclust(dist(t(b)), method = "single")
      b <- b[, hc$order, drop = F]
    }

    return (as.data.frame(b))
  }

#' @export
combine_maaslin_heatmap <- function() {
  #library(pheatmap)
  output_file = "pheatmap.pdf"
  cell_value = 'qval'
  title = ""
  if (cell_value == "pval") {
    value <- -log(df$pval) * sign(df$coef)
    value <- pmax(-20, pmin(20, value))
    if (is.null(title))
      title <- "Significant associations (-log(pval)*sign(coeff))"
  } else if (cell_value == "qval") {
    value <- -log(df$qval) * sign(df$coef)
    value <- pmax(-20, pmin(20, value))
    if (is.null(title))
      title <- "Significant associations (-log(qval)*sign(coeff))"
  } else if (cell_value == "coef") {
    value <- df$coef
    if (is.null(title))
      title <- "Significant associations (coeff)"
  }
  cell_value = "Q.value"
  data_label = 'Data'
  metadata_label = 'Metadata'
  border_color = "grey93"
  color = colorRampPalette(c("blue", "whitesmoke", "red"))(51) #whitesmoke
  # read MaAsLin output
  gaps <- c(0, 0)
  mat <- NA
  site_colours <- set_eric_coloures()
  bodysite_color_order = list(
    "Body site" = c(
      "Buccal mucosa" = site_colours$Buccal_mucosa,
      "Supragingival plaque" = site_colours$Supragingival_plaque,
      "Tongue dorsum" = site_colours$Tongue_dorsum,
      "Stool" = site_colours$Stool,
      "Anterior nares" = site_colours$Anterior_nares,
      "Posterior fornix" = site_colours$Posterior_fornix
    )
  )
  #print(bodysite_color_order)
  df_1 <-
    read.table(
      '/Users/rah/Documents/Metabolomics/Projects/MESA/C8-pos_known/linear_model_output_exam1/significant_results.tsv',
      header = TRUE,
      sep = "\t",
      fill = TRUE,
      comment.char = "" ,
      check.names = FALSE
    )
  mat <- convert_maaslin_output2matrix(df_1)
  #mat <- mat[, colnames(mat)%in%names(bugs_to_show_assoc[["Anterior nares"]])]
  print ("C8_pos")
  #colnames(mat) <- names(bugs_to_show_assoc[["Anterior nares"]])
  gaps[1] <- length(colnames(mat))
  method_order = data.frame(ID = factor(rep(c("C8-pos"), each = length(colnames(
    mat
  )))))

  df_2 <-
    read.table(
      '/Users/rah/Documents/Metabolomics/Projects/MESA/HILIC-pos_known/linear_model_output_exam1/significant_results.tsv',
      header = TRUE,
      sep = "\t",
      fill = TRUE,
      comment.char = "" ,
      check.names = FALSE
    )
  mat2 <- convert_maaslin_output2matrix(df_2)
  if (!is.na(mat2) && dim(mat2)[1] > 0) {
    #mat2 <- mat2[, colnames(mat2)%in%names(bugs_to_show_assoc[["Buccal mucosa"]])]
    print ("HILIC-pos")
    #colnames(mat2) <- bugs_to_show_assoc[["Buccal mucosa"]]
    #mat <- cbind(mat, mat2)
    mat <- combine2df(mat, mat2)
    gaps[2] <- length(colnames(mat))
    method_order = rbind(method_order, data.frame(ID = factor(rep(
      c("HILIC-pos"), each = length(colnames(mat2))
    ))))
  }

  df_3 <-
    read.table(
      '/Users/rah/Documents/Metabolomics/Projects/MESA/HILIC-neg_known/linear_model_output_exam1/significant_results.tsv',
      header = TRUE,
      sep = "\t",
      fill = TRUE,
      comment.char = "" ,
      check.names = FALSE
    )
  mat3 <- convert_maaslin_output2matrix(df_3)
  if (!is.na(mat3) && dim(mat3)[1] > 0) {
    #mat3 <- mat3[, colnames(mat3)%in%names(bugs_to_show_assoc[["Supragingival plaque"]])]
    print("HILIC-neg")
    #print(length(colnames(mat3)))
    #colnames(mat3) <- bugs_to_show_assoc[["Supragingival plaque"]]
    #mat <- cbind(mat, mat3)
    mat <- combine2df(mat, mat3)
    gaps[3] <- length(colnames(mat))
    method_order = rbind(method_order, data.frame(ID = factor(rep(
      c("HILIC-neg"), each = length(colnames(mat3))
    ))))
  }
  #df_4 <- read.table(  '/Users/rah/Documents/Metabolomics/Projects/MESA/C8-pos_known/linear_model_output_exam1/significant_results.tsv',
  #                    header = TRUE, sep = "\t", fill = TRUE, comment.char = "" , check.names = FALSE)
  #mat4 <- convert_maaslin_output2matrix(df_4)
  #mat4 <- mat4[, colnames(mat4)%in%names(bugs_to_show_assoc[["C18-neg"]])]

  #print(length(colnames(mat4)), names(bugs_to_show_assoc[["Tongue dorsum"]]))
  #colnames(mat4) <- bugs_to_show_assoc[["C18-neg"]]
  #mat <- cbind(mat, mat4)
  #print("C18-neg")
  #gaps[4] <- length(colnames(mat))
  #method_order = rbind(method_order, data.frame(ID = factor(rep(c("Tongue dorsum"), each=length(colnames(mat4))))))
  colnames(method_order)[1] <- "Profiling method"
  #rownames(bodysite_color_order)[1] <- "Body site"
  #rownames(bodysite_order) <- colnames(mat)
  labels_col <- colnames(mat)
  colnames(mat) <- rownames(method_order)
  #bodysite_order <- as.matrix(bodysite_order)
  #rownames(bodysite_order) <- colnames(mat)
  #bodysite_order <- as.data.frame(bodysite_order)
  #print (rownames(bodysite_order))
  #pdf(file=output_file, height = length(colnames(b))/5+8, width = length(rownames(b))/5+5)
  #pdf(file=output_file, height = 25, width = 10)
  mylegend = TRUE
  my_show_rownames = TRUE
  #if (!grepl('nares', output_file)){
  #  mylegend = FALSE
  #  my_show_rownames = FALSE
  #}
  mat[is.na(mat)] <- 0
  plot_result <-
    pheatmap(
      mat,
      cellwidth = 5,
      cellheight = 5,
      # changed to 3
      main = title,
      annotation_col = method_order["Profiling method"],
      #annotation_colors = bodysite_color_order["Profiling method"],
      annotation_legend = T,
      fontsize = 6,
      kmeans_k = NA,
      #border=TRUE,
      labels_col = labels_col,
      show_rownames = T,
      show_colnames = T,
      scale = "none",
      #clustering_method = "complete",
      cluster_rows = FALSE,
      cluster_cols = F,
      clustering_distance_rows = "euclidean",
      clustering_distance_cols = "euclidean",
      legend = TRUE,
      border_color = border_color,
      color = color,
      treeheight_row = 0,
      treeheight_col = 0,
      gaps_col = gaps,
      display_numbers = matrix(ifelse(mat > 0.0, "+", ifelse(mat < 0.0, "-", "")),  nrow(mat))
    )
  dev.off()
  ggsave(
    filename = 'associations_overview.pdf',
    plot = plot_result$gtable,
    width = 12,
    limitsize = FALSE,
    height = 3,
    units = "in",
    dpi = 300
  )
  return (plot_result)
}

#' @export
combine2df <- function(df1, df2) {
  all_rows <- dplyr::union(rownames(df1), rownames(df2))
  df1 <- df1[all_rows,]
  df2 <- df2[all_rows,]
  result <- cbind(df1, df2)
  rownames(result) <- all_rows
  return(result)
}
#' @export
combine_QI_TF <- function(QI_file, TF_file, output_name){
  # read the Progenesis file
  QI_data <- read.table(
    QI_file,
    header = F,
    row.names = NULL,
    sep = ",",
    fill = FALSE,
    comment.char = "" ,
    check.names = FALSE
  )

  # name the columns
  colnames(QI_data) <- lapply(QI_data[3,], as.character)

  # find the sart of raw abundance data
  data_indx <- find_indx(QI_data, word = 'Raw abundance')

  # remove the first two columns
  QI_data <- QI_data[4:nrow(QI_data), ]

  # use raw abundance
  QI_data <- QI_data[, c(1,3,5,data_indx[2]:ncol(QI_data))]

  # re-order columns and remove columns with name ""
  indx_last_sample <- grep("Accepted Compound ID", colnames(QI_data))
  cols <- colnames(QI_data)
  cols <- cols[1:indx_last_sample]
  samples_cols <-
    cols[!(cols %in% c(
      "Compound",
      "m/z",
      "Retention time (min)",
      "Accepted Compound ID",
      ""
    ))]
  samples_cols <- samples_cols[order(samples_cols)]
  cols_in_order <-
    c(c("Compound", "m/z", "Retention time (min)", "Accepted Compound ID"),
      samples_cols)
  QI_data <- QI_data[, cols_in_order]


  #####read first file for TR profiles #########################
  #1) on the second tab, I calculate the average retention time

  RT_profile_data <-
    readxl::read_excel(
      TF_file,
      sheet = 2,
      col_names = FALSE
    )
  colnames(RT_profile_data) <- sapply(RT_profile_data[1,], as.factor)
  RT_profile_data_rownames <- sapply(RT_profile_data[,  1], as.factor)
  RT_profile_data <- RT_profile_data[-1, -1]
  rownames(RT_profile_data) <- RT_profile_data_rownames[-1]
  RT_profile_data <- as.data.frame(RT_profile_data)

  ##### Calculate average RT
  RT_profile_data <- omicsArt:::numeric_dataframe(RT_profile_data)
  RT_profile_data[RT_profile_data <= 0.0] <- NA
  RT <- colMeans(x = RT_profile_data, na.rm = T)
  RT_profile_data <- rbind(RT = RT, RT_profile_data)

  #2) I add the RT to the first tab, then copy & transpose to a new tab

  #####read first file for TR profiles #########################
  Intensity_profile_data <-
    readxl::read_excel(
      TF_file,
      sheet = 1,
      col_names = FALSE
    )
  colnames(Intensity_profile_data) <-
    sapply(Intensity_profile_data[1,], as.factor)
  Intensity_profile_data_rownames <-
    sapply(Intensity_profile_data[,  1], as.factor)
  Intensity_profile_data <- Intensity_profile_data[-1, -1]
  rownames(Intensity_profile_data) <-
    Intensity_profile_data_rownames[-1]
  Intensity_profile_data <- as.data.frame(Intensity_profile_data)

  ##### add the average RT
  Intensity_profile_data <- omicsArt:::numeric_dataframe(Intensity_profile_data)
  Intensity_profile_data <- rbind(RT = RT, Intensity_profile_data)

  # clean data
  ## remove "Standards" mm, Bpp, FFA, miniMM
  clean_rows <- row.names(Intensity_profile_data)
  clean_rows <- clean_rows[!grepl("*miniMM*|*FFA*|*Bpp*|*mm*|*Standards*", clean_rows)]
  Intensity_profile_clean_data <- Intensity_profile_data[clean_rows,]
  TF_data <- as.data.frame(t(Intensity_profile_clean_data))
  TF_data[,c("Compound", "m/z", "Retention time (min)", "Accepted Compound ID")] <- NA
  TF_data[, "Retention time (min)"] <-  TF_data[, "RT"]

  # rename TF metabolites if they are in QI data
  row.names(TF_data) <- ifelse(row.names(TF_data) %in% row.names(QI_data),
                               paste(row.names(TF_data), "TF", sep = '_'), row.names(TF_data))

  TF_data[, "Software"] <- "TF"
  TF_data$`Accepted Compound ID` <- rownames(TF_data)
  QI_data[, "Software"] <- "QI"
  #rownames(QI_data) <- QI_data$Compound
  combined <- rbind(TF_data[, colnames(QI_data)], QI_data)

  # order the final column names
  cols_in_order <-
    c(c("Compound", 'Software', "m/z", "Retention time (min)", "Accepted Compound ID"),
      samples_cols)
  combined <- combined[, cols_in_order]

  colnames(combined) <- c(c("Compound_ID", 'Software', "MZ", "RT", "Metabolite"),
                          samples_cols)
  combined[combined == ""] <- NA
  combined <- with(combined,  combined[order(Metabolite, na.last = TRUE), ])
  hs <- createStyle(textDecoration = "BOLD", fontColour = "#FFFFFF", fontSize = 12,
                    fontName = "Arial Narrow", fgFill = "#4F80BD")
  options("openxlsx.borderColour" = "#4F80BD") ## set default border colour
  write.xlsx(combined,
             file = output_name,
             colNames = TRUE)
  #borders = "rows",
  #headerStyle = hs)
}

#' @export
lollipop_plot <-
  function(stats_table,
           threshold = 0.0,
           method = 'nominal',
           x = "x",
           y = "coef",
           pvalue_col = "P.Value",
           fdr_col = "fdr",
           orderby = "coef",
           y_label = "Y label",
           x_label = 'X label',
           feature_col = "feature",
           point_color = "darkolivegreen4",
           label = T) {
    if(x %in% colnames(stats_table) ){
      colnames(stats_table) <- gsub(x, "x", colnames(stats_table))
    }
    if(y %in% colnames(stats_table)){
      colnames(stats_table) <- gsub(y, "y", colnames(stats_table))
    }
    stats_table$x <- as.numeric(stats_table$x)
    stats_table$y <- as.numeric(stats_table$y)
    p <- ggplot(stats_table, aes(x=x, y=y, label = feature)) +
    geom_segment( aes(x=x, xend=x, y=0, yend=y), size = 0.005) +
    geom_point( fill = point_color,
                color = 'black',
                alpha = .5,
                shape = 21,
                size = .65,
                stroke = 0.05) +
    #theme_light() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.border = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    xlab(x_label) +
    ylab(y_label) +
    geom_hline(yintercept=threshold, color="red", size = 0.05) +
      scale_x_continuous(n.breaks = 10) +
      scale_y_continuous(n.breaks = 5)+
    theme_omicsEye()
    if(label)
      p <- p +
      geom_text_repel(
        size = 2,
        force = 1,
        fontface =  "italic")
    return (p)
  }
# Adapted form: https://rstudio-pubs-static.s3.amazonaws.com/455435_30729e265f7a4d049400d03a18e218db.html

#' @export
entropy <- function(target) {
  #if(all(is.na(target)))  0
  freq <- table(target)/length(target)
  # vectorize
  vec <- as.data.frame(freq)[,2]
  #drop 0 to avoid NaN resulting from log2
  vec<-vec[vec>0]
  #compute entropy
  -sum(vec * log2(vec))
}

#----------------------------------

#' @export
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

