nature_theme <- function(x_axis_labels, y_label) {
  # set default text format based on categorical and length
  angle = NULL
  hjust = NULL
  size = 8
  if (max(nchar(x_axis_labels), na.rm=TRUE) > 5) {
    angle = 45
    hjust = 1
    size = 6
  }
  axis_title_size = 10
  if (nchar(y_label) > 15) {
    axis_title_size = 8
  }
  if (nchar(y_label) > 25) {
    axis_title_size = 6
  }
  return ( ggplot2::theme_bw() + ggplot2::theme(
    axis.text.x = ggplot2::element_text(size = size, vjust = 1, hjust = hjust, angle = angle),
    axis.text.y = ggplot2::element_text(size = 8, hjust = 1),
    axis.title = ggplot2::element_text(size = axis_title_size),
    plot.title = ggplot2::element_text(size = 7, face = 'bold'),
    #legend.title = ggplot2::element_text(size = 6, face = 'bold'),
    #legend.text = ggplot2::element_text(size = 6),
    axis.line = ggplot2::element_line(colour = 'black', size = .25),
    axis.line.x = ggplot2::element_line(colour = 'black', size = .25),
    axis.line.y = ggplot2::element_line(colour = 'black', size = .25),
    panel.border = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank())
  )
}
alpha_boxplot <- function(stats_table, meta = "group", value = "value", pvalue= NA, xlabel = meta, ylabel= NA  ){
  if(meta %in% colnames(stats_table))
    colnames(stats_table) <- gsub(meta, "group", colnames(stats_table))
  if(value %in% colnames(stats_table)){
    colnames(stats_table) <- gsub(value, "value", colnames(stats_table))
    #if (!is.na(orderby) && orderby == value)
    #  orderby= "value"
  }
  alpha_plot <- ggplot2::ggplot(
    data = subset(stats_table, !is.na(group)), ggplot2::aes(factor(group), y=value)) +
    ggplot2::geom_boxplot(
      ggplot2::aes(fill = group),
      outlier.alpha = 0.0,
      na.rm = TRUE,
      alpha = .5,
      show.legend = FALSE

    ) +
    ggplot2::geom_point(
      ggplot2::aes(fill = group),
      alpha = 0.75 ,
      size = 1,
      shape = 21,
      stroke = 0.15,
      color = 'black',
      position = ggplot2::position_jitterdodge()) +
    ylab(ylabel) +  xlab(xlabel) +
    ggplot2::annotate(
      geom = "text",
      x = Inf,
      y = 0,
      hjust = 1,
      vjust = 1,
      label = sprintf(
        "Kruskal-Wallis p-value: %.4f",
        pvalue
      ) ,
      color = "black",
      size = 2,
      fontface = "italic"
    )+
    theme_omicsEye() + nature_theme(stats_table["value"],ylabel)+theme(legend.position = "none")
  return(alpha_plot)
}

alpha_scatterplot <- function(stats_table, meta = "group", value = "value", pvalue= NA, xlabel = meta, ylabel= NA  ){
  if(meta %in% colnames(stats_table))
    colnames(stats_table) <- gsub(meta, "group", colnames(stats_table))
  if(value %in% colnames(stats_table)){
    colnames(stats_table) <- gsub(value, "value", colnames(stats_table))
    #if (!is.na(orderby) && orderby == value)
    #  orderby= "value"
  }
  #pvalue <- NA
  #res <- cor.test(sapply(stats_table["group"], as.numeric), sapply(stats_table["value"], as.numeric), method = "spearman")
  #pvalue <- res$p.value
  print (pvalue)
  alpha_plot <- ggplot2::ggplot(
    data = subset(stats_table, !is.na(group)), ggplot2::aes(group, y=value)) +
    ggplot2::stat_smooth(
      method = "glm",
      size = 0.5,
      color = 'grey',
      na.rm = TRUE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(fill = value),
      alpha = 0.75 ,
      size = 1,
      shape = 21,
      stroke = 0.15,
      color = 'black',
      position = ggplot2::position_jitterdodge()) +
    ggplot2::guides(alpha = 'none') +
    ylab(ylabel) +  xlab(xlabel) +
    ggplot2::annotate(
      geom = "text",
      x = Inf,
      y = 0,
      hjust = 1,
      vjust = 1,
      label = sprintf(
        "Spearman p-value: %.4f",
        pvalue
      ) ,
      color = "black",
      size = 2,
      fontface = "italic"
    )+
    theme_omicsEye() +theme(legend.position = "none")
  return(alpha_plot)
}

alpha_diversity_all <- function(data, metadata){
  alpha_diversity_data <- metadata
  alpha_diversity_data$alpha <- diversity(t(data), "shannon")

  alpha_diversity_test <-
    setNames(data.frame(matrix(
      nrow = length(unique(colnames(metadata))),  ncol = 2,
    )), c("Metadata", "P.value"))
  rownames(alpha_diversity_test) <- unique(colnames(metadata))
  alpha_diversity_plots <- vector('list', length(unique(colnames(metadata))))
  names(alpha_diversity_plots) <- unique(colnames(metadata))
  i = 1
  pdf(
    paste('analysis', '/alpaha_diversity.pdf', sep = ''),
    width = 2.4,
    height = 2.25,
    onefile = TRUE
  )
  for (meta in unique(colnames(metadata))){
    temp_pval <- NA

    tryCatch({
      print(meta)
      #print(tem_kruskal.test$p.value)
      #alpha_diversity_data[meta,"P.value" ] <-tem_kruskal.test$p.value
      if (is.numeric(alpha_diversity_data[colnames(data), meta]) &
          length(unique(alpha_diversity_data[colnames(data), meta])) > 2) {
        res <- cor.test(sapply(alpha_diversity_data[colnames(data), meta], as.numeric),
                        sapply(alpha_diversity_data[colnames(data), "alpha"], as.numeric), method = "spearman")
        temp_pval <- res$p.value
        print(temp_pval)
        alpha_diversity_plots[[meta]] <- alpha_scatterplot(stats_table=alpha_diversity_data ,
                                                           meta = meta, value = "alpha",
                                                           pvalue=temp_pval, xlabel = meta, ylabel= 'Shannon index')
      }else{
        tem_kruskal.test <- kruskal.test(alpha_diversity_data[colnames(data), "alpha"], alpha_diversity_data[colnames(data), meta] )
        temp_pval <- tem_kruskal.test$p.value
        alpha_diversity_plots[[meta]] <- alpha_boxplot(stats_table=alpha_diversity_data ,
                                                       meta = meta, value = "alpha",
                                                       pvalue=temp_pval, xlabel = meta, ylabel= 'Shannon index')
      }
      alpha_diversity_test[meta, "P.value"] <- temp_pval
      alpha_diversity_test[meta, "Metadata"] <- meta
      i <- i + 1
      stdout <-
        capture.output(print(alpha_diversity_plots[[meta]]), type = "message")
    }, error = function(e) {
      print(meta)
      print(paste('error:', e))
    })
  }
}
