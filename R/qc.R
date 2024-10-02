#!/usr/bin/env Rscript

###############################################################################
# omicsArt

# Copyright (c) 2020 Ali Rahnavard

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
###############################################################################

# load in the required libraries, report an error if they are not installed
for (lib in c("tibble", "logging", "data.table", "dplyr", "ggrepel", "ggplot2", "scales", "gridExtra", "future", "plyr", "cowplot")) {
  if (!suppressPackageStartupMessages(require(lib, character.only = TRUE))) stop(paste("Please install the R package: ", lib))
}
###############################################################
# If running on the command line, load other qc modules #
###############################################################

# this evaluates to true if script is being called directly as an executable
if (identical(environment(), globalenv()) &&
  !length(grep("^source\\(", sys.calls()))) {
  # source all R in qc package, relative to this folder
  script_options <- commandArgs(trailingOnly = FALSE)
  script_path <- sub("--file=", "", script_options[grep("--file=", script_options)])
  script_dir <- dirname(script_path)
  script_name <- basename(script_path)

  for (R_file in dir(script_dir, pattern = "*.R"))
  {
    if (!(R_file == script_name)) source(file.path(script_dir, R_file))
  }
}



massPattern <- function(input, r2_threshold = 0.90, q_threshold = 0.0001, coef_threshold = 0.1,
                        min_prevalence = 2, output_path = "./", outputname = "results.txt") {
  dir.create(file.path(output_path), showWarnings = FALSE)
  # print (c(r2_threshold , q_threshold, coef_threshold, min_prevalence))
  loaded_data <- load_data(input = input, type = "all", name = "Compound")
  data <- loaded_data$data
  data <- numeric_dataframe(data)
  sample_metadata <- loaded_data$sample_metadata
  feature_metadata <- loaded_data$feature_metadata


  data <- as.data.frame(t(data))
  metadata <- as.data.frame(t(sample_metadata))
  metadata <- as.data.frame(metadata["injection order"], drop = F)
  colnames(metadata) <- "injection_order"

  metadata <- numeric_dataframe(metadata)
  # metadata$injection_order <- as.numeric(metadata$injection_order) * 1.0

  data <- filter_by_prevelance(data, prev = min_prevalence, axis = 2)

  # fit linear model
  fit_data <- fit.data(features = data2, metadata = metadata, model = "SLM")
  r2_distribution <- as.data.frame(fit_data$results$r2, drop = F)
  r2_distribution$rank <- NA
  colnames(r2_distribution)[1] <- "R2"
  r2_distribution <- r2_distribution[order(r2_distribution$R2), ]
  r2_distribution$rank <- seq.int(nrow(r2_distribution))
  colnames(r2_distribution) <- c("y", "x")
  pdf(paste(output_path, "/r2_distribution.pdf", sep = ""), width = 2.75, height = 2.25, onefile = TRUE)
  r2_distribution_plot <- ggplot2::ggplot(
    data = r2_distribution,
    ggplot2::aes(x = as.numeric(as.character(x)), y = as.numeric(as.character(y)))
  ) +
    ggplot2::geom_point(color = "black", alpha = .1, shape = 21, size = 2, stroke = 0.1) +
    ggplot2::scale_x_continuous(limits = c(min(r2_distribution["x"]), max(r2_distribution["x"]))) +
    ggplot2::scale_y_continuous(limits = c(min(r2_distribution["y"]), max(r2_distribution["y"]))) +
    ggplot2::guides(alpha = "none") +
    ggplot2::labs("") +
    ggplot2::labs(title = "Distribution of coefficient of determination", x = "rank of R2", y = "R2") +
    theme_nature()

  stdout <- capture.output(print(r2_distribution_plot), type = "message")
  if (length(stdout) > 0) logging::logdebug(stdout)
  dev.off()
  print("r2_distribution_plot done!")
  # filter for high correlated features with injection order
  to_be_removed <- fit_data$results[which(fit_data$results$r2 >= r2_threshold &
    fit_data$results$qval <= q_threshold &
    abs(fit_data$results$coef) >= coef_threshold), "feature"]

  filtered_data <- data[!rownames(data) %in% to_be_removed, ]
  removed_data <- data[rownames(data) %in% to_be_removed, ]

  df <- readxl::read_excel(input, sheet = 1, col_names = FALSE)
  df2 <- tibble::add_column(as.data.frame(df), to_be_removed = NA, .after = 1)

  df2$to_be_removed <- ifelse(df2[, 1] %in% to_be_removed, "Yes", "No")
  write.table(df2, paste(output_path, "/", outputname, sep = ""),
    sep = "\t", eol = "\n", quote = F, col.names = F, row.names = F
  )
  write.table(fit_data$results, paste(output_path, "/stats.txt", sep = ""),
    sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T
  )
  # plot one example
  pool <- ifelse(startsWith(rownames(data2), "PREF"), 2, 1) #' Refrence pool', 'Sample')
  if (length(to_be_removed) > 0) {
    pdf(paste(output_path, "/features_associated_w_injection_order_r2.pdf", sep = ""), width = 2.75, height = 2.25, onefile = TRUE)
    for (i in 1:length(to_be_removed)) {
      input_df <- as.data.frame(cbind(metadata$injection_order, data2[to_be_removed[i]], pool))
      x_label <- "Injection order"
      y_label <- "Intensity"
      title <- to_be_removed[i]
      colnames(input_df) <- c("x", "y", "pool")
      temp_plot <- ggplot2::ggplot(
        data = input_df,
        ggplot2::aes(as.numeric(as.character(x)), as.numeric(as.character(y)))
      ) +
        ggplot2::geom_point(fill = pool, color = "black", alpha = .75, shape = 21, size = 2, stroke = 0.1) +
        ggplot2::scale_x_continuous(limits = c(min(input_df["x"]), max(input_df["x"]))) +
        ggplot2::scale_fill_manual(values = c("darkolivegreen4", "orange")) +
        ggplot2::scale_y_continuous(limits = c(min(input_df["y"]), max(input_df["y"]))) +
        ggplot2::stat_smooth(method = "glm", color = "blue", na.rm = T) +
        ggplot2::guides(alpha = "none") +
        ggplot2::labs("") +
        ggplot2::labs(title = title, x = x_label, y = y_label) +
        theme_nature()
      stdout <- capture.output(print(temp_plot), type = "message")
      if (length(stdout) > 0) logging::logdebug(stdout)
    }
    dev.off()
    print("diagnostics plots done!")
  }
}
