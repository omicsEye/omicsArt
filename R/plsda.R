mplsda <- function(data, metadata, column, output) {
  if (!(column %in% colnames(metadata))) {
    stop("The specified column does not exist in the metadata")
  }
  common_rows <- intersect(rownames(data), rownames(metadata))
  data <- data[common_rows, ]
  metadata <- metadata[common_rows, ]
  data[is.na(data)] <- 0
  metadata[, column] <- factor(metadata[, column])
  if (length(levels(metadata[, column])) == 1) {
    stop("the selected column in metadata needs at least two unique labels")
  }
  dir.create(file.path(output), showWarnings = FALSE)
  label_pairs <- t(utils::combn(levels(metadata[, column]), 2))
  for (j in 1:nrow(label_pairs)) {
    label1 <- label_pairs[j, 1]
    label2 <- label_pairs[j, 2]
    filtered_metadata <- metadata[metadata[, column] %in% c(label1, label2), ]
    filtered_metadata[, column] <- factor(filtered_metadata[, column])
    filtered_data <- subset(data, rownames(data) %in% rownames(filtered_metadata))
    filtered_data <- filtered_data[, colSums(filtered_data) > 0]
    model <- ropls::opls(filtered_data, filtered_metadata[, column], predI = 1, orthoI = NA, fig.pdfC = NULL)
    dir.create(file.path(output, paste(label1, "vs", label2)), showWarnings = FALSE)
    for (plot_type in c("x-loading", "correlation", "summary")) {
      ropls::plot(model, typeVc = plot_type, fig.pdfC = file.path(output, paste(label1, "vs", label2), paste(plot_type, ".pdf", sep = "")))
    }
    utils::write.table(t(methods::slot(model, "vipVn")), file = file.path(output, paste(label1, "vs", label2), "VIP_values.csv"), sep = ",", col.names = NA)
    utils::write.table(t(methods::slot(model, "loadingMN")), file = file.path(output, paste(label1, "vs", label2), "X_loadings.csv"), sep = ",", col.names = NA)
    utils::write.table(t(methods::slot(model, "scoreMN")), file = file.path(output, paste(label1, "vs", label2), "X_scores.csv"), sep = ",", col.names = NA)
    label1_metadata <- filtered_metadata[filtered_metadata[, column] == label1, ]
    label2_metadata <- filtered_metadata[filtered_metadata[, column] == label2, ]
    label1_data <- subset(filtered_data, rownames(filtered_data) %in% rownames(label1_metadata))
    label2_data <- subset(filtered_data, rownames(filtered_data) %in% rownames(label2_metadata))
    label1_predictions <- stats::predict(model, label1_data)
    label1_correct <- sum(label1_predictions == label1)
    label1_incorrect <- sum(label1_predictions == label2)
    label2_predictions <- stats::predict(model, label2_data)
    label2_correct <- sum(label2_predictions == label2)
    label2_incorrect <- sum(label2_predictions == label1)
    predictions <- matrix(c(label1_correct, label1_incorrect, label2_incorrect, label2_correct), ncol = 2, byrow = TRUE)
    colnames(predictions) <- c(paste("thought to be", label1), paste("thought to be", label2))
    rownames(predictions) <- c(paste("truly", label1), paste("truly", label2))
    predictions <- as.table(predictions)
    utils::write.table(predictions, file = file.path(output, paste(label1, "vs", label2), "training_predictions.csv"), sep = ",")
  }
}
