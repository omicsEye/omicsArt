library(vegan)
library(permute)

#' @export
metadataCorr <- function(metadata, entropy_threshold = 0.5, p_threshold = 0.05, cluster = F) {
  metadata[metadata == "null"] <- NA
  metadata[metadata == "NULL"] <- NA
  metadata[metadata == ""] <- NA

  # remove columns with all NA
  metadata <- metadata[, colSums(is.na(metadata)) < nrow(metadata)]
  # metadata <- omicsArt::numeric_dataframe(metadata)

  # remove features with entropy <= entropy_threshold
  metadata <- metadata[, apply(metadata, 2, entropy) >= entropy_threshold]

  # create a matrix for p-values and test statistics
  M_perm <- matrix(NA, nrow = ncol(metadata), ncol = ncol(metadata))

  rownames(M_perm) <- colnames(metadata)
  colnames(M_perm) <- rownames(M_perm)
  R_perm <- P_perm <- M_perm

  for (meta1 in colnames(metadata)) {
    # meta1 is  Continuous
    if (is.numeric(metadata[1, meta1]) &
      length(unique(metadata[[meta1]])) >= 2) {
      # logging::loginfo(
      #   "Continuous data: %s",
      #   meta1)
      for (meta2 in colnames(metadata)) {
        # if (meta == meta2) next
        if (is.numeric(metadata[1, meta2]) &
          length(unique(metadata[[meta2]])) >= 2) {
          # meta2 is  Continuous
          # logging::loginfo(
          #   "Continuous data: %s",
          #   meta2)
          tryCatch(
            {
              res <- NA
              res <- cor.test(sapply(metadata[meta1], as.numeric), sapply(metadata[meta2], as.numeric), method = "spearman")
              P_perm[meta1, meta2] <- res$p.value
              R_perm[meta1, meta2] <- res$estimate
            },
            error = function(e) {
              print(meta1)
              print(meta2)
              print(paste("error:", e))
            }
          )
        } else if (length(unique(metadata[[meta2]])) >= 2) {
          # meta2 is  categorical
          # logging::loginfo(
          #   "Categorical data: %s",
          #   meta2)
          tryCatch(
            {
              res <- NA
              res <- kruskal.test(get(meta1) ~ get(meta2), data = metadata)
              P_perm[meta1, meta2] <- res$p.value
              R_perm[meta1, meta2] <- res$statistic
              print(meta1, meta2, res$statistic)
            },
            error = function(e) {
              print(meta1)
              print(meta2)
              print(paste("error:", e))
            }
          )
        }
      }
    } else if (length(unique(metadata[[meta1]])) >= 2) {
      # meta1 is  categorical
      # logging::loginfo(
      #   "Categorical data: %s",
      #   meta1)
      for (meta2 in colnames(metadata)) {
        # if (meta == meta2) next
        if (is.numeric(metadata[1, meta2]) &
          length(unique(metadata[[meta2]])) >= 2) {
          # meta2 is  Continuous
          # logging::loginfo(
          #   "Continuous data: %s",
          #   meta2)
          tryCatch(
            {
              res <- NA
              res <- kruskal.test(get(meta2) ~ get(meta1), data = metadata)
              P_perm[meta1, meta2] <- res$p.value
              R_perm[meta1, meta2] <- res$statistic
            },
            error = function(e) {
              print(meta1)
              print(meta2)
              print(paste("error:", e))
            }
          )
        } else if (length(unique(metadata[[meta2]])) >= 2) {
          # meta2 is  categorical
          # logging::loginfo(
          #   "Categorical data: %s",
          #   meta2)
          tryCatch(
            {
              res <- NA
              res <- chisq.test(metadata[[meta1]], metadata[[meta2]], correct = FALSE)
              P_perm[meta1, meta2] <- res$p.value
              R_perm[meta1, meta2] <- res$statistic
            },
            error = function(e) {
              print(meta1)
              print(meta2)
              print(paste("error:", e))
            }
          )
        }
      }
    }
  }
  # R_perm[R_perm>1]<- NA
  test_heatmap <- omicsArt::ps2heatmap(R_perm, P_perm, FDR = T)
  # P_perm <- P_perm[ , colSums(is.na(P_perm)) < 60]
  # P_perm<- P_perm[colnames(P_perm), ]

  # heatmap3::heatmap3(P_perm)
  # library(RColorBrewer)
  tryCatch(
    {
      pval_hetamap <- pheatmap::pheatmap(-log(P_perm + 0.00000001),
        # cellwidth = 5,
        # cellheight = 5,
        # changed to 3
        # main = title,
        # fontsize = 6,
        kmeans_k = NA,
        border = TRUE,
        show_rownames = TRUE,
        show_colnames = TRUE,
        scale = "none",
        cluster_rows = cluster,
        cluster_cols = cluster,
        clustering_distance_rows = "euclidean",
        clustering_distance_cols = "euclidean",
        legend = TRUE,
        border_color = "grey93",
        color = colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100), # color(range_value*100),
        # breaks = breaks,
        treeheight_row = 0,
        treeheight_col = 0,
        display_numbers = matrix(ifelse(
          P_perm < p_threshold, "*", ""
        ), nrow(P_perm))
      )
    },
    error = function(e) {
      print("Might need to have option  cluster = F")
      print(meta2)
      print(paste("error:", e))
    }
  )

  result <- list()
  result$pval_hetamap <- pval_hetamap
  result$stat_pval_hetamap <- test_heatmap
  result$P_perm <- P_perm
  result$R_perm <- R_perm
  return(result)
}
