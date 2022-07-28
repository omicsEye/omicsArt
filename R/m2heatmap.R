# suitable for mantel.rtest
#by Jason Lloyd-price and Ali Rahnavard
library(RColorBrewer)

#' @export
m2heatmap <- function(O, P, Ocil, Ociu, title, stars=F) {

  library(reshape2)
  library(viridis)
  df <- melt(O)
  colnames(df) <- c("V1", "V2", "obs")
  df$V2 <- factor(df$V2, levels=rev(levels(df$V2)))
  df$obs <- df$obs^2
  if (!missing(P)) {
    Padj <- matrix(p.adjust(P, method="fdr"), nrow=nrow(P))
    df$padj <- melt(Padj)$value
    print(Padj)
  }
  if (!missing(Ocil)) {
    df$cil <- melt(Ocil)$value^2
    df$cil[melt(Ocil)$value < 0] <- 0
    df$ciu <- melt(Ociu)$value^2
  }

  # Increase the contrast of the color scale where it matters
  alpha <- 1.9
  beta <- 0.1
  colorvalues <- pbeta(seq(0, 1, length=101), alpha, beta)
  df$lblcolor <- ifelse(qbeta(df$obs, alpha, beta) < 0.8, "black", "white")

  ggp <- ggplot(data=df, aes(x=V1, y=V2)) +
    geom_tile(aes(fill=obs)) +
    scale_fill_gradientn(colours=colorRampPalette(brewer.pal(n = 9, name = "Blues"))(100),
                         values=colorvalues,
                         na.value="white", limits=c(0, 1), name="Variance Explained")
  #ggp
  nudge_y <- 0
  if (!missing(P)) {
    ggp <- ggp + if (stars) {
      nudge_y <- -0.14
      geom_text(aes(label=ifelse(padj<=0.001, "***", ifelse(padj<=0.01, "**", ifelse(padj<=0.05, "*", ""))),
                    color=lblcolor),
                size=2, nudge_y=0.175, fontface="bold")
    } else {
      geom_text(aes(label=ifelse(is.na(padj), "", sprintf("%.2g", padj)),
                    color=lblcolor), size=1.7)
    }
  }
  ggp <- ggp + if (!missing(Ocil)) {
    geom_text(aes(label=ifelse(is.na(obs), "", sprintf("%.1f%%\n[%.1f%% - %.1f%%]", 100*obs, 100*cil, 100*ciu)),
                  color=lblcolor),
              size=1.7, nudge_y=nudge_y)
  } else {
    geom_text(aes(label=ifelse(is.na(obs), "", sprintf("%.1f%%", 100*obs)),
                  color=lblcolor),
              size=1.7, nudge_y=nudge_y)
  }
  ggp <- ggp +
    geom_text(aes(label=ifelse(V1!=V2, ifelse(is.na(obs), "", ""), "")), size=1.7, color="gray") +
    theme_omicsEye() +
    theme(axis.text.x=element_text(angle=45, hjust=1),
          legend.position = "left", axis.ticks.y = element_blank()) +
    scale_x_discrete(expand=c(0, 0)) + scale_y_discrete(expand=c(0, 0), position = "right") +
    scale_color_manual(values=c(white="white", black="black")) +
    guides(fill=guide_colourbar(title=NULL, barheight=unit(0.65,"npc"), label.position = "left"), color="none") +
    xlab(NULL) + ylab(NULL)
  # if (!missing(title)) {
  #     ggp <- ggp + ggtitle(title)
  # }

  return (ggp)
}
