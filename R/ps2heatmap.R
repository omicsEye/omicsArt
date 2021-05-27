
# suitable for PERMENOVA test
# by Jason Lloyd-price and Ali Rahnavard

#' Combined p-value and test statistics score "heatmap"
#'
#' @param R2 A matrix of R2 values
#' @param P A matrix of P-values
#' @param fontsize The desired font size of the % R2 values in the heatmap.
#' @param FDR If \code{T}, P-values are first BH-adjusted, and significance
#' stars are shown as *** < 0.001, ** < 0.01, * < 0.05.
#' @return ggplot object.


#' @export
ps2heatmap <- function(R2, P, fontsize=5, FDR=T, alpha=NA, beta=NA) {
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
  df$NAtext[is.na(df$R2)] <- "" #"N/A"
  df$stars <- cut(df$P, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), c("***", "**", "*", ""))
  df$stars[is.na(df$R2)] <- ""

  # Try to make a reasonable color scheme that has contrast where needed
  colors <- colorRampPalette(brewer.pal(n = 9, name = "Blues"))(100)
  R2Q <- quantile(R2, c(0.25, 0.75), na.rm=T)
  if (is.na(alpha) || is.na(beta)) {
    labhat <- optim(par=c(0, 0), method="Nelder-Mead",
                    fn=function(lab) sum((pbeta(c(0.25, 0.75), exp(lab[1]), exp(lab[2])) - R2Q)^2))
    abhat <- exp(labhat$par)
    colorvalues <- pbeta(seq(0, 1, length=length(colors)), abhat[1], abhat[2])

    # Label colors flip to white when the color is too dark
    df$lblcolor <- ifelse(qbeta(df$R2, abhat[1], abhat[2]) < 0.8, "black", "white")
    #cat(sprintf("Best-fit alpha = %g, beta = %g\n", abhat[1], abhat[2]))
  } else {
    colorvalues <- pbeta(seq(0, 1, length=length(colors)), alpha, beta)

    # Label colors flip to white when the color is too dark
    df$lblcolor <- ifelse(qbeta(df$R2, alpha, beta) < 0.8, "black", "white")
  }


  ggp <- ggplot(data=df, aes(x=Dataset, y=Feature)) +
    geom_tile(aes(fill=R2)) +
    geom_text(aes(label=varExpPct, color=lblcolor), size=fontsize/(14/5), nudge_y=-0.15) +
    geom_text(aes(label=NAtext), color="grey", size=fontsize/(14/5), nudge_y=-0.15) +
    geom_text(aes(label=stars, color=lblcolor), fontface="bold", size=1.25*fontsize/(14/5), nudge_y=0.12) +
    scale_fill_gradientn(colors=colors, values=colorvalues, limits=c(0, 1), na.value="white") +
    scale_color_manual(values=c(black="black", white="white")) +
    scale_x_discrete(expand=c(0,0)) + xlab(NULL) +
    scale_y_discrete(expand=c(0,0), position = "right", limits = rev(levels(df$Feature))) + ylab(NULL) +
    guides(color="none",
           fill=guide_colourbar(title=NULL, barheight=unit(40,"mm"), label.position = "left")) +
    theme_omicsEye() +
    theme(axis.text.x = element_text(angle=-17, hjust=0),
          panel.border=element_rect(fill=NA),
          legend.position = "left", axis.ticks.y = element_blank())

  return (ggp)
}
