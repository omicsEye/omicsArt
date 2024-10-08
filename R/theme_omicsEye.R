# library(cowplot)
# library(ggplot2)

#' @export
theme_omicsEye <- function() {
  list(
    cowplot::theme_cowplot(),
    ggplot2::theme(
      text               = ggplot2::element_text(size = 6),
      axis.text          = ggplot2::element_text(size = 5),
      axis.title.x       = ggplot2::element_text(margin = ggplot2::margin(1, 0, 0.5, 0)),
      axis.title.x.top   = ggplot2::element_text(margin = ggplot2::margin(0, 0, 2, 0)),
      axis.title.y       = ggplot2::element_text(margin = ggplot2::margin(0, 1, 0, 0.5)),
      axis.title.y.right = ggplot2::element_text(margin = ggplot2::margin(0, 0, 0, 2)),
      axis.text.x        = ggplot2::element_text(margin = ggplot2::margin(1, 0, 0, 0)),
      axis.text.x.top    = ggplot2::element_text(margin = ggplot2::margin(0, 0, 1, 0)),
      axis.text.y        = ggplot2::element_text(margin = ggplot2::margin(0, 1, 0, 0)),
      axis.text.y.right  = ggplot2::element_text(margin = ggplot2::margin(0, 0, 0, 1)),
      axis.ticks         = ggplot2::element_line(size = 0.3),
      axis.ticks.length  = ggplot2::unit(2, "pt"),
      axis.line          = ggplot2::element_line(size = 0.3),
      axis.line.x        = ggplot2::element_line(size = 0.3),
      axis.line.y        = ggplot2::element_line(size = 0.3),
      line               = ggplot2::element_line(size = 0.3),
      legend.margin      = ggplot2::margin(4, 4, 4, 4),
      legend.key.size    = ggplot2::unit(8, "pt"),
      legend.box.spacing = ggplot2::unit(4, "pt"),
      panel.spacing      = ggplot2::unit(1.5, "pt"),
      plot.title         = ggplot2::element_text(size = 8),
      plot.margin        = ggplot2::margin(1, 1, 1, 1),
      strip.background   = ggplot2::element_blank(),
      strip.text         = ggplot2::element_text(size = 6),
      strip.text.x       = ggplot2::element_text(margin = ggplot2::margin(3, 0, 3, 0)),
      strip.text.y       = ggplot2::element_text(margin = ggplot2::margin(0, 3, 0, 3))
    )
  )
}

#' @export
theme_omicsEye_presentation <- function() {
  list(
    cowplot::theme_cowplot(),
    ggplot2::theme(
      text               = ggplot2::element_text(size = 12),
      axis.text          = ggplot2::element_text(size = 7),
      axis.title.x       = ggplot2::element_text(margin = ggplot2::margin(1, 0, 0.5, 0)),
      axis.title.x.top   = ggplot2::element_text(margin = ggplot2::margin(0, 0, 2, 0)),
      axis.title.y       = ggplot2::element_text(margin = ggplot2::margin(0, 1, 0, 0.5)),
      axis.title.y.right = ggplot2::element_text(margin = ggplot2::margin(0, 0, 0, 2)),
      axis.text.x        = ggplot2::element_text(margin = ggplot2::margin(1, 0, 0, 0)),
      axis.text.x.top    = ggplot2::element_text(margin = ggplot2::margin(0, 0, 1, 0)),
      axis.text.y        = ggplot2::element_text(margin = ggplot2::margin(0, 1, 0, 0)),
      axis.text.y.right  = ggplot2::element_text(margin = ggplot2::margin(0, 0, 0, 1)),
      axis.ticks         = ggplot2::element_line(size = .5),
      axis.ticks.length  = ggplot2::unit(2, "pt"),
      axis.line          = ggplot2::element_line(size = .5),
      axis.line.x        = ggplot2::element_line(size = .5),
      axis.line.y        = ggplot2::element_line(size = .5),
      line               = ggplot2::element_line(size = .5),
      legend.margin      = ggplot2::margin(4, 4, 4, 4),
      legend.key.size    = ggplot2::unit(10, "pt"),
      legend.box.spacing = ggplot2::unit(6, "pt"),
      panel.spacing      = ggplot2::unit(1.5, "pt"),
      plot.title         = ggplot2::element_text(size = 12),
      plot.margin        = ggplot2::margin(1, 1, 1, 1),
      strip.background   = ggplot2::element_blank(),
      strip.text         = ggplot2::element_text(size = 10),
      strip.text.x       = ggplot2::element_text(margin = ggplot2::margin(3, 0, 3, 0)),
      strip.text.y       = ggplot2::element_text(margin = ggplot2::margin(0, 3, 0, 3))
    )
  )
}
