# Helper function : Before and After PCoA plots
library(ggplot2)

pcoa_plot <- function (dist, colors, fills, shapes, mysize=3, color_col, fill_col, shape_col,
                       maintitle="", legendtitle="", legendtitle2="", flipPC1=F, flipPC2=F, discrete=T) {
  pc <- cmdscale(dist, k=3, eig=T)
  colnames(pc$points) <- paste0("PC", 1:ncol(pc$points))
  var <- (pc$eig / sum(pc$eig)) * 100
  df <- cbind(pc$points, data.frame(color=color_col, fill=fill_col, shape=shape_col))
  if (flipPC1) { df$PC1 <- df$PC1 * -1 }
  if (flipPC2) { df$PC2 <- df$PC2 * -1 }
  lim <- max(abs(as.matrix(df[,1:2]))) # to make output a square
  p <- ggplot(df, aes(x=PC1, y=PC2, color=color, fill=fill, shape=shape)) +
    geom_point(size = mysize) +
    scale_color_manual(name=legendtitle, values=colors) +
    scale_fill_manual(name=legendtitle, values=fills) +
    xlim(-lim, lim) +
    ylim(-lim, lim) +
    labs(title=maintitle,
         x=paste0("PC 1 [",round(var[1],1),"%]"),
         y=paste0("PC 2 [",round(var[2],1),"%]")) +
    coord_fixed(ratio = 1) +
    theme_bw() +
    theme(title = element_text(face="bold", size=14),
          axis.title = element_text(size=12),
          legend.title = element_text(face="bold", size=12),
          legend.text = element_text(size=10),
          legend.spacing = unit(0.2,"cm"))
  # Condition for shape included or not
  if (legendtitle2 == "") {
    p <- p + scale_shape_manual(name=legendtitle, values=shapes)
  } else {
    p <- p + scale_shape_manual(name=legendtitle2, values=shapes)
  }
  # Condition for discrete or continuous color scale
  if (discrete == F) {
    suppressWarnings(suppressMessages( p <- p +
      scale_color_viridis(name=legendtitle, option="plasma") +
      scale_fill_viridis(name=legendtitle, option="plasma")))
      #  + scale_shape_manual(guide="none", values=shapes)
  }
  return(p)
}
