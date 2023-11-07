# Figure 4 : Swiss Roll dataset & Comparison to popular Manifold Learning
# note: best to run within Rstudio to see the 3D plot rendering

library(vegan)
library(Rtsne)
library(scales)
library(ggplot2)
library(plotly)
library(cowplot)
source("lib/lmdist_source.r")
set.seed(125)

# Loading Dataset
swiss_data <- read.delim("data/swiss_roll/swiss_roll.txt", sep=" ", header=F)
swiss_color <- c(read.delim("data/swiss_roll/swiss_roll_colors.txt", sep=" ", header=F)[,1])
swiss_euc <- vegdist(swiss_data, method="euclidean")
## note the denisty of distances
plot(density(swiss_euc), lwd=3, xlab="Distances", ylab="Density", main="Density of Swiss Dists")


# Helper function - 2D plot
plot2d <- function (df_2d, mytitle="", axes="PC", make_square=T) {
  df_2d <- as.data.frame(df_2d)
  colnames(df_2d) <- c("V1","V2")
  p <- ggplot(df_2d, aes(x=V1, y=V2)) +
    geom_point(size=3, pch=16, col=alpha(swiss_color, 0.8)) +
    labs(x=paste0(axes," 1"), y=paste0(axes," 2"), title=mytitle) +
    theme_bw() +
    theme(legend.position = "none",
          title = element_text(face="bold", size=12))
  if (make_square == TRUE) {
    mylim <- max(abs(as.matrix(df_2d)))
    p <- p +
      xlim(-mylim, mylim) +
      ylim(-mylim, mylim) + 
      coord_fixed(ratio=1)
  } else {
    mylim <- (max(df_2d[,1]) - min(df_2d[,1])) / 2
    p <- p +
      ylim(-(2/3)*mylim, (2/3)*mylim) + 
      coord_fixed(ratio=1.5)
  }
  return(p)
}


# Display in 3D
df_3d <- as.data.frame(swiss_data)
df_3d$cols <- factor(swiss_color, levels=unique(swiss_color))
fig5a <- plot_ly(data=df_3d, x=~V1, y=~V2, z=~V3, color=~cols, colors=unique(swiss_color))
fig5a <- fig5a %>% add_markers(marker=list(size=10))
fig5a <- fig5a %>% layout(scene = list(xaxis = list(title = 'Dim. 1'),
                                   yaxis = list(title = 'Dim. 2'),
                                   zaxis = list(title = 'Dim. 3')))
fig5a <- fig5a %>% layout(showlegend = FALSE)
fig5a

# Traditional PCA - 2D
swiss_pc <- as.data.frame(cmdscale(swiss_euc, k=2, eig=F))
fig5b <- plot2d(swiss_pc, "Swiss Roll PCA")

# Manifold Learning: t-SNE
swiss_tsne <- Rtsne(swiss_euc, dims=2, perplexity=60, initial_dims=3)
fig5c <- plot2d(swiss_tsne$Y, "Swiss Roll t-SNE", "Axis")

# Manifold Learning : Isomap
swiss_isomap <- isomap(swiss_euc, ndim=2, k=5)
fig5d <- plot2d(swiss_isomap$points, "Swiss Roll Isomap", "Axis", F)

# PCA w/ LMdist
swiss_lmdist <- lm.dist(swiss_euc, 6.2, phi=0)
swiss_lmd <- cmdscale(swiss_lmdist, k=2, eig=F)
fig5e <- plot2d(swiss_lmd, "Swiss Roll LMdist", "PC", F)

# Statistics : PC1 & manifold
manifold <- c(read.delim("data/swiss_roll/swiss_roll_manifold.txt", sep=" ", header=F)[,1])
adonis2(swiss_euc ~ manifold)   # F score: 5.51, p = 0.002
adonis2(swiss_lmdist ~ manifold) # F score: 3475.9, p = 0.001

# FIGURE 5 : Swiss roll dataset
## saving the 3d original plot bby hand in R viewer
warning("3D plot saved from python script for data creation.")
library(patchwork)
## arranging plots
s <- ggplot() + theme_void()
toprow <- cowplot::plot_grid(s, s, fig5b, s, labels=c("","A","B",""), nrow=1, rel_widths=c(0.5,1,1,0.5))
bottomrow <- cowplot::plot_grid(fig5c, fig5d, fig5e, labels=c("C","D","E"), nrow=1)
## high res
# tiff("figures/tif_files/figure5_swiss_roll.tif", width=8, height=5.5, unit="in", res=1200)
# cowplot::plot_grid(toprow, bottomrow, nrow=2)
# dev.off()
## low res
tiff("figures/tif_files_low_res/figure5_swiss_roll_lowres.tif", width=8, height=5.5, unit="in", res=350)
cowplot::plot_grid(toprow, bottomrow, nrow=2)
dev.off()


