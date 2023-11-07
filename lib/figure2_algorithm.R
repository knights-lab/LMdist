# Figure 2: Algorithm visually explained
warning("Note that this script must be run in Rstudio to produce igraph figures properly!")

library(igraph)
library(plotrix)
library(vegan)
library(ggplot2)
library(cowplot)
library(ggpubr)
source('lib/lmdist_source.r')
set.seed(125)

# Small set of simulated data points
mini_sim <-read.table("data/simgradient/fake_rel_abun_long_n100.txt", row=1, header=T, sep="\t")
mini_sim <- mini_sim[c(1,12,19,31,43,50,61,72,80,92,100),]
mini_col_num <- factor(c(1,12,19,31,43,50,61,72,80,92,100))
mini_cols <- viridis::viridis(11, alpha=0.8)


# Before & after plots (PCoA)
mini_d <- vegdist(mini_sim, method="euclidean")
## FIGURE 2A: Before plot
mini_pc <- cmdscale(mini_d, k=2, eig=F)
mini_pc_df <- data.frame(PC1=mini_pc[,1], PC2=mini_pc[,2], mycolor=mini_col_num)
lim <- max(abs(mini_pc)) + 0.02
figure2a <- ggplot(mini_pc_df, aes(x=PC1, y=PC2, color=mycolor)) +
  geom_point(size = 8, pch = 16) +
  scale_color_manual(values=mini_cols) +
  xlim(-lim, lim) +
  ylim(-lim, lim) +
  labs(title="",x="PC 1",y="PC 2") +
  coord_fixed(ratio = 1) +
  theme_bw() +
  theme(axis.title = element_text(size=16),
        legend.position = "none")


## FIGURE 2F: After plot
mini_lmd <- lm.dist(mini_d, neighborhood.radius=0.18)
mini_pc_lmd <- cmdscale(mini_lmd, k=2, eig=F)
mini_lmd_df <- data.frame(PC1=mini_pc_lmd[,1], PC2=mini_pc_lmd[,2], mycolor=mini_col_num)
lim <- max(abs(mini_pc_lmd)) + 0.02
figure2f <- ggplot(mini_lmd_df, aes(x=PC1, y=PC2, color=mycolor)) +
  geom_point(size = 8, pch = 16) +
  scale_color_manual(values=mini_cols) +
  xlim(-lim, lim) +
  ylim(-lim, lim) +
  labs(title="",x="PC 1",y="PC 2") +
  coord_fixed(ratio = 1) +
  theme_bw() +
  theme(axis.title = element_text(size=16),
        legend.position = "none")


# iGraph plots of distance space and the connections made
## Make graph object from the points
##    note: using different radius then best radius for a more informative example
adj <- matrix(0, nrow(mini_sim), nrow(mini_sim))
for(i in 1:nrow(adj)) {
  adj[i,which(as.matrix(mini_d)[i,] <= 0.24823)] <- 1 
}
adj[adj>0] <- as.matrix(mini_d)[adj>0]
g <- graph.adjacency(adj, weighted=T, mode='undirected')
E(g)

## Prepare edges for one node
edge_colors3 <- rep("white",22)
edge_colors3[c(7,8,10:12)] <- "darkgrey"
edge_widths3 <- rep(0,22); edge_widths3[c(7,8,10:12)] <-3

## FIGURE 2B : graph w/ radius circle (no connections)
plot(g, vertex.label=NA, vertex.color=mini_cols, vertex.frame.color=mini_cols, vertex.size=22,
     edge.color="white", edge.width=edge_widths3)
draw.circle(x=0.02, y=-0.07, radius=0.7, lwd=3, lty="dashed")
text(x=0.5, y=-0.88, labels="radius = 0.25", cex=0.8)
figure2b <- recordPlot()

## FIGURE 2C : graph w/ radius circle & connections
set.seed(125)
plot(g, vertex.label=NA, vertex.color=mini_cols, vertex.frame.color=mini_cols, vertex.size=22,
     edge.color=edge_colors3, edge.width=edge_widths3)
draw.circle(x=0.02, y=-0.07, radius=0.7, lwd=3, lty="dashed")
text(x=0.5, y=-0.88, labels="radius = 0.25", cex=0.8)
figure2c <- recordPlot()

## FIGURE 2D : graph with all paths
set.seed(125)
plot(g, vertex.label=NA, vertex.color=mini_cols, vertex.frame.color=mini_cols, vertex.size=22,
     edge.color="darkgrey", edge.width=3)
figure2d <- recordPlot()

## FIGURE 2E : graph with path comparisons
edge_colors5 <- rep("white",22)
edge_colors5[c(4,7,12,19)] <- "darkgrey"
edge_widths5 <- rep(0,22); edge_widths5[c(4,7,12,19)] <- 3
print(paste0("original dist: ",round(as.matrix(mini_d)[2,11],3))) # original: 0.284
print(paste0("new dist: ",round(shortest.paths(g)[2,11],3)))      # LMdist: 0.783
set.seed(125)
plot(g, vertex.label=NA, vertex.color=mini_cols, vertex.frame.color=mini_cols, vertex.size=22,
     edge.color=edge_colors5, edge.width=edge_widths5)
lines(x=c(-0.9,0.71), y=c(-0.75,0.9), lwd=3, lty="dashed", col="black")
legend(x=-1.3, y=1.21, legend=c("original: 0.28", "LMdist: 0.78"), title="Inter-node Distance",
       col=c("black","darkgrey"), lty=c("dashed","solid"), lwd=2, bty="n", cex=0.8, y.intersp=1.2)
figure2e <- recordPlot()


# Putting it all together in one (high res tiff)
figure2b_gg <- ggdraw(figure2b, xlim=c(0.1,0.9), ylim=c(0.1,0.9), clip = "on")
figure2c_gg <- ggdraw(figure2c, xlim=c(0.1,0.9), ylim=c(0.1,0.9), clip = "on")
figure2d_gg <- ggdraw(figure2d, xlim=c(0.1,0.9), ylim=c(0.1,0.9), clip = "on")
figure2e_gg <- ggdraw(figure2e, xlim=c(0.1,0.9), ylim=c(0.1,0.9), clip = "on")
# tiff("figures/tif_files/figure2_algorithm.tif", width=9, height=6, unit="in", res=1200)
# cowplot::plot_grid(figure2a, figure2b_gg, figure2c_gg, figure2d_gg, figure2e_gg, figure2f,
#                   nrow = 2, labels="AUTO")
# dev.off()

# Low resolution version
tiff("figures/tif_files_low_res/figure2_algorithm_lowres.tif", width=9, height=6, unit="in", res=350)
cowplot::plot_grid(figure2a, figure2b_gg, figure2c_gg, figure2d_gg, figure2e_gg, figure2f,
                   nrow = 2, labels="AUTO")
dev.off()



