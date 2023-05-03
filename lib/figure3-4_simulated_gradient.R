# Figure 3: Creation of Simulated dataset & Figure 4: application of LMdist
library(gtools)
library(scales)
library(vegan)
library(cowplot)
library(data.table)
source('lib/lgd_source.r')
source('lib/figurehelper_pcoaplot.R')
set.seed(125)

# Making of a simulated dataset (see create_simulation_data from previous code versions)
## parameters
g <- 1000     # length of gradient to sample (from 0)d
sd <- 150     # standard deviation of OTU coenoclines
skip <- sd/15 # size of jump between OTUs
n <- 50       # number of samples to get from gradient
conf <- 1     # confidence factor for dirichlet (larger -> closer to data)
## coenoclines
m <- 1000                                     # max gradient length
pad <- sd/skip
xvals <- seq(-(pad*sd), m+(pad*sd), by=skip)[-1]  # ensures 0 to m in figure is flat
x <- sapply(xvals, function (x) {dnorm(seq(0,m,1), x, sd)*100})  # create null distributions
x[abs(x) < 0.00001] <- 0                      # remove otu entirely at tails
x <- sweep(x, 1, rowSums(x), "/")             # normalize rows (i.e. possible sample combinations)
x <- x[,-which(colSums(x) == 0)]              # remove buffer otus
## sample coenoclines
inds <- round(seq(0, g, length.out=(n+2)))[-c(1,n+2)] # evenly selects samples along distribution
dat <- as.data.frame(x[inds,])
dat <- dat[,-which(colSums(dat) == 0)]
rownames(dat) <- paste0("sample.", rownames(dat))
colnames(dat) <- paste0("otu.", 1:ncol(dat))
## add dirichlet noise
noise <- as.data.frame(t(apply(dat, 1, function (x) colMeans(rdirichlet(100, x*conf))))) # use each sample as its own prior
colnames(noise) <- colnames(dat)
rownames(noise) <- paste0("noise.",seq(1,nrow(noise),1))
dat <- rbind(dat, noise)


# Coenocline visualization
nv <- ncol(x)/18   # every nv-th otu coenocline
colors <- c("#440154FF",viridis::viridis(((ncol(dat)-pad) / nv)-2))
rx <- x[,colSums(dat) > 0]
rx[g:nrow(rx),] <- NA
rx <- rx[,seq(pad, ncol(rx)-pad, by=nv)]
rx[rx == 0] <- NA
par(mgp=c(1,0,0))
matplot(rx, type="l", col=colors, lty="solid", lwd=10, xlim=c(0,m),
        ylab="relative abundances", xlab="gradient",
        cex.lab=1.1, xaxt="n", yaxt="n")
rect(200,-0.0001,250,0.0272, border="black", lwd=4)
text(x=225, y=0.029, labels="sample.10", xpd=T, cex=0.8)
figure3a <- recordPlot()
par(mgp=c(3,1,0))


# Sample heatmap visualization
h_idx <- rep(50:1, each=2)
h_idx <- h_idx + c(0,50)
# mat <- as.matrix(dat[h_idx,])
# heatmap(mat, Rowv=NA, Colv=NA, revC=F,
#         col=c("white",viridis::viridis(44)))
tmp <- dat[h_idx,]
tmp$id <- rownames(tmp)
suppressWarnings(mat_long <- melt(tmp))
mat_long$id <- factor(mat_long$id, levels=unique(mat_long$id))
mat_long$color <- as.numeric(gsub("otu.", "", mat_long$variable))
mat_long$value_alpha <- mat_long$value
mat_long$value_alpha[mat_long$value_alpha <= .Machine$double.eps/10] <- NA
mat_long$value_alpha <- mat_long$value_alpha / max(mat_long$value_alpha, na.rm=T)
mat_long$value_alpha <- cut(mat_long$value_alpha, breaks=c(0,0.01,0.05,0.1,0.2,0.3,0.4,0.5,1))
figure3b <- ggplot(mat_long, aes(x=variable, y=id, fill=color))+
  geom_tile(aes(alpha=value_alpha)) +
  labs(x="", y="", fill="gradient") +
  scale_x_discrete(breaks = c("otu.1","otu.50","otu.100","otu.150","otu.200", paste0("otu.", max(mat_long$color)))) +
  scale_y_discrete(breaks = c("sample.1","sample.10","sample.20","sample.30","sample.40","sample.50")) +
  scale_fill_viridis_c(breaks=c(0,max(mat_long$color)), labels=c("","")) +
  scale_alpha_discrete(range=c(0.5,1), na.value=0) +
  guides(alpha = "none",
         fill = guide_colourbar(barwidth = 8, barheight = 1)) +
  theme_classic() +
  theme(plot.background=element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1, vjust=1),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        legend.position = "bottom",
        legend.title = element_text(vjust=0.85),
        legend.margin=margin(grid::unit(0, "cm")),
        plot.margin=margin(0.5, 0.5, 0.1, 0.1, "cm"))
figure3b <- figure3b + coord_cartesian(xlim=c(1,231), clip="off") +
  annotate("rect", xmin = -60, xmax = 235, ymin = 78, ymax = 84, col="black", fill=alpha("white",0))

# FIGURE 3 : Combine Coenoclines visualization & heatmap
figure3a_gg <- ggdraw(figure3a)
## high resolution
tiff("figures/tif_files/figure3_simdata_creation.tif", width=8, height=4, units="in", res=1200)
cowplot::plot_grid(figure3a_gg, figure3b, nrow=1, labels="AUTO", rel_widths=c(1.2,1))
dev.off()
## low resolution
tiff("figures/tif_files_low_res/figure3_simdata_creation_lowres.tif", width=8, height=4, units="in", res=350)
cowplot::plot_grid(figure3a_gg, figure3b, nrow=1, labels="AUTO", rel_widths=c(1.2,1))
dev.off()



# PCoA with a typical distance metric (Bray Curtis for now)
d <- vegdist(dat, method="bray")
pc <- cmdscale(d, k=10, eig=T)
df <- as.data.frame(pc$points)
colnames(df) <- paste0("PC",1:ncol(df))
df$gradient <- rep(1:50, 2)
var <- (pc$eig / sum(pc$eig)) * 100
# plot(1:10, var[1:10], type="b", main="Stress Plot", xlab="PC", ylab="Stress")
mylim <- max(abs(as.matrix(df[,1:2])))
figure4a <- ggplot(df, aes(x=PC1, y=PC2, col=gradient)) +
  geom_point(size=5, pch=20) +
  scale_color_viridis_c(alpha = 0.8, breaks=c(0,50), labels=c("","")) +
  xlim(-mylim, mylim) +
  ylim(-mylim, mylim) +
  labs(title="Simulated Gradient\nPCoA (Bray-Curtis)",
       x=paste0("PC 1 [", round(var[1],1), "%]"),
       y=paste0("PC 2 [", round(var[2],1), "%]")) +
  guides(col = guide_colorbar(barwidth = 8, barheight = 1)) +
  coord_fixed(ratio = 1) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_text(vjust=0.85, face="bold", size=12),
        title = element_text(face="bold", size=12),
        axis.title = element_text(size=10))
adonis2(d ~ gradient, data = data.frame(id=rownames(dat), gradient=rep(1:50, 2))) # F score: 120, p = 0.001
cor.test(pc$points[,1], rep(1:50, 2)) # cor = 0.97, p < 2.2e-16


# Compare community distances with gradient distances
cg <- data.frame(comm=c(d), grad=c(vegdist(df$gradient, "euclidean")))
figure4b <- ggplot(cg, aes(x=grad, y=comm)) +
  geom_point(size=2, color=alpha("blue", 0.2)) +
  geom_smooth(color="black", size=1, method="gam") +
  labs(x="Gradient Distances\n(Euclidean)", y="Community Distances\n(Bray-Curtis)", color="") +
  theme_bw() +
  theme(axis.title = element_text(size=12),
        axis.text = element_text(size=10))
xdens4b <- axis_canvas(figure4b, axis = "x") +
  geom_density(data = cg, aes(grad),
               alpha = 0.4, size = 0.8, color = "darkgrey", fill="darkgrey") +
  theme(plot.margin = margin(0.2,0,0,0, "cm"))
ydens4b <- axis_canvas(figure4b, axis = "y", coord_flip = TRUE) +
  geom_density(data = cg, aes(comm),
               alpha = 0.4, size = 0.8, color="blue", fill="blue") +
  coord_flip()
figure4b <- insert_xaxis_grob(figure4b, xdens4b, grid::unit(.2, "null"), position = "top")
figure4b <- insert_yaxis_grob(figure4b, ydens4b, grid::unit(.2, "null"), position = "right")
plot(figure4b)


# After LMdist : PCoA
lgd <- lg.dist(d, 0.4)
pc_lgd <- cmdscale(lgd, k=10, eig=T)
df_lgd <- as.data.frame(pc_lgd$points)
colnames(df_lgd) <- paste0("PC",1:ncol(df_lgd))
df_lgd$gradient <- rep(1:50, 2)
var_lgd <- (pc_lgd$eig / sum(pc_lgd$eig)) * 100
# plot(1:10, var_lgd[1:10], type="b", main="Stress Plot (After LMdist)", xlab="PC", ylab="Stress")
mylim <- max(abs(as.matrix(df_lgd[,1:2])))
figure4c <- ggplot(df_lgd, aes(x=PC1, y=PC2, col=gradient)) +
  geom_point(size=5, pch=20) +
  scale_color_viridis_c(alpha = 0.8, breaks=c(0,50), labels=c("","")) +
  xlim(-mylim, mylim) +
  ylim(-mylim, mylim) +
  labs(title="Simulated Gradient\nPCoA (LMdist-adjusted)",
       x=paste0("PC 1 [", round(var_lgd[1],1), "%]"),
       y=paste0("PC 2 [", round(var_lgd[2],1), "%]")) +
  guides(col = guide_colorbar(barwidth = 8, barheight = 1)) +
  coord_fixed(ratio = 1) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_text(vjust=0.85, face="bold", size=12),
        title = element_text(face="bold", size=12),
        axis.title = element_text(size=10))
adonis2(lgd ~ gradient, data = data.frame(id=rownames(dat), gradient=rep(1:50, 2))) # F score: 1778.6, p = 0.001
cor.test(pc_lgd$points[,1], rep(1:50, 2))  # cor = 0.999, p < 2.2e-16


# After LGD : Compare Community & Gradient Distances
cg_lgd <- data.frame(comm=c(lgd), grad=c(vegdist(df_lgd$gradient, "euclidean")))
figure4d <- ggplot(cg_lgd, aes(x=grad, y=comm)) +
  geom_point(size=2, color=alpha("orange", 0.2)) +
  geom_smooth(color="black", size=1, method="gam") +
  labs(x="Gradient Distances\n(Euclidean)", y="Community Distances\n(LMdist-adjusted)", color="") +
  theme_bw() +
  theme(axis.title = element_text(size=12))
xdens4d <- axis_canvas(figure4d, axis = "x") +
  geom_density(data = cg_lgd, aes(grad),
               alpha = 0.4, size = 1, color = "darkgrey", fill="darkgrey")
ydens4d <- axis_canvas(figure4d, axis = "y", coord_flip = TRUE) +
  geom_density(data = cg_lgd, aes(comm),
               alpha = 0.4, size = 1, color="orange", fill="orange") +
  coord_flip()
figure4d <- insert_xaxis_grob(figure4d, xdens4d, grid::unit(.2, "null"), position = "top")
figure4d <- insert_yaxis_grob(figure4d, ydens4d, grid::unit(.2, "null"), position = "right")
plot(figure4d)


# FIGURE 4: Applying LMdist to simulated gradient
## high resolution
tiff("figures/tif_files/figure4_simdata_lmdist.tif", width=8, height=8, units="in", res=1200)
cowplot::plot_grid(figure4a, figure4b, figure4c, figure4d, nrow=2, labels="AUTO")
dev.off()
## low resolution
tiff("figures/tif_files_low_res/figure4_simdata_lmdist_lowres.tif", width=8, height=8, units="in", res=350)
cowplot::plot_grid(figure4a, figure4b, figure4c, figure4d, nrow=2, labels="AUTO")
dev.off()




