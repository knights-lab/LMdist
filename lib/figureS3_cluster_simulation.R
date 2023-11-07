# Supplemental Figure 3: Clusters simulation & more whittaker data

library(stats)
library(plot.matrix)
library(vegan)
library(scales)
library(gtools)
library(ggplot2)
source("lib/lmdist_source.r")
set.seed(125)


##### Scratch Work: Playing around with clusters #####
# what do real samples look like? Let's consider the Soil dataset
soil <- read.delim("data/soil/44766_clean_otus_norm.txt", sep="\t", row=1, header=T) # loading soil dataset
soil <- soil[,runif(10, min=1, max=ncol(soil))]
soil <- soil[rowSums(soil) > 0,]
plot(density(soil[,1]), main="Soil Samples", type="l", lwd=1, col=alpha("blue",0.4), ylim=c(0,800))
for (c in 2:10) {
  lines(density(soil[,c]), lwd=1, col=alpha("blue",0.4))
}
legend(x=0.01, y=450, legend=c("Sample Distribution"), lwd=1, col="blue")
# can we emulate this with a beta distribution?
len_x <- seq(0,1,length.out=50)
prior <- dbeta(len_x, 1, 20)
plot(len_x, prior, type="l", lwd=3, main="Beta Distribution & Dirichlet Samples",
     xlab="Proportion of Sample", ylab="Frequency") # curve of prior distribution used
mysample <- rbeta(len_x, 1, 20) # sample the same curve
lines(density(mysample), lty="dashed", lwd=3, col="purple") # distribution of sample from beta
# what about other samples if we use this as input to dirichlet?
conf <- 10
mydirichlet <- rdirichlet(20, mysample*conf) # 20 dirichlet samples
mydirichlet[mydirichlet < .Machine$double.eps] <- 0
rowSums(mydirichlet)
for (row in 1:20) {
  lines(density(mydirichlet[row,],n=16), pch=16, cex=2, col=alpha("purple",0.2))
}
legend(x=0.4, y=12, legend=c("Beta(1,20)","Random Beta Sample","Dirichlet Samples"),
       lty=c("solid","dashed","solid"), lwd=c(3,3,1), col=c("black","purple","purple"))
# make a sample and visualize it
myset <- rbind(mysample/sum(mysample), mydirichlet)
heatmap(myset, col=topo.colors(254), Rowv=NA)
legend(x="right", legend=c("low","  |","  |","  |","high"), fill=topo.colors(5))



##### Function to create clusters dataset #####
# Create nclust clusters from nsamp samples
## note: balanced group sizes, but could add another parameter for inbalanced groups
create_clusters <- function (minsamp=90, nclust=3, n_otu=500, beta_ab=c(1,20)) {
  # What to do if overlap = T vs overlap = F??
  # Initialize variables needed
  nsamp_per <- ceiling(minsamp / nclust)    # number of samples per cluster
  len_x <- seq(0, 1, length.out=n_otu)      # x values for the beta function
  conf <- 10                                # confidence factor for dirichlet, currently built-in
  matrix <- matrix(NA, nrow=0, ncol=n_otu)  # empty matrix to fill with samples
  # Use beta distribution as prior to dirichlet of samples
  for (c in 1:nclust) {
    prior <- rbeta(len_x, beta_ab[1], beta_ab[2])
    dir_samples <- rdirichlet(nsamp_per, prior * conf)
    dir_samples[dir_samples < .Machine$double.eps] <- 0 # below machine epsilon set to 0
    matrix <- rbind(matrix, dir_samples)
  }
  # Clean up matrix and return
  matrix <- matrix[rowSums(matrix) > 0,]  # remove empty rows (should not happen, but safety check)
  matrix <- matrix[,colSums(matrix) > 0]  # remove empty columns
  return(matrix)
}
# Plotting output
mymat <- create_clusters()
heatmap(mymat, Rowv=NA, col=c("white",topo.colors(49)))
supp3a <- recordPlot()

# Ordination with these clusters
my_distances <- c("bray", "jaccard", "aitchison", "euclidean", "manhattan", "gower", "chisq", "canberra")
par(mfrow=c(2,4),mar=c(2,2,5,1), mgp=c(0.5,0,0))
## Default parameters, all distances
for (m in my_distances) {
  if (m == "aitchison") d <- vegdist(mymat+0.000000001, method=m)
  else d <- vegdist(mymat, method=m)
  pc <- cmdscale(d, k=2)
  plot(pc[,1], pc[,2], pch=16, cex=2, col=alpha(rep(c("purple", "orange", "blue"), each=30),0.6),
       xlab="PC 1", ylab="PC 2", asp=1, xaxt="n", yaxt="n")
  title(m, adj = 0.5, line = 1)
}
mtext("Default parameters", side = 3, line = -2, outer = TRUE, cex=1.5)
par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1), mgp=c(3,1,0))


# Distribution of distances
par(mfrow=c(2,4),mar=c(3,3,5,1), mgp=c(1.5,0.5,0))
for (m in my_distances) {
  if (m == "aitchison") {
    d <- vegdist(mymat+0.000000001, method=m)
  } else {
    d <- vegdist(mymat, method=m)
  }
  plot(density(d), lwd=3, col="black", main="", xlab="Value", ylab="Freq")
  title(m, adj=0.5, line=1)
}
mtext("Distance Distributions (default params)", side = 3, line = -2, outer = TRUE, cex=1.5)
par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1), mgp=c(3,1,0))


# Compare just Aitchison distance and LMdist
aitch <- vegdist(mymat+0.000000001, method="aitchison")
pc_aitch <- cmdscale(aitch, k=2)
lmd <- lm.dist(aitch, smooth=T) # chose not to adjust
pc_lmd <- cmdscale(lmd, k=2)
## plotting
simple_2dplot <- function (mypc, mytitle) {
  df <- as.data.frame(mypc)
  colnames(df) <- c("PC1", "PC2")
  df$Group <- rep(c("A", "B", "C"), each=30)
  p <- ggplot(df, aes(x=PC1, y=PC2, col=Group)) +
    geom_point(size=3, pch=16) +
    scale_color_manual(values=c("purple","orange","blue")) +
    labs(x="PC 1", y="PC 2", title=mytitle) +
    coord_fixed(ratio=1) +
    theme_bw() +
    theme(title = element_text(face="bold", size=12),
          axis.title = element_text(size=10),
          legend.position = "none")
  return(p)
}
supp3b <- simple_2dplot(pc_aitch, "Sim. Clusters (Aitchison)")
supp3c <- simple_2dplot(pc_lmd, "Sim. Clusters (LMdist)")



##### Real Cluster examples #####
# Use figure helper for plots
source("lib/figurehelper_pcoaplot.R")

# Whitaker Table 3
whit <- as.matrix(read.delim("data/smokey/Whittaker - Table3.tsv", sep="\t", row=1, header=T))
whit[is.na(whit)] <- 0
whit_d <- vegdist(t(whit), "bray")
moisture <- data.frame(Station = paste0("S", 1:12), moisture_gradient = 1:12)
## Before (Bray-Curtis)
whit2_col <- viridis::plasma(12)
whit2_shp <- rep(16, 12)
mg <- factor(moisture$moisture_gradient, levels=c(1,2,3,4,5,6,7,8,9,10,11,12))
whit2_before <- pcoa_plot(whit_d, whit2_col, whit2_col, whit2_shp, 4,
                          mg, mg, mg,
                          "Whittaker Table 3 (Bray-Curtis)", "Moisture", discrete=T)
## After (LMdist)
whit2_lmd <- lm.dist(whit_d) # chooses 0.633
whit2_after <- pcoa_plot(whit2_lmd, whit2_col, whit2_col, whit2_shp, 4,
                         mg, mg, mg,
                         "Whittaker Table 3 (LMDist-adjusted)", "Moisture", discrete=T)


# Western Immigration (Vangay et al)
## loading dataset
imp <- read.delim('data/imp/unweighted_unifrac_dm.txt', sep="\t", row=1, header=T)
meta_imp <- read.delim('data/imp/map.txt', sep="\t", header=T)
rownames(meta_imp) <- meta_imp$SampleID
meta_imp$Sample.Group <- as.factor(meta_imp$Sample.Group)
meta_imp <- meta_imp[rownames(meta_imp) %in% rownames(imp),]
imp <- imp[rownames(meta_imp),rownames(meta_imp)]
meta_imp$Sample.Group <- factor(meta_imp$Sample.Group, levels=c("Control","Hmong1st","Hmong2nd", "HmongThai", "Karen1st", "KarenThai"),
       labels=c("U.S. Control", "Hmong 1st", "Hmong 2nd", "Hmong Thai", "Karen 1st", "Karen Thai"))
## colors and shapes
imp_col <- c("#e38cbd", "#7e7e7e", "#d65ca0", "#d65ca0", "#4ea99b", "#4ea99b")
names(imp_col) <- c("Hmong 2nd", "U.S. Control", "Hmong Thai", "Hmong 1st", "Karen Thai", "Karen 1st")
imp_shp <- c(17, 17, 16, 1, 16, 1)
names(imp_shp) <- c("Hmong 2nd", "U.S. Control", "Hmong Thai", "Hmong 1st", "Karen Thai", "Karen 1st")
## Before (Unw. Unifrac)
imp_before <- pcoa_plot(imp, imp_col, imp_col, imp_shp, 2,
                        meta_imp$Sample.Group, meta_imp$Sample.Group, meta_imp$Sample.Group,
                        "Immigration (Unw. UniFrac)", "Group", discrete=T)
## After (LMdist)
imp_lmd <- lm.dist(imp, epsilon=0.1)
imp_after <- pcoa_plot(imp, imp_col, imp_col, imp_shp, 2,
                       meta_imp$Sample.Group, meta_imp$Sample.Group, meta_imp$Sample.Group,
                       "Immigration (LMdist-adjusted)", "Group", discrete=T)


# Global Gut study (Yatsunenko et al)
## Load the global gut dataset
dat_hg <- read.delim("data/global_gut/421_clean_otus.txt", header=T, sep="\t", quote="\"", row.names=1, comment.char="", as.is=T)
meta_hg <- read.delim("data/global_gut/clean_map.txt",header=T, sep="\t", quote="\"", comment.char="", as.is=T)
rownames(meta_hg) <- meta_hg$Sample_ID
pop_cols <- c("red", "blue", "green"); names(pop_cols) <- c("Malawi", "USA", "Venezuela");
## Load the unifrac distances
d_hg <- read.delim("data/global_gut/clean_uw_unifrac_dist.txt", header=T, row.names=1, sep="\t", as.is=T)
d_hg <- as.dist(as.matrix(d_hg))
identical(rownames(as.matrix(d_hg)), meta_hg$Sample_ID)
## Before (unw. UniFrac)
gg_before <- pcoa_plot(d_hg, pop_cols, pop_cols, rep(16, 3), 2,
                       meta_hg$geo_loc_name, meta_hg$geo_loc_name, meta_hg$geo_loc_name,
                       "Global Gut (Unw. UniFrac)", "Population", discrete=T)
## After (LMdist)
hg_lmd <- lm.dist(d_hg)
gg_after <- pcoa_plot(hg_lmd, pop_cols, pop_cols, rep(16, 3), 2,
                      meta_hg$geo_loc_name, meta_hg$geo_loc_name, meta_hg$geo_loc_name,
                      "Global Gut (LMdist-adjusted)", "Population", discrete=T)


##### FIGURE SUPPLEMENTAL 3 #####
library(cowplot)
library(gridGraphics)
supp3a <- ggdraw(supp3a)
## combined plot
left <- cowplot::plot_grid(supp3a, supp3b, supp3c, nrow=3, rel_heights=c(2,1,0.5), labels="AUTO")
right <- cowplot::plot_grid(imp_before, imp_after, gg_before, gg_after, whit2_before, whit2_after,
                            nrow=3, ncol=2, labels=c("D","","E","","F",""))
## high res
# tiff("figures/tif_files/suppfig3_clusters.tif", width=14, height=12, unit="in", res=1200)
# cowplot::plot_grid(left, right, nrow=1, rel_widths=c(1,2))
# dev.off()
## low res
tiff("figures/tif_files_low_res/suppfig3_clusters_lowres.tif", width=14, height=12, unit="in", res=350)
cowplot::plot_grid(left, right, nrow=1, rel_widths=c(1,2))
dev.off()


# Scratch : Beta distribution plots
# --- Set up x values
# len_x <- seq(0,1,by=0.025)    # x values for the beta distributions
# --- Plotting the beta distribution with various values
# plot(len_x, dbeta(len_x, 0.5,1), type="l", col="blue", lwd=3, ylab="Prob Dist", xlab="X", main="PDF of Beta", ylim=c(0,3.5))
# text(0.17, 1.7, label="beta(0.5,1)", col="blue")
# lines(len_x, dbeta(len_x, 0.5,0.5), col="red", lwd=3)
# text(0.58, 0.5, label="beta(0.5,0.5)", col="red")
# lines(len_x, dbeta(len_x, 1,1), col="purple", lwd=3)
# text(0.82, 1.2, label="beta(1,1)", col="purple")
# lines(len_x, dbeta(len_x, 3,3), col="darkgreen", lwd=3)
# text(0.61, 1.95, label="beta(3,3)", col="darkgreen")
# lines(len_x, dbeta(len_x, 3,1.5), col="cyan2", lwd=3)
# text(0.8, 2, label="beta(3,1.5)", col="cyan2")
# lines(len_x, dbeta(len_x, 6,8), col="orange", lwd=3)
# text(0.55, 2.75, label="beta(6,8)", col="orange")
# lines(len_x, dbeta(len_x, 3,10), lwd=3, lty="dashed")
# text(0.32, 3.2, label="beta(3,10)")
# --- Output of the the beta functions
# par(mar=c(4,2,2,2), mfrow=c(2,2))
# plot(dbeta(len_x, 6,8), main="Density: dbeta(6,8)", pch=16, cex=2, xlab="", ylab="")
# plot(pbeta(len_x, 6,8), main="Distribution: pbeta(6,8)", pch=16, cex=2, xlab="", ylab="")
# plot(qbeta(len_x, 6,8), main="Quantile Func: qbeta(6,8)", pch=16, cex=2, xlab="", ylab="")
# plot(density(rbeta(len_x, 6,8)), main="Random Deviates: rbeta(6,8)", lwd=3, xlab="", ylab="")
# par(mar=c(5.1,4.1,4.1,2.1), mfrow=c(1,1))

