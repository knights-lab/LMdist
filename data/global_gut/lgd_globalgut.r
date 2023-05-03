# lgd_humangut.r
# Local Gradient Distance on Human Global Gut Data
# Knights Lab - University of Minnesota
# September 2019
# usage : lgd_humangut.r

##### Set Up #####
source('/project/flatiron2/suzie/detrending/fake/lgd_source.r')
current_dir <- "/project/flatiron2/suzie/detrending/humangut/"
library(geiger)
library(phyloseq)


##### Load Data #####
dat_hg <- read.delim(paste0(current_dir, "processed_data/421_clean_otus.txt"),
                     header=T, sep="\t", quote="\"", row.names=1, comment.char="", as.is=T)
meta_hg <- read.delim(paste0(current_dir, "processed_data/clean_map.txt"),
                   header=T, sep="\t", quote="\"", comment.char="", as.is=T)
pop_cols <- c("red", "blue", "green") # malawians, usa, venezuelans


##### Get Distances #####
# Bring in unweighted UniFrac matrix
d_hg <- read.delim(paste0(current_dir, "processed_data/clean_uw_unifrac_dist.txt"),
                   header=T, row.names=1, sep="\t", as.is=T)
# Convert to distance object
d_hg <- as.dist(as.matrix(d_hg))

##### Local Gradient Distance #####
# note: try just U.S. to uncover age
# note: try community with different age groups
# determine reasonable neighborhood radius
plot(density(d_hg))
r <- 0.4
# calculate local gradient distances
cat('Calculating local gradient distance...\n')
lgd_hg <- lg.dist(d_hg, neighborhood.size = 15, weighted=TRUE)
cat('Calculating PCoA of original distances...\n')
pc.d_hg <- cmdscale(d_hg)
cat('Calculating PCoA of transformed distances...\n')
pc.lgd_hg <- cmdscale(lgd_hg)


##### Plots with Meta Data #####
# Geography groupings
cat('Plotting PCoA of original distances by geography...\n')
df_o <- data.frame(X=pc.d_hg[,1], Y=pc.d_hg[,2], age=meta_hg$age, body_site=meta_hg$body_site,
                   sex=meta_hg$sex, env_biome=meta_hg$env_biome, geo_loc=meta_hg$geo_loc_name)
ggplot(df_o, aes(x=X, y=Y, color=factor(geo_loc))) +
  geom_point(alpha=0.7, size=4) +
  scale_color_manual(values=pop_cols) +
  labs(title="Original Distances",
       x="PC1",
       y="PC2") +
  guides(color = guide_legend("population")) +
  xlim(range(pc.d_hg)) + ylim(range(pc.d_hg)) +
  theme_classic() + NULL
cat('Plotting PCoA of transformed distances by geography...\n')
df_t <- data.frame(X=pc.lgd_hg[,1], Y=pc.lgd_hg[,2], age=meta_hg$age, body_site=meta_hg$body_site,
                   sex=meta_hg$sex, env_biome=meta_hg$env_biome, geo_loc=meta_hg$geo_loc_name)
ggplot(df_t, aes(x=X, y=Y, color=factor(geo_loc))) +
  geom_point(alpha=0.7, size=4) +
  scale_color_manual(values=pop_cols) +
  labs(title="Transformed Distances",
       x="PC1",
       y="PC2") +
  guides(color = guide_legend("population")) +
  xlim(range(pc.lgd_hg)) + ylim(range(pc.lgd_hg)) +
  theme_classic() + NULL
# Age Groupings
cat('Plotting PCoA of original distances by age...\n')
ggplot(df_o, aes(x=X, y=Y, color=age)) +
  geom_point(alpha=0.7, size=4) +
  scale_color_viridis_c() +
  labs(title="Original Distances",
       x="PC1",
       y="PC2",
       color="age") +
  xlim(range(pc.d_hg)) + ylim(range(pc.d_hg)) +
  theme_classic() + NULL
cat('Plotting PCoA of transformed distances by age...\n')
ggplot(df_t, aes(x=X, y=Y, color=as.numeric(age))) +
  geom_point(alpha=0.7, size=4) +
  scale_color_viridis_c() +
  labs(title="Transformed Distances",
       x="PC1",
       y="PC2",
       color="age") +
  xlim(range(pc.lgd_hg)) + ylim(range(pc.lgd_hg)) +
  theme_classic() + NULL

##### Break into groups and re-plot #####