# Figure 5: Case Studies

library(ggplot2)
library(vegan)
library(reshape2)
library(viridis)
source("lib/lgd_source.r")
library(egg)
set.seed(125)

##### Loading Datasets #####
# Soil dataset (88 soils)
meta_soil <- read.table("data/soil/clean_map.txt", header=T, sep="\t")
rownames(meta_soil) <- meta_soil$SampleID
soil_n <- read.table("data/soil/44766_clean_otus_norm.txt", row=1, header=T, sep="\t")
meta_soil$pH <- cut(meta_soil$ph, breaks = c(0,4:8,14))
levels(meta_soil$pH) <- c("<4","4-5","5-6","6-7","7-8",">8")

# Guerrero Negro dataset
meta_gn <- read.table("data/guerrero_negro/clean_map.txt", header=T, sep="\t")
rownames(meta_gn) <- meta_gn$SampleID
meta_gn$depth_mm <- c("0-1","1-2","2-3","3-4","4-5","5-6","6-10","10-22","22-34")[as.factor(meta_gn$depth)]
gn_n <- read.table("data/guerrero_negro/47908_clean_otus_norm.txt", row=1, header=T, sep="\t")
meta_gn$depth <- as.factor(meta_gn$depth)
meta_gn$depth_mm <- factor(meta_gn$depth_mm, levels=unique(meta_gn$depth_mm))

# Cecum dataset
cecum_d <- readRDS("data/cecum/cecum_bray_dist.rds") # Bray-Curtis distances
meta_cecum <- read.table("data/cecum/cecum_meta.csv", sep=",", row=1, header=T)
## remove 7 samples (which were removed for the paper)
remove <- c("WPC2FHBD15B06C","WPC2FHBD18B07C","WPC2FHBD15B05C","WPC2FD15B04C",
            "WPC2FD15B05C","WPC2FD08B08C","WPC2FD01B10C")
remove_idx <- which(meta_cecum$SeqID %in% remove)
## also remove Day 1 samples since they were so different
remove_idx <- c(remove_idx, which(meta_cecum$Collection == "D01"))
## reformatting
cecum_d <- as.dist(as.matrix(cecum_d)[-remove_idx,-remove_idx])
meta_cecum <- meta_cecum[-remove_idx,]
meta_cecum$Collection <- as.factor(meta_cecum$Collection)
levels(meta_cecum$Collection) <- gsub("D", "Day ", levels(meta_cecum$Collection))
meta_cecum$System <- as.factor(meta_cecum$System)
levels(meta_cecum$System) <- c("Hatch Brood", "Pen")

# Whittaker : Herbs and Shrubs
## Loading Herbs & Locations dataset suggested by Peter Kennedy
locs <- read.csv("data/herbs/1_Sampling_Location_List.csv", header=T)
herbs <- read.csv("data/herbs/3_Herb_Shrub_Data.csv", header=T)
herb_c <- dcast(herbs, Sample.Number ~ Species.Name, value.var="First_Value", fill=0, fun.aggregate=mean)
rownames(herb_c) <- herb_c$Sample.Number
herb_c$Sample.Number <- NULL
herb_c <- herb_c[rownames(herb_c) %in% locs$Station_Number,]
locs <- locs[match(rownames(herb_c), locs$Station_Number),]
removeherb <- which(is.na(locs$Elevation..m.)) # remove NA values for elevation
removeherb <- c(removeherb, which(locs$Soil_Type != "D")) # remove non-diorite soils (so elevation is just within diorite soils)
removeherb <- c(removeherb, which(locs$Station_Number %in% c(311,241,238,245,140,254,239,154,237,240,89))) # remove some uncharacteristic outliers (under/oversampled)
herb_c <- herb_c[-removeherb,]
locs <- locs[-removeherb,]
herb_d <- vegdist(herb_c, "bray")
## Create simplified site description
site_desc_groups <- tolower(locs$Site_Description)
site_desc_groups[site_desc_groups %in% c("ravine sheltered slope","sheltered slope")] <- "sheltered slope"
site_desc_groups[site_desc_groups %in% c("open clope","open slope")] <- "open slope"
site_desc_groups[site_desc_groups == ""] <- NA
site_desc_groups[!(site_desc_groups %in% c("sheltered slope","open slope"))] <- "other"
locs$Site_Group <- site_desc_groups



##### Helper Functions : Before & After plots #####
source("lib/figurehelper_pcoaplot.R")
pcoa_plotf6 <- function (dist, colors, fills, shapes, mysize=3, color_col, fill_col, shape_col,
                         maintitle="", legendtitle="", legendtitle2="", flipPC1=F, flipPC2=F, discrete=T) {
  p <- pcoa_plot(dist, colors, fills, shapes, mysize, color_col, fill_col, shape_col,
                 maintitle, legendtitle, legendtitle2, flipPC1, flipPC2, discrete)
  p <- p + theme(title = element_text(face="bold", size=12),
                 axis.title = element_text(size=10),
                 legend.title = element_text(face="bold", size=10),
                 legend.text = element_text(size=9))
  return(p)
}



##### Generate Plots #####
# Soil Before
soil_col <- c("black", "#002dd5", "#1951f4", "#6e8df7", "#aabcf9", "black")
soil_fil <- c("black", "#002dd5", "#1951f4", "#6e8df7", "#aabcf9", "white")
soil_shp <- c(22, 25, 24, 21, 23, 23)
soil_jacc <- vegdist(t(soil_n), method="jaccard")
soil_before <- pcoa_plotf6(soil_jacc, soil_col, soil_fil, soil_shp, 3,
                         meta_soil$pH, meta_soil$pH, meta_soil$pH,
                         "88 Soils (Jaccard)", "pH level   ", flipPC1=T, flipPC2=T)

# Soil After (chose r = 0.865)
soil_lgd <- lg.dist(soil_jacc)
soil_after <- pcoa_plotf6(soil_lgd, soil_col, soil_fil, soil_shp, 3,
                        meta_soil$pH, meta_soil$pH, meta_soil$pH,
                        "88 Soils (LMdist-adjusted)", "pH level   ", flipPC1=T)

# Soil Smoothed
soil_slgd <- lg.dist(soil_jacc, smooth=T)
soil_smooth <- pcoa_plotf6(soil_slgd, soil_col, soil_fil, soil_shp, 3,
                         meta_soil$pH, meta_soil$pH, meta_soil$pH,
                         "88 Soils (LMdist-smoothed)", "pH level   ", flipPC1=T)


# Microbial mats Before
gn_col <- c("#F2182F", "#FF6A4F", "#FFA96A", "#FFDE95", "#FCF7C2",
            "#D6F0F6", "#88D9E8", "#23AED1", "#0076B4")
gn_jacc <- vegdist(t(gn_n), method="jaccard")
gn_before <- pcoa_plotf6(gn_jacc, gn_col, gn_col, rep(16,length(gn_col)), 4,
                       meta_gn$depth_mm, meta_gn$depth_mm, meta_gn$depth_mm,
                       "Microbial Mats (Jaccard)", "depth (mm)")

# Microbial mats After (chose r = 0.783)
gn_lgd <- lg.dist(gn_jacc)
gn_after <- pcoa_plotf6(gn_lgd, gn_col, gn_col, rep(16,length(gn_col)), 4,
                      meta_gn$depth_mm, meta_gn$depth_mm, meta_gn$depth_mm,
                      "Microbial Mats (LMdist-adjusted)", "depth (mm)")

# Microbial mats Smoothed
gn_slgd <- lg.dist(gn_jacc, smooth=T)
gn_smooth <- pcoa_plotf6(gn_slgd, gn_col, gn_col, rep(16,length(gn_col)), 4,
                       meta_gn$depth_mm, meta_gn$depth_mm, meta_gn$depth_mm,
                       "Microbial Mats (LMdist-smoothed)", "depth (mm)")


# Turkey Cecum Before
cecum_col <- c("#8e0013","#fa4d3f","#f99314","#00ae03","#57f89f","#0e2dfc","#6033f7") # day 1: "#fd2beb"
cecum_shp <- c(16, 17)
cecum_before <- pcoa_plotf6(cecum_d, cecum_col, cecum_col, cecum_shp, 3,
                          meta_cecum$Collection, meta_cecum$Collection, meta_cecum$System,
                          "Turkey Cecum (Bray-Curtis)", "Collection", "System")

# Turkey Cecum After (chose r = 0.618, using r = 0.7)
cecum_lgd <- lg.dist(cecum_d, 0.7)
cecum_after <- pcoa_plotf6(cecum_lgd, cecum_col, cecum_col, cecum_shp, 3,
                         meta_cecum$Collection, meta_cecum$Collection, meta_cecum$System,
                         "Turkey Cecum (LMdist-adjusted)", "Collection", "System")

# Turkey Cecum Smooth
cecum_slgd <- lg.dist(cecum_d, 0.7, smooth=T)
cecum_smooth <- pcoa_plotf6(cecum_slgd, cecum_col, cecum_col, cecum_shp, 3,
                          meta_cecum$Collection, meta_cecum$Collection, meta_cecum$System,
                          "Turkey Cecum (LMdist-smoothed)", "Collection", "System")


# Whittaker Before
whit_col <- viridis::plasma(20, alpha=0.8)
whit_shp <- c(16, 17, 15)
whit_before <- pcoa_plotf6(herb_d, whit_col, whit_col, whit_shp, 3,
                         locs$Elevation..m., locs$Elevation..m., locs$Site_Group,
                         "Whittaker Herbs (Bray-Curtis)", "Elevation", "Site",
                         discrete=F, flipPC2=T)

# Whittaker After (chose r = 0.675)
whit_lgd <- lg.dist(herb_d)
whit_after <- pcoa_plotf6(whit_lgd, whit_col, whit_col, whit_shp, 3,
                        locs$Elevation..m., locs$Elevation..m., locs$Site_Group,
                        "Whittaker Herbs (LMdist-adjusted)", "Elevation", "Site",
                        discrete=F)

# Whittaker Smooth
whit_slgd <- lg.dist(herb_d, smooth=T)
whit_smooth <- pcoa_plotf6(whit_slgd, whit_col, whit_col, whit_shp, 3,
                         locs$Elevation..m., locs$Elevation..m., locs$Site_Group,
                         "Whittaker Herbs (LMdist-smoothed)", "Elevation", "Site",
                         discrete=F)



##### GGPUBR: Organize plots #####
# FIGURE 6 : Plot Before and After with labels
fig6 <- ggarrange(soil_before, soil_after, gn_before, gn_after,
                  cecum_before, cecum_after, whit_before, whit_after,
                  nrow=2, ncol=4, labels=c("A","","B","","C","","D",""),
                  label.args=list(gp=grid::gpar(font=2, cex=1.5)))
## high res
tiff("figures/tif_files/figure6_case_studies.tif", width=19, height=7.2, unit="in", res=1200)
fig6
dev.off()
## low res
tiff("figures/tif_files_low_res/figure6_case_studies_lowres.tif", width=19, height=7.2, unit="in", res=350)
fig6
dev.off()


# SUPP FIGURE 2 : Smooth plots as well (exported as 1750 x 2000)
supp2 <- ggarrange(soil_before, soil_after, soil_smooth,
                   gn_before, gn_after, gn_smooth,
                   cecum_before, cecum_after, cecum_smooth,
                   whit_before, whit_after, whit_smooth,
                   nrow=4, ncol=3, labels=c("A","","","B","","","C","","","D","",""),
                   label.args=list(gp=grid::gpar(font=2, cex=1.5)))
## high res
tiff("figures/tif_files/suppfig2_case_studies_smoothing.tif", width=15, height=14, unit="in", res=1200)
supp
dev.off()
## low res
tiff("figures/tif_files_low_res/suppfig2_case_studies_smoothing_lowres.tif", width=15, height=14, unit="in", res=350)
supp2
dev.off()



##### STATS : Testing before vs. after (PERMANOVA) #####
# Soil Dataset
## PCoAs
soil_pc_before <- cmdscale(soil_jacc, k=2)
soil_pc_after <- cmdscale(soil_lgd, k=2)
## from original publication (ANOSIM, mantel test spearman)
anosim(soil_jacc, meta_soil$pH) # R: 0.56, p = 0.001
anosim(soil_lgd, meta_soil$pH)  # R: 0.51, p = 0.001
mantel(soil_jacc, vegdist(meta_soil$ph, "euclidean"), method="spearman") # r: 0.76, p = 0.001
mantel(soil_lgd, vegdist(meta_soil$ph, "euclidean"), method="spearman")  # r: 0.70, p = 0.001
## pH gradient: F-score increases! Significant before & after though
adonis2(soil_jacc ~ ph, data=meta_soil) # F score: 8.9, p = 0.001
adonis2(soil_lgd ~ ph, data=meta_soil)  # F score: 125.9, p = 0.001
cor.test(meta_soil$ph, soil_pc_before[,1]) # cor = 0.924, p < 2.2e-16
cor.test(meta_soil$ph, soil_pc_after[,1])  # cor = 0.929, p < 2.2e-16
## other variables (silt_clay, annual_season_temp, latitude/longitude are most correlated w/ PC2)
sort(abs(apply(meta_soil[,unlist(lapply(meta_soil, is.numeric))], 2, function (col) {cor(col, soil_pc_after[,2])})), decreasing=T)
##    silt clay
adonis2(soil_jacc ~ silt_clay, data=meta_soil)  # F score: 2.23, p = 0.001
adonis2(soil_lgd ~ silt_clay, data=meta_soil)   # F score: 5.3, p = 0.008
cor.test(meta_soil$silt_clay, soil_pc_before[,2]) # cor = -0.03, p = 0.781
cor.test(meta_soil$silt_clay, soil_pc_after[,2])  # cor = 0.60, p < 0.001
##    annual season temp
adonis2(soil_jacc ~ annual_season_temp, data=meta_soil)  # F score: 3.119, p = 0.001
adonis2(soil_lgd ~ annual_season_temp, data=meta_soil)   # F score: 9.507, p = 0.001
cor.test(meta_soil$annual_season_temp, soil_pc_before[,2]) # cor = 0.379, p = 0.0002
cor.test(meta_soil$annual_season_temp, soil_pc_after[,2])  # cor = 0.545, p < 0.0001


# Microbial Mats Dataset
## PCoAs
gn_pc_before <- cmdscale(gn_jacc, k=2)
gn_pc_after <- cmdscale(gn_lgd, k=2)
## depth gradient
adonis2(gn_jacc ~ end_depth, data=meta_gn)  # F score: 4.55, p = 0.001
adonis2(gn_lgd ~ end_depth, data=meta_gn)   # F score: 17.99, p = 0.002
cor.test(meta_gn$end_depth, gn_pc_before[,1], method="spearman") # rho = 0.941, p < 0.001
cor.test(meta_gn$end_depth, gn_pc_after[,1], method="spearman")  # cor = 0.987, p < 0.001
## other variables (arylsulfataseaandrelenzymes, sequencing)
sort(abs(apply(meta_gn[,unlist(lapply(meta_gn, is.numeric))], 2, function (col) {cor(col, gn_pc_after[,2])})), decreasing = T)
adonis2(gn_jacc ~ arylsulfataseaandrelenzymes, data=meta_gn)  # F score: 4.61, p = 0.001
adonis2(gn_lgd ~ arylsulfataseaandrelenzymes, data=meta_gn)   # F score: 12.05, p = 0.003
cor.test(meta_gn$arylsulfataseaandrelenzymes, gn_pc_before[,2]) # cor = 0.512, p = 0.029
cor.test(meta_gn$arylsulfataseaandrelenzymes, gn_pc_after[,2])  # cor = 0.610, p = 0.007
## linear regression with depth & PC1
tmp_gn <- meta_gn
gn_lm_bef <- lm(gn_pc_before[,1] ~ depth, data=tmp_gn) # first depth is not considered significant
summary(gn_lm_bef)
gn_lm_aft <- lm(gn_pc_after[,1] ~ depth, data=tmp_gn) # all depths significant
summary(gn_lm_aft)


# Turkey Cecum Dataset
## PCoAs
cecum_pc_before <- cmdscale(cecum_d, k=2)
cecum_pc_after <- cmdscale(cecum_lgd, k=2)
## time gradient
coll_day <- as.numeric(gsub("Day ","", meta_cecum$Collection))
adonis2(cecum_d ~ coll_day)   # F score: 23.78, p = 0.001
adonis2(cecum_lgd ~ coll_day) # F score: 47.61, p = 0.001
cor.test(coll_day, cecum_pc_before[,1]) # cor = 0.704, p < 2.2e-16
cor.test(coll_day, cecum_pc_after[,1])  # cor = 0.689, p < 2.2e-16
## "system" separation?
adonis2(cecum_d ~ coll_day + meta_cecum$System) # both signif
adonis2(cecum_lgd ~ coll_day + meta_cecum$System) # both signif
mantel(cecum_d, vegdist(coll_day, "euclidean"), method="pearson")   # r=0.38, p=0.001
mantel(cecum_lgd, vegdist(coll_day, "euclidean"), method="pearson") # r=0.37, p=0.001


# Whittaker Dataset
## PCoAs
whit_pc_before <- cmdscale(herb_d, k=2)
whit_pc_after <- cmdscale(whit_lgd, k=2)
## elevation gradient
adonis2(herb_d ~ Elevation..m., data=locs)   # F score: 36.8, p = 0.001
adonis2(whit_lgd ~ Elevation..m., data=locs) # F score: 166.5, p = 0.001
cor.test(locs$Elevation..m., whit_pc_before[,1]) # cor = 0.814, p < 2.2e-16
cor.test(locs$Elevation..m., whit_pc_after[,1])  # cor = 0.858, p < 2.2e-16
## site description, aspect
adonis2(herb_d ~ Elevation..m. + Site_Group + Aspect, data=locs)   # Site_Group: F score: 4.93, p = 0.001
adonis2(whit_lgd ~ Elevation..m. + Site_Group + Aspect, data=locs) # Site_Group: F score: 6.76, p = 0.001
kruskal.test(whit_pc_before[,2], factor(locs$Site_Group)) # chi-sq: 33.3, p < 0.001
kruskal.test(whit_pc_after[,2], factor(locs$Site_Group))  # chi-sq: 49.1, p < 0.001
## linear tests



