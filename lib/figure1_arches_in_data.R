# Figure 1 : Example arch in data (soil)
library(vegan)
source("lib/figurehelper_pcoaplot.R")
set.seed(125)

# Load soil dataset (normalized data)
## normalized relative abundances
soil_n <- read.delim("data/soil/44766_clean_otus_norm.txt", row=1, header=T, sep="\t")
## metadata, reformatted pH column
meta_soil <- read.table("data/soil/clean_map.txt", header=T, sep="\t")
rownames(meta_soil) <- meta_soil$SampleID
meta_soil$pH <- cut(meta_soil$ph, breaks = c(0,4:8,14))
levels(meta_soil$pH) <- c("<4","4-5","5-6","6-7","7-8",">8")

# Generate soil PCoA plot
soil_col <- c("black", "#002dd5", "#1951f4", "#6e8df7", "#aabcf9", "black")
soil_fil <- c("black", "#002dd5", "#1951f4", "#6e8df7", "#aabcf9", "white")
soil_shp <- c(22, 25, 24, 21, 23, 23)
soil_jacc <- vegdist(t(soil_n), method="jaccard")
figure1 <- pcoa_plot(soil_jacc, soil_col, soil_fil, soil_shp, 3,
                     meta_soil$pH, meta_soil$pH, meta_soil$pH,
                     "88 Soils (Jaccard)", "pH level", flipPC1=F, flipPC2=T)

# Save a full size tif (350 dpi)
tiff("figures/tif_files/figure1_soil_pcoa.tif", width=4, height=4, units="in", res=1200)
figure1
dev.off()

# Lower resolution for initial submission
tiff("figures/tif_files_low_res/figure1_soil_pcoa_lowres.tif", width=4, height=4, units="in", res=400)
figure1
dev.off()
