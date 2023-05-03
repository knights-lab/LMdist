# Possible supplemental figure : coenoclines of the soil dataset
library(scales)
set.seed(125)

soil <- read.delim("data/soil/44766_clean_otus.txt", row=1, header=T)
soil_meta <- read.delim("data/soil/clean_map.txt", row=1, header=T)
identical(rownames(soil_meta), colnames(soil))
ph_gradient <- soil_meta$ph


# Reduce to most abundant OTUs
mostabun <- soil[head(order(rowSums(soil), decreasing=T),50),]


# Plotting a smoothed curve : What method to use?
# ord <- order(ph_gradient)
# lo <- loess(unlist(mostabun[1,])[ord] ~ ph_gradient[ord])
# ssp <- smooth.spline(x=ph_gradient[ord], y=unlist(mostabun[1,])[ord], spar=0.8)
# plot(x=ph_gradient[ord], y=unlist(mostabun[1,])[ord], type="l", lwd=2, main="Soil Coenoclines", xlab="pH gradient", ylab="OTU abundance")
# lines(predict(lo), col="blue", lwd=2)
# lines(ssp, col="green", lwd=2) # prefer this smoothing spline


# Set up colors
# PCoA colors (low to high pH): "black", "#002dd5", "#1951f4", "#6e8df7", "#aabcf9", "white"
mycolors <- alpha(viridis::mako(nrow(mostabun)),0.8)
avg_peakph <- apply(mostabun, 1, function (row) {mean(ph_gradient[tail(order(row),5)])})
mostabun <- mostabun[order(avg_peakph),] # re-order the data frame to be ordered by colors desired

# Plot the top 50 most abundant
ord <- order(ph_gradient)
spl <- smooth.spline(x=ph_gradient[ord], y=unlist(mostabun[1,])[ord], spar=0.8)
plot(spl, type="l", col=mycolors[1], lwd=3, ylim=c(0,50),
     main="Soil Coenoclines: Top 50 OTUs", xlab="pH gradient", ylab="OTU abundance (smoothed spline)")
for (i in 2:nrow(mostabun)) {
  spline <- smooth.spline(x=ph_gradient[ord], y=unlist(mostabun[i,])[ord], spar=0.8)
  lines(spline, col=mycolors[i], lwd=3)
}
supp1 <- recordPlot()

# Export to TIFF
## high res
tiff("figures/tif_files/suppfig1_soil_coenoclines.tif", width=5, height=4, unit="in", res=1200)
supp1
dev.off()
## low res
tiff("figures/tif_files_low_res/suppfig1_soil_coenoclines_lowres.tif", width=5, height=4, unit="in", res=350)
supp1
dev.off()



# Plot all OTUs
## ordering the dataframe for colors
# allcolors <- alpha(viridis::mako(nrow(soil)), 0.8)
# all_avg_peak <- apply(soil, 1, function (row) {mean(ph_gradient[tail(order(row),5)])})
# df <- soil[order(all_avg_peak),]
# ## plotting
# ord <- order(ph_gradient)
# spl <- smooth.spline(x=ph_gradient[ord], y=unlist(df[1,])[ord], spar=0.8)
# plot(spl, type="l", col=allcolors[1], lwd=3, ylim=c(0,50),
#      main="Soil Coenoclines: All OTUs", xlab="pH gradient", ylab="OTU abundance (smoothed spline)")
# for (i in 2:nrow(soil)) {
#   spline <- smooth.spline(x=ph_gradient[ord], y=unlist(df[i,])[ord], spar=0.8)
#   lines(spline, col=allcolors[i], lwd=3)
# }
