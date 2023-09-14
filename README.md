# LMdist
LMdist is an unsupervised algorithm for correcting oversaturated pairwise distance and dissimilarity measures. Compatible with any pairwise beta diversity measure, LMdist will uncover underlying environmental gradients and resolve arches/horseshoes in ordination.


---


## Installation

As of now, the functions for using LMdist exist as R [^1] source code here on GitHub. To get started, download the source code file `lib/lmdist_source.r` and place this file in the desired location on your machine.

[^1]: The source code was made using R 4.3.1, other versions of R may be incompatible. To check your R version, enter `R.version` in the console.


---


## Usage

To use the relevant functions for LMdist, simply call the source file within your script using
```r
source('/path/to/file/lmdist_source.r')
```

Then you can call the LMdist function using
```r
lm.dist(distance_object)
```
Continue below for descriptions of functions/parameters along with a more detailed tutorial.


---


## Introduction

Description of functions and parameters.

`lm.dist()` : Primary function for the LMdist algorithm, takes a distance object or matrix and returns an adjusted distance object of the same size, with values adjusted according to the local manifold. This is the only function a typical user will utilize, other functions included are helper functions called within `lm.dist()`.

| **Parameter** | **Description** |
| ---------- | ---------- |
| **d** | [required] Distance object or matrix to be adjusted |
| **neighborhood.radius** | [optional, default: NULL] Radius value to determine neighborhoods (NULL if algorithm should pick) |
| weighted | [optional, default: True] True/False use weighted edges in graph |
| epsilon | [optional, default: 0.05] Amount by which a smaller radius must be better correlated with the PCoA distances |
| phi | [optional, default: 0.10] Minimum graph degree:n ratio for which a radius is considered valid. |
| **smooth** | [optional, default: F] Smooth results of multiple radii around the chosen "best" radius value. Best radius value is set to neighborhood.radius if a single value is provided, and 25 radii around this value are added to create a weighted average for smoothed results. |

`lm.evaluate()` : Helper function called within `lm.dist()` which returns a graph & relevant information for a particular neighborhood radius.

`lm.graph()` : Helper function called within `lm.evaluate()` to generate the graph of nodes and edges using a particular neighborhood radius.

`greedy.connect()` : Helper function called within `lm.evaluate()` to greedily connect disconnected graph components using the minimum distance between these components, if needed.

`lm.smooth()` : Helper function called within `lm.dist()` to optionally smooth the results of multiple neighborhood radii.


---


## Tutorials

### Dune Tutorial

Using the publicly available dune dataset [^2] from the vegan [^3] package, we demonstrate how the LMdist algorithm adjusts distances to more accurately depict sample relationships.

We can begin by sourcing the LMdist source code and loading the dune dataset.

```r
source("lmdist_source.r")
library(vegan)
set.seed(25)

# Loading datasets
data(dune)
data(dune.env)
```

Next, we compute Bray-Curtis distances, visualizing the outcome using PCoA.

```r
# Compute Bray-Curtis distances for the dune dataset, then visualize with PCoA
dune.d <- vegdist(dune, method="bray")
dune.pc <- cmdscale(dune.d, k=2, eig=F)
plot(dune.pc, pch=16, cex=2, col=c("#1E7D7D","#319DDC","#E4AE54","#F5674E")[dune.env$Moisture], xlab="PC 1", ylab="PC 2", main="Original PCoA (dune)")
legend("topright", pch=16, col=c("#1E7D7D","#319DDC","#E4AE54","#F5674E"), legend=levels(dune.env$Moisture), title="Moisture")
```
![PCoA of dune dataset with Bray-Curtis dissimilarity. Samples arranged in a sparse arch formation.](tutorial_figs/dune_original_pcoa.png "PCoA (dune data)")


We can see that the dataset may be presenting an arch-like geometry, largely driven by the moisture gradient. Since LMdist works best in correcting oversaturated distances, we can look at the distribution of distances to see if there is a bounded effect here.

```r
# Distribution of distances
plot(hist(dune.d, breaks=12), main="Distribution of Dissimilarities (dune)", xlab="Bray-Curtis Dissimilarity")
```
![Distribution of Bray-Curtis dissimilarities present a slight left skew.](tutorial_figs/dune_original_distribution.png "Distribution of Bray-Curtis dissimilarities (dune data)")

There appears to be a very slight left skew in this dataset, so we can try using LMdist to adjust distances. We are using the default algorithm which will try 50 radii values for adjusting distances and pick the one which represents true pairwise distances the best.

```r
# Adjust the distances using LMdist
dune.lmd <- lm.dist(dune.d)
```
```
    [1] "setting initial r as best --> 1"
    [1] "new best r found --> 0.543"
```

We can see the default algorithm chooses to adjust distances, picking a best radius of 0.543 for creating the underlying manifold graph. Let's visualize these adjusted distances with PCoA.

```r
# Visualize the LMdist adjusted values with PCoA
dune.pc.lmd <- cmdscale(dune.lmd, k=2, eig=F)
plot(dune.pc.lmd, pch=16, cex=2, col=c("#1E7D7D","#319DDC","#E4AE54","#F5674E")[dune.env$Moisture], xlab="PC 1", ylab="PC 2", main="LMdist PCoA (dune, defaults)")
legend("bottomleft", pch=16, col=c("#1E7D7D","#319DDC","#E4AE54","#F5674E"), legend=levels(dune.env$Moisture), title="Moisture")

# Distribution of LMdist-adjusted distances
plot(hist(dune.lmd, breaks=12), main="Distribution of LMdist-adjustment (dune)", xlab="LMdist-adjusted (r=0.543)")
```

![PCoA of dune dataset with Lmdist-adjustment. Samples are more linearly arranged.](tutorial_figs/dune_lmdist_pcoa.png "LMdist PCoA (dune data)")

The adjusted plot appears to have flattened the curve from the original PCoA along the x-axis (PC 1). Other variation in the samples is now represented by the second axis (PC 2), notably the variation in low moisture samples.

![LMdist distribution of dissimilarities in the dune dataset reduces left skew.](tutorial_figs/dune_lmdist_distribution.png "LMdist distribution of dissimilarities (dune data)")

The distribution of distances has reduced left skew after adjustment with LMdist.


While a simple solution, it may be difficult to trust one radius value alone. Researchers may therefore choose to apply an optional smoothing to adjust distances according to more than one radius value by a Gaussian weighting.

```r
# Use the "smooth" option to get combined results from multiple radii, thereby not placing complete trust in just one radius.
dune.lmd <- lm.dist(dune.d, smooth=T)
```
```
    [1] "setting initial r as best --> 1"
    [1] "new best r found --> 0.543"
```

The same best radius is picked because the algorithm remains the same in optimizing the best radius value. Now, though, a Gaussian weight centered at the best radius value is applied to the results of multiple radii around this radius, returning a weighted average.

```r
# Visualize results with smoothing applied
dune.pc.lmd <- cmdscale(dune.lmd, k=2, eig=F)
plot(dune.pc.lmd, pch=16, cex=2, col=c("#1E7D7D","#319DDC","#E4AE54","#F5674E")[dune.env$Moisture], xlab="PC 1", ylab="PC 2", main="LMdist PCoA (dune, smoothed)")
legend("bottomleft", pch=16, col=c("#1E7D7D","#319DDC","#E4AE54","#F5674E"), legend=levels(dune.env$Moisture), title="Moisture")

# Distribution of LMdist-adjusted distances (with smoothing)
plot(hist(dune.lmd, breaks=12), main="Distribution of LMdist-adjustment (dune, smooth)", xlab="LMdist-adjusted (r=0.543, smooth)")
```

![PCoA of dune dataset with Lmdist-adjustment and smoothing. Samples are more linearly arranged.](tutorial_figs/dune_lmdist_pcoa.png "LMdist smooth PCoA (dune data)")

The results are similar to the single best radius, but you may not small adjustments for individual samples, an effect induced by this weighted averaging of other radii.

![LMdist distribution of dissimilarities in the dune dataset reduces left skew.](tutorial_figs/dune_lmdist_distribution.png "LMdist smooth distribution of dissimilarities (dune data)")

The distribution of distances also resembles the default algorithm.

Smoothing can be applied to any LMdist output, including a user-provided radius value. The next tutorial with the `iris` dataset walks through applying a user-defined radius.


---


### Iris tutorial

Using the public iris dataset [^4], we see that LMdist does not adjust distances if pairwise distances are not oversaturated. However, we can override the default LMdist algorithm to adjust distances using a provided neighborhood radius value.

```r
source("lmdist_source.r")
set.seed(25)

# Loading the iris dataset
data(iris)
```

Visualizing Euclidean distances of the iris dataset, we can see samples cluster according to 3 species groups.

```r
# Compute Euclidean distances of the iris dataset and visualize these sample relationships using PCA.
iris.d <- dist(iris[,1:4])
iris.pc <- cmdscale(iris.d, k=2, eig=F)
plot(iris.pc, pch=16, cex=1.5, col=c("purple","orange","blue")[factor(iris$Species)], xlab="PC 1", ylab="PC 2", main="Original PCA (iris)")
legend("bottomright", pch=16, col=c("purple","orange","blue"), legend=levels(factor(iris$Species)), title="Species")
```

![PCoA of iris dataset with Euclidean distances forms three clusters by species.](tutorial_figs/iris_original_pcoa.png "PCoA (iris data)")

Theses distances are also not oversaturated, no left skew is present.

```r
# Density of distances
plot(hist(iris.d, breaks=20), main="Distribution of Distances (iris)", xlab="Euclidean Distances")
```

![Distribution of Euclidean distances in iris dataset looks normally distributed.](tutorial_figs/iris_original_distribution.png "Distribution of Euclidean distances (iris data)")

It is therefore unsurprising that the LMdist default algorithm does not choose to adjust the distances, using a radius value of the maximum distance and therefore trusting all provided distances.

```r
# Apply the default LMdist function
## Note that no adjustment is made because the original pairwise distances are already best representing the sample relationships.
iris.lmd <- lm.dist(iris.d)
```

```
    [1] "setting initial r as best --> 7.085"
```

```
# Visualizing the output with PCA
iris.pc.lmd <- cmdscale(iris.lmd, k=2, eig=F)
plot(iris.pc.lmd, pch=16, cex=1.5, col=c("purple","orange","blue")[factor(iris$Species)], xlab="PC 1", ylab="PC 2", main="LMdist PCA (iris, defaults)")
legend("bottomright", pch=16, col=c("purple","orange","blue"), legend=levels(factor(iris$Species)), title="Species")
```

![PCoA of iris dataset with LMdist-adjusted Euclidean distances is unchanged from original.](tutorial_figs/iris_lmdist_defaults_pcoa.png "LMdist PCoA (iris data, defaults)")

So, the pairwise relationships and subsequent PCA plot are unchanged. 

However, if we want to explore the dataset further by applying a particular radius value, we can do so with the "neighborhood.radius" parameter. Here, we apply a radius value of 2.

```r
# We can force the LMdist algorithm to adjust the pairwise distances by including a specific radius value (or values) to be used.
iris.lmd <- lm.dist(iris.d, 2)
```

```
    [1] "setting initial r as best --> 2"
```

```r
# Visualize these adjusted values
iris.pc.lmd <- cmdscale(iris.lmd, k=2, eig=F)
plot(iris.pc.lmd, pch=16, cex=1.5, col=c("purple","orange","blue")[factor(iris$Species)], xlab="PC 1", ylab="PC 2", main="LMdist PCA (iris, radius 2)")
legend("bottomleft", pch=16, col=c("purple","orange","blue"), legend=levels(factor(iris$Species)), title="Species")

# Visualize the distribution of adjusted values
plot(hist(iris.lmd, breaks=20), main="Distribution of LMdist (iris, r=2)", xlab="LMdist adjusted (r=2)")
```

![PCoA of iris dataset with LMdist-adjusted Euclidean distances using radius 2 causes the setosa cluster to become tighter.](tutorial_figs/iris_lmdist_radius2_pcoa.png "LMdist PCoA (iris data, r=2)")

The smaller radius seems to accentuate sample distances within the *Iris versicolor* and *Iris virginica* species, making the *Iris setosa* species appear very tightly clustered by comparison in PC1 and PC2.

![Distribution of LMdist-adjusted Euclidean distances in iris dataset is similar to original.](tutorial_figs/iris_lmdist_radius2_distribution.png "Distribution of LMdist-adjusted Euclidean distances (iris, r=2)")

The distribution of distances is largely unchanged because the distinctive species groups persisted.

This example with the iris dataset shows how the LMdist algorithm may be used for exploration but is also an example of when the LMdist algorithm is not needed. The lack of oversaturation of distances means Euclidean distances are already representing sample relationships well (as expected from a dataset with only 4 initial dimensions).


---


[^2]: Batterink M. & Wijffels G. (1983): Een vergelijkend vegetatiekundig onderzoek naar de typologie en invloeden van het beheer van 1973 tot 1982 in de duinweilanden op Terschelling. Report Agricultural University, Department of Vegetation Science, Plant Ecology and Weed Science, Wageningen.

[^3]: Oksanen J, Simpson G, Blanchet F, Kindt R, Legendre P, Minchin P, O'Hara R, Solymos P, Stevens M, Szoecs E, Wagner H, Barbour M, Bedward M, Bolker B, Borcard D, Carvalho G, Chirico M, De Caceres M, Durand S, Evangelista H, FitzJohn R, Friendly M, Furneaux B, Hannigan G, Hill M, Lahti L, McGlinn D, Ouellette M, Ribeiro Cunha E, Smith T, Stier A, Ter Braak C, Weedon J (2022). *vegan: Community Ecology Package*. R package version 2.6-4, [https://CRAN.R-project.org/package=vegan](https://CRAN.R-project.org/package=vegan).

[^4] Anderson, Edgar (1935). The irises of the Gaspe Peninsula, Bulletin of the American Iris Society, 59, 2–5.
