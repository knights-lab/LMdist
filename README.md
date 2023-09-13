# LMdist
LMdist is an unsupervised algorithm for correcting oversaturated pairwise distance and dissimilarity measures. Compatible with any pairwise beta diversity measure, LMdist will uncover underlying environmental gradients and resolve arches/horseshoes in ordination.



## Installation

As of now, the functions for using LMdist exist as R [^1] source code here on GitHub. To get started, download the source code file `lib/lmdist_source.r` and place this file in the desired location on your machine.

[^1]: The source code was made using R 4.3.1, other versions of R may be incompatible. To check your R version, enter `R.version` in the console.



## Usage

To use the relevant functions for LMdist, simply call the source file within your script using
```
source('/path/to/file/lmdist_source.r')
```

Below you will find descriptions of functions/parameters along with a more detailed tutorial.



## Introduction

Description of functions and parameters.

**lm.dist() :** Primary function for the LMdist algorithm, takes a distance object or matrix and returns an adjusted distance object of the same size, with values adjusted according to the local manifold. This is the only function a typical user will utilize, other functions included are helper functions called within `lm.dist()`.

| **Parameter** | **Description** |
| ---------- | ---------- |
| d | [required] Distance object or matrix to be adjusted |
| neighborhood.radius | [optional, default: NULL] Radius value to determine neighborhoods (NULL if algorithm should pick) |
| weighted | [optional, default: True] True/False use weighted edges in graph |
| epsilon | [optional, default: 0.05] Amount by which a smaller radius must be better correlated with the PCoA distances |
| phi | [optional, default: 0.10] Minimum graph degree:n ratio for which a radius is considered valid. |

**lm.evaluate() :** Helper function called within `lm.dist()` which returns a graph & relevant information for a particular neighborhood radius.

**lm.graph() :** Helper function called within `lm.evaluate()` to generate the graph of nodes and edges using a particular neighborhood radius.

**greedy.connect() :** Helper function called within `lm.evaluate()` to greedily connect disconnected graph components using the minimum distance between these components, if needed.

**lm.smooth() :** Helper function called within `lm.dist()` to optionally smooth the results of multiple neighborhood radii.



## Tutorial

### Dune Tutorial

Using the publicly available dune dataset [^2] from the vegan [^3] package, we demonstrate how the LMdist algorithm adjusts distances to more accurately depict sample relationships.

```r
source("lmdist_source.r")
library(vegan)
set.seed(25)

# Loading datasets
data(dune)
data(dune.env)

# Compute Bray-Curtis distances for the dune dataset, then visualize with PCoA
dune.d <- vegdist(dune, method="bray")
dune.pc <- cmdscale(dune.d, k=2, eig=F)
plot(dune.pc, pch=16, cex=2, col=c("#1E7D7D","#319DDC","#E4AE54","#F5674E")[dune.env$Moisture], xlab="PC 1", ylab="PC 2", main="Original PCoA (dune)")
legend("topright", pch=16, col=c("#1E7D7D","#319DDC","#E4AE54","#F5674E"), legend=levels(dune.env$Moisture), title="Moisture")
```

```r
# Adjust the distances using LMdist, then visualize the adjusted distances with PCoA.
dune.lmd <- lm.dist(dune.d)
dune.pc.lmd <- cmdscale(dune.lmd, k=2, eig=F)
plot(dune.pc.lmd, pch=16, cex=2, col=c("#1E7D7D","#319DDC","#E4AE54","#F5674E")[dune.env$Moisture], xlab="PC 1", ylab="PC 2", main="LMdist PCoA (dune, defaults)")
legend("bottomleft", pch=16, col=c("#1E7D7D","#319DDC","#E4AE54","#F5674E"), legend=levels(dune.env$Moisture), title="Moisture")
```

```r
# Use the "smooth" option to get combined results from multiple radii, that way not just one radius value is used.
dune.lmd <- lm.dist(dune.d, smooth=T)
dune.pc.lmd <- cmdscale(dune.lmd, k=2, eig=F)
plot(dune.pc.lmd, pch=16, cex=2, col=c("#1E7D7D","#319DDC","#E4AE54","#F5674E")[dune.env$Moisture], xlab="PC 1", ylab="PC 2", main="LMdist PCoA (dune, smoothed)")
legend("bottomleft", pch=16, col=c("#1E7D7D","#319DDC","#E4AE54","#F5674E"), legend=levels(dune.env$Moisture), title="Moisture")
```


### Iris tutorial

Using the public iris dataset [^4], we see that LMdist does not adjust distances if pairwise distances are not necessarily oversaturated. However, we can override the default LMdist algorithm to adjust distances using a provided neighborhood radius value.

```r
source("lmdist_source.r")
set.seed(25)

# Loading the iris dataset
data(iris)

# Compute Euclidean distances of the iris dataset and visualize these sample relationships using PCA.
iris.d <- dist(iris[,1:4])
iris.pc <- cmdscale(iris.d, k=2, eig=F)
plot(iris.pc, pch=16, cex=1.5, col=c("purple","orange","blue")[factor(iris$Species)], xlab="PC 1", ylab="PC 2", main="Original PCA (iris)")
legend("bottomright", pch=16, col=c("purple","orange","blue"), legend=levels(factor(iris$Species)), title="Species")

# The default LMdist algorithm tries 50 neighborhood radii and returns the best fit for the data.
# In the iris dataset, you'll note the algorithm does not adjust distances, because the distances are not oversaturated.
iris.lmd <- lm.dist(iris.d)
iris.pc.lmd <- cmdscale(iris.lmd, k=2, eig=F)
plot(iris.pc.lmd, pch=16, cex=1.5, col=c("purple","orange","blue")[factor(iris$Species)], xlab="PC 1", ylab="PC 2", main="LMdist PCA (iris, defaults)")
legend("bottomright", pch=16, col=c("purple","orange","blue"), legend=levels(factor(iris$Species)), title="Species")


# We can force the LMdist algorithm to adjust the pairwise distances by including a specific radius value to be used.
iris.lmd <- lm.dist(iris.d, 2)
iris.pc.lmd <- cmdscale(iris.lmd, k=2, eig=F)
plot(iris.pc.lmd, pch=16, cex=1.5, col=c("purple","orange","blue")[factor(iris$Species)], xlab="PC 1", ylab="PC 2", main="LMdist PCA (iris, radius 2)")
legend("bottomleft", pch=16, col=c("purple","orange","blue"), legend=levels(factor(iris$Species)), title="Species")
```


[^2]: Batterink M. & Wijffels G. (1983): Een vergelijkend vegetatiekundig onderzoek naar de typologie en invloeden van het beheer van 1973 tot 1982 in de duinweilanden op Terschelling. Report Agricultural University, Department of Vegetation Science, Plant Ecology and Weed Science, Wageningen.

[^3]: Oksanen J, Simpson G, Blanchet F, Kindt R, Legendre P, Minchin P, O'Hara R, Solymos P, Stevens M, Szoecs E, Wagner H, Barbour M, Bedward M, Bolker B, Borcard D, Carvalho G, Chirico M, De Caceres M, Durand S, Evangelista H, FitzJohn R, Friendly M, Furneaux B, Hannigan G, Hill M, Lahti L, McGlinn D, Ouellette M, Ribeiro Cunha E, Smith T, Stier A, Ter Braak C, Weedon J (2022). *vegan: Community Ecology Package*. R package version 2.6-4, [https://CRAN.R-project.org/package=vegan](https://CRAN.R-project.org/package=vegan).

[^4] Anderson, Edgar (1935). The irises of the Gaspe Peninsula, Bulletin of the American Iris Society, 59, 2–5.
