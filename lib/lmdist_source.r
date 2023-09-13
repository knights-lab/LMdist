# lmdist_source.r
# usage : source('lmdist_source.r')


##### Functions #####
"lm.dist" <- function(d, neighborhood.radius=NULL, weighted=TRUE, smooth=FALSE, epsilon=0.05, phi=0.10) {
  # This master function calculates distances based on connected graph of neighborhoods
  # d : distance object to be adjusted
  # neighborhood.radius : radius to determine neighborhoods (NULL if algorithm should pick)
  # weighted : True/False use weighted edges in graph (recommended: TRUE)
  # smooth : True/False smooth results of multiple radii, specifically all valid radii in neighborhood.radius (default: False)
  # epsilon : (default: 0.05)
  # phi : Percentage of samples which must fall within (default: 0.10)
  
  # Verify parameter validity
  ## distance object or distance matrix must be provided
  if (class(d) != "dist") {
    if (isSymmetric(as.matrix(d)) == TRUE) {
      d <- as.dist(d)
    } else {
      stop("'d' must be a distance object or symmetric distance matrix.")
    }
  }
  if (dim(as.matrix(d))[1] <= 10 ) {
    stop("'d' is too small, recommended to use a larger dataset, minimum of 10 samples.")
  }
  ## neighborhood radius must be a double or list of doubles or NULL (default)
  if (typeof(neighborhood.radius) != "double" & suppressWarnings(any(is.na(as.double(neighborhood.radius)) == TRUE))) {
    stop("'neighborhood.radius' must be numeric type 'double' or 'NULL' (default) to automatically find best neighborhood.")
  } else if (is.null(neighborhood.radius)) {
    neighborhood.sizes = seq(max(d), min(d), length.out=50) # try 50 radii in range of d
  } else {
    neighborhood.sizes = sort(as.double(neighborhood.radius), decreasing=T)
  }
  ## weighted should be True (default) or False
  if (typeof(weighted) != "logical") {
    stop("'weighted' must be logical type TRUE (default) or FALSE.")
  }
  ## make sure epsilon does not make it impossible to improve
  tmp <- cmdscale(d, k=5)
  tmp <- round(cor(dist(tmp), as.dist(d), method="pearson"),3)
  if (((epsilon + tmp) > 0.995) & tmp < 0.99 ) {
    warning(paste0("Epsilon too large for dataset, updated to ", round((1-tmp)/2, 3)))
    epsilon <- round((1-tmp)/2, 3) # smaller epsilon (for smaller datasets with high correlation but arch present)
  }

  #  Create connected graph with provided neighborhood range(s)
  ix <- 1         # index
  lmd <- NULL     # distances
  curr_corr <- -1 # optimization function score
  curr_ratio <- 1 # ratio degree:n
  curr_radius <- neighborhood.sizes[1]
  mylist <- list()
  while(ix <= length(neighborhood.sizes)){
    ns <- neighborhood.sizes[ix]
    out <- lm.evaluate(d, neighborhood.size=ns, weighted=weighted) # list(lmd, is.valid, ratio, corr)
    if (!out[[2]]) {
      stop(paste0('Graph is invalid (r = ',ns,'), check all parameters are valid or try a different radius.'))
    }
    if (out[[3]] < phi) break
    if ((out[[3]] >= phi) & (out[[4]] > (curr_corr + epsilon))) {
      if (ix == 1) {print(paste0("setting initial r as best --> ", round(ns,3)))}
      else {print(paste0("new best r found --> ", round(ns,3)))}
      lmd <- out[[1]]
      curr_ratio <- out[[3]]
      curr_corr <- out[[4]]
      curr_radius <- ns
    }
    mylist <- append(mylist, list(out[[1]]))
    ix <- ix + 1
  }
  # Smooth results if desired
  if (smooth == TRUE) {
    if (length(mylist) > 1) {
      lmd <- lm.smooth(mylist, which(neighborhood.sizes == curr_radius), neighborhood.sizes[1:length(mylist)])
    } else {
      # adding up to 12 on each side of the one provided neighborhood
      tmp_rvals <- round(seq(max(d), curr_radius, length.out=12), 3)
      tmp_rvals <- c(tmp_rvals, round(seq(curr_radius, min(d), length.out=12), 3))
      val_rvals <- c(curr_radius)
      for (r in tmp_rvals) {
        tmp_out <- lm.evaluate(d, neighborhood.size=r, weighted=weighted)
        if (tmp_out[[3]] >= phi) {
          mylist <- append(mylist, list(tmp_out[[1]]))
          val_rvals <- c(val_rvals, r)
        }
      }
      lmd <- lm.smooth(mylist, 1, val_rvals)
    }
  }

  # Return dist object
  if (is.null(lmd)) {stop(paste0("Radius ", round(ns, 3)," resulted in null graph. Try larger radius.")); return(NULL);}
  return(as.dist(lmd))
}

"lm.evaluate" <- function(d, neighborhood.size=NULL, weighted=TRUE) {
  # For a given radius value, returns graph & relevant information
  is.valid <- FALSE         # default, must be proven otherwise
  n <- dim(as.matrix(d))[1] # number of samples
  g <- lm.graph(d, neighborhood.size=neighborhood.size, weighted=weighted)
  if (clusters(g)$no == 1) {
    is.valid <- TRUE #eigs <- eigen(graph.laplacian(g), only.values=TRUE)$values
  } else {
    # if invalid graph, use minimum spanning tree to add minimum number of edges
    g <- greedy.connect(d, g, neighborhood.size)
    if (clusters(g)$no == 1) {
      is.valid <- TRUE
    }
  }
  ratio <- mean(degree(g)) / n # record degree:n ratio
  if(is.valid) {
    # Compute the LMDist-adjusted distances
    lm_d <- shortest.paths(g)
    # Determine the best number of dimensions of PCoA with which to compare, the fewest number of dimensions (between 2 and 10) to describe >80% variance
    pc <- suppressWarnings(cmdscale(lm_d, k=10, eig=T))
    scree <- (pc$eig[1:10] / sum(pc$eig))
    scree_tot <- c()
    for (i in 1:length(scree)) { scree_tot <- c(scree_tot, sum(scree[1:i])) }
    ndim <- min(which(scree_tot >= 0.8))
    if (ndim == 1) { ndim <- 2 }
    corr <- round(cor(dist(pc$points[,1:ndim]), as.dist(lm_d), method="pearson"),3) # use correlation b/w PC & lmd dists as optimization func
  } else {
    lm_d <- NULL
    corr <- -1
  }
  ## Debugging feature : Print igraph
  # png(paste0("results/igraphs_gif/r_",round(neighborhood.radius,4),".png"), width=600, height=600)
  # plot(g, vertex.label=NA, vertex.color="purple", vertex.size=8,
  #      main=paste0("r = ", round(neighborhood.radius,4)))
  # dev.off()
  
  return(list(lm_d, is.valid, ratio, corr))
}

"lm.graph" <- function(d, neighborhood.size=NULL, weighted=TRUE) {
  # Creates graph with provided neighborhood metric
  require('igraph', warn.conflicts=FALSE, quietly=TRUE)
  d <- as.matrix(d)
  adj <- matrix(0, nrow(d), nrow(d))
  for(i in 1:nrow(adj)) {
    # retired neighborhood size n: adj[i,order(d[i,])[1:(neighborhood.size+1)]] <- 1
    adj[i,which(d[i,] <= neighborhood.size)] <- 1
  }
  if(weighted) adj[adj>0] <- d[adj>0] else weighted <- NULL
  rownames(adj) <- rownames(d); colnames(adj) <- colnames(d);
  g <- graph.adjacency(adj, weighted=weighted, mode='undirected')
  return(g)
}

"greedy.connect" <- function (d, g, ns) {
  # Greedily connects disconnected components with minimum distance between components
  # note: prints warning to output to notify user of this operation
  # retired method : dc <- split(1:nrow(as.matrix(d)), components(g)$membership)
  og <- graph.adjacency(as.matrix(d), weighted=T, mode='undirected')
  mstree <- mst(og) # retired: graph.adjacency(mst(d), weighted=TRUE, mode='undirected')
  me <- split(as_edgelist(mstree), seq(gsize(mstree))); ge <- split(as_edgelist(g), seq(gsize(g)));
  new_edges <- unname(unlist(me[!(me %in% ge)]))
  d_idx <- sapply(seq(1,length(new_edges),2), function (i) {
    tmpd <- as.matrix(d)
    a <- which(rownames(tmpd) == new_edges[i])    # get row idx
    b <- which(colnames(tmpd) == new_edges[i+1])  # get column index
    return(c(a, b))
  })
  w <- sapply(1:ncol(d_idx), function (i) {as.matrix(d)[d_idx[1,i], d_idx[2,i]]})
  warning(paste0("Disconnected graph with given neighborhood size (", round(ns,3), "). Adding ", length(new_edges)/2, " edges."), call.=FALSE)
  if (any(w > 2*max(E(g)$weight))) warning("Some added edge weights exceed twice length of neighborhood edges.")
  g <- add_edges(g, new_edges, attr = list("weight"=w))
  return(g)
}

"lm.smooth" <- function (lmdlist, bestidx, rvals) {
  require('rdd', warn.conflicts=FALSE, quietly=TRUE)
  if (is.null(lmdlist)) warning("Smoothing not possible, null list.")
  # Determine new distances as weighted mean of other r distance outputs
  ## (1) get weights centered at best r index
  ##     note: using a Gaussian with bandwidth ~15% of length
  wts <- kernelwts(1:length(rvals), center=bestidx, bw=ceiling(length(rvals)*0.15), kernel="gaussian")
  ## (2) per pairwise distance, get weighted mean
  means <- lmdlist[[1]]    # temporarily set means to first distance object
  for (i in 1:length(lmdlist[[1]])) {
    x <- sapply(lmdlist, function (l) { l[[i]] })
    means[i] <- weighted.mean(x, wts)
  }
  ## (3) return in a distance object
  return(as.dist(means))
}
