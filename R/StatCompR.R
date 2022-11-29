#' @title data beaver1
#' @name beaver1
#' @description The dataset containing beaver1, iris, WorldPhones
NULL

#' @title data iris
#' @name iris
#' @description The dataset containing beaver1, iris, WorldPhones
NULL

#' @title data WorldPhones
#' @name WorldPhones
#' @description The dataset containing beaver1, iris, WorldPhones
NULL


#' @title Calculate the Euclidean distance between two vectors
#'
#' @param x A vector
#' @param y A vector
#'
#' @return the Euclidean distance between two array of vectors
#' 
#' @examples
#' \dontrun{
#' x <- c(1,2,0)
#' y <- c(2,1,1)
#' euclidean(x,y)
#' }
#' @export
euclidean <- function(x, y){
  sqrt(sum((x-y)^2))
}

#' @title Calculate the pairwise Euclidean distance between two vector groups
#'
#' @param x An array of multiple vectors, which different rows are different vectors 
#' @param y An array of multiple vectors, which different rows are different vectors 
#'
#' @return Adjacency matrix between two sets of vectors
#' @examples
#' \dontrun{
#' x <- iris[,1:4][iris$Species == 'setosa',]
#' y <- iris[,1:4][iris$Species == 'versicolor',]
#' dis <- pairwise_distance(x,y)
#' }
#' @export
pairwise_distance <- function(x, y){
  dis = apply(x, 1, function(p) apply(y, 1, function(q) euclidean(p, q)))
  t(dis)
}


#' @title K-means cluster and K-medians cluster for a group of multiple vectors
#'
#' @param x An array of multiple vectors, which different rows are different vectors 
#' @param K Cluster numbers
#' @param center The initial center coordinates (default: NULL)
#' @param method str. "mean" (K-means) or "median" (K-medians) (default: "mean")
#' @param eps Conditions for Judging Convergence (default: 1e-5)
#' @param i_max Maximum number of iterations (default: 1000)
#'
#' @return A list containing clustering results of different vectors and cluster centers
#' @examples
#' \dontrun{
#' x <- iris[(iris$Species == 'setosa') | (iris$Species == 'versicolor'),]
#' result <- cluster.prototype(x[1:4], K=2)
#' x['cluster'] = result[[1]]
#' pca_iris = princomp(x[,1:4])$scores[,1:2]
#' plot(pca_iris, t='n')
#' text(pca_iris, labels=x$Species,col=x$cluster) 
#’ #text represents the annotation of the dataset, and the color is the clustering result, they are consistent
#' }
#'
#' @export
cluster.prototype <- function(x, K=3, center=NULL, method='mean', eps=1e-5, i_max=1000){
  if(is.null(center)) center = x[1:K, ]
  i=1
  i_max=i_max # The upper limit of the number of iterations
  new_center = center
  while(i < i_max){
    dis = pairwise_distance(x, center)
    cluster_index = apply(dis,1,order)[1,]
    
    # K-means method
    if(method == 'mean'){
      for(k in 1:K){
        if(sum(cluster_index==k)==1){
          new_center[k,] = x[cluster_index==k,]
          next
        }
        if(sum(cluster_index==k)==0) new_center[k,] = x[1,]
        new_center[k,] = colMeans(x[cluster_index==k,])
      }
    }
    
    # K-median method
    if(method == 'median'){
      for(k in 1:K){
        if(sum(cluster_index==k)==1){
          new_center[k,] = x[cluster_index==k,]
          next
        }
        if(sum(cluster_index==k)==0) new_center[k,] = x[1,]
        group_dis = pairwise_distance(x[cluster_index==k,],x[cluster_index==k,])
        center_index = order(rowSums(group_dis))[1]
        new_center[k,] = x[cluster_index==k,][center_index,]
      }
    }
    
    delta = abs(sum(sqrt(rowSums((new_center - center)^2))))
    cat('The MSE between',i-1,'and', i,'epochs = ',delta,'\n')
    if(delta < eps) {
      print('convenient')
      break
    }
    center = new_center
    i = i+1
  }
  return(list(cluster_index, new_center))
}


#' @title Calculate the distance between two clusters
#'
#' @param P An array of multiple vectors, which different rows are different vectors
#' @param Q An array of multiple vectors, which different rows are different vectors
#' @param method the method using to calculate the distance between two clusters. 
#' ('min', 'max', 'mean', default: 'min')
#'
#' @return The distance between two vector groups
#' @export
cluster_distanse <- function(P, Q, method='min'){
  dis <- pairwise_distance(P, Q)
  if(method == 'min')  return(min(dis))
  if(method == 'max')  return(max(dis))
  if(method == 'mean')  return(mean(dis))
}

#' @title Bottom-up Hierarchical Clustering using AGNES algorithm, 
#' including single-linkage, complete-linkage and average-linkage method
#'
#' @param x An array of multiple vectors, which different rows are different vectors
#' @param K Cluster numbers
#' @param method the method using to calculate the distance between two clusters. 
#' ('min', 'max', 'mean', default: 'min')
#'
#' @return A vector of the clusters id for each input vectors
#' @examples
#' \dontrun{
#' x <- iris[1:100,1:4]
#' label = cluster.prototype(x, K=2,method='mean')
#' }
#' @export
cluster.hierarchical <- function(x, K=3, method='min'){
  n <- dim(x)[1]
  C <- 1:n
  M <- pairwise_distance(x, x)
  for(i in 1:n){
    for(j in 1:i){
      M[i,j] = Inf
    }
  }
  k = n
  while(k > K){
    merge_ind <- c(which(M==min(M), arr.ind = TRUE)[1,])

    C[C==merge_ind[2]] <- merge_ind[1]
    C[C>merge_ind[2]] <- C[C>merge_ind[2]] - 1
    k = k-1
    M = M[-merge_ind[2], -merge_ind[2]]
    for(i in 1:k){ # update the distance matrix
      t = cluster_distanse(x[which(C==merge_ind[1]),], x[C==i,], method=method)
      M[merge_ind[1],i] = t
      M[i,merge_ind[1]] = t
    }
    for(i in 1:k){
      for(j in 1:i){
        M[i,j] = Inf
      }
    }
  }
  return(C)
}















