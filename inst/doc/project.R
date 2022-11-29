## -----------------------------------------------------------------------------
library(StatComp22023)
data(data)

## -----------------------------------------------------------------------------
x <- iris[1:100,]
result <- cluster.prototype(x[1:4], K=2)
x['cluster'] = result[[1]]

## -----------------------------------------------------------------------------
result[[2]]

## -----------------------------------------------------------------------------
pca_iris = princomp(x[,1:4])$scores[,1:2]
plot(pca_iris, t='n')
text(pca_iris, labels=x$Species,col=x$cluster) 

## -----------------------------------------------------------------------------
x <- iris[(iris$Species == 'setosa') | (iris$Species == 'versicolor'),]
result <- cluster.prototype(x[1:4], K=2, method='median')
x['cluster'] = result[[1]]
pca_iris = princomp(x[,1:4])$scores[,1:2]
plot(pca_iris, t='n')
text(pca_iris, labels=x$Species,col=x$cluster) 

## -----------------------------------------------------------------------------
result[[2]]

## -----------------------------------------------------------------------------
x <- iris[1:100,]
label = cluster.hierarchical(x[,1:4], K=2,method='mean')

## -----------------------------------------------------------------------------
x['cluster'] = label
pca_iris = princomp(x[,1:4])$scores[,1:2]
plot(pca_iris, t='n')
text(pca_iris, labels=x$Species,col=x$cluster)

