##Problem 1

#Plot the given points
x = c(0,3,0,21,23)
y = c(0,0,6,2,2)
plot(x,y)

#Create a matrix to hold all the points, with x and y as rows
m = matrix(nrow=2, ncol=5)
m[1,] = x
m[2,] = y

#We've created our matrix with the x and y coordinates ("genes") as rows
#and the different points ("samples") as columns. This is how most biological
#matrices are organized, but many functions expect it to be the other way
#around, so let's also get the transverse of our matrix
trm = t(m)

#Heatmap and dendrogram. If you haven't already, load the package.
library("pheatmap")
#calculate Euclidean distance between the points
pointsDist = dist(trm, method = "euclidean")
#Draw a heatmap of the distances, with dendrograms from heirarchical clustering
pheatmap(pointsDist, labels_col=seq(from=1, to=5))

##Problem 2

#Run kmeans with k=3. We'll tell the algorithm to choose 5 sets of random
#starting points and give us the most common answer
kclusters = kmeans(trm, centers=3, nstart=5)
#Plot the data, with different colors for the 3 clusters
clus = kclusters$cluster
plot(x, y, col=clus)

## Problem 3

#Generate a matrix of random data for 15 samples, where we've measured the 
#differential RNA levels of 100 genes in every sample. We'll set a seed so
#that everyone's data looks the same.
set.seed(20)
exSample = rnorm(100)
exData = replicate(15, rnorm(100))
#Hierarchical clustering & heatmap
pheatmap(exData, labels_col=seq(from=1, to=15))

#PCA - we'll use a package called FactoMineR for this
library(FactoMineR)
result = PCA(t(exData))