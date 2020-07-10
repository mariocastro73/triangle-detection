##------ Sat May 23 11:38:00 2020 ------##
#
# This code takes csv files created with the python scripts 
# and make a "batch" Voronoi tesselation analysis

############## Libraries 
library(dplyr)
library(shape)
library(RColorBrewer)
library(MASS)
library(RTriangle)

############## Auxiliary functions 

# Takes a dataframe with coordinates x and y and compute all the metrics related to a Voronoi tesselation (distribution of neighbors and areas)
datavoronoi <- function(data) {
  colnames(data) <- c('x','y') # Rename columns
  n <- length(data$x)
  p <- pslg(data)
  t <- triangulate(p)
  E <- t$E
  nn <-matrix(0,nrow=n,ncol=n)
  for(i in 1:dim(t$E)[1]) {
    nn[E[i,1],E[i,2]] = 1 
    nn[E[i,2],E[i,1]] = 1 
  }
  neigh <- colSums(nn)
  r <- range(neigh)
  h <- hist(neigh,r[1]:r[2],plot = FALSE)
  p6 <- h$density[which(h$breaks==6)-1]
  mu2 <- var(neigh)
  
  lem <- c(p6,mu2) # Store lemaitre coefficients (height of distribution at n=6 and variance of distribution)
  p61 <- seq(0.1,.7,.01) # Auxiliary vector for plotting 
  p62 <- seq(0.1,1,.01) # Auxiliary vector for plotting 
  mu2 <- seq(0,4,.01) # Auxiliary vector for plotting 
  
  output <- list(lem=lem,neigh=neigh,mu2=mu2,p61=p61,p62=p62,data=data,triangles=t,nn=nn) # Store data into list
  return(output)  # Return list
}



# Find clusters: If flag=0 it creates randomly N coordinates to illustrate how the code works.
# The thresh value is fixed according to the experimental scale
mycluster <- function(dataset,N=10,thresh=1,flag=1) {
  if(flag==0) {
    x <- rt(N,2)
    y <- rt(N,2)
    data <- data.frame(x,y)
  }
  else {
    data <- dataset$data
    N <- dim(data)[1]
    x <- data$x
    y <- data$y
    }
  d <- dist(data)
  h <- hclust(d,method='single')
  c <- cutree(h,h=thresh)
  sm <- min(c(x,y))
  sM <- max(c(x,y))
  if(flag==0) {
    plot(data$x,data$y,pch=19,col=c,xlim=c(sm,sM),ylim=c(sm,sM),cex=.5,xlab='X coordinate',ylab='Y coordinate')
    for(i in 1:N) draw.circle(data$x[i],data$y[i],thresh/2)
  }
  else {
    plot(dataset$triangles,pch=19,cex=.1,lwd=.5)
    points(data$x,data$y,pch=19,col=c,cex=.7)
  }
  return(list(x=x,y=y,c=c,h=h,dist=d))
}



########### File preprocessing

f <- read.csv('Ca2Na.csv')
data <- f[,2:3]
output <- datavoronoi(data)
saveRDS(output,'output-CaCl2.rda')

f <- read.csv('Ca2Na3nM.csv')
data <- f[,2:3]
output <- datavoronoi(data)
saveRDS(output,'output-CaCl2-3nM.rda')

f <- read.csv('MgK.csv')
data <- f[,2:3]
output <- datavoronoi(data)
saveRDS(output,'output-Mg2+K+.rda')

f <- read.csv('MgLi.csv')
data <- f[,2:3]
output <- datavoronoi(data)
saveRDS(output,'output-Mg2+Li+.rda')

f <- read.csv('MgNa.csv')
data <- f[,2:3]
output <- datavoronoi(data)
saveRDS(output,'output-Mg2+Na+.rda')



##################### COMPARISON PLOTS
CaCl <- readRDS('output-CaCl2.rda')
CaCl3nM <- readRDS('output-CaCl2-3nM.rda')
MgK <- readRDS('output-Mg2+K+.rda')
MgLi <- readRDS('output-Mg2+Li+.rda')
MgNa <- readRDS('output-Mg2+Na+.rda')
Hex <- readRDS('output-Hex.rda')
par(mar=c(5.1,5.1,4.1,2.1))
# Color pallette
colCaCl <- rgb(red = 0.65, green = 0.8, blue = .9, alpha = 0.8)
colCaCl3 <-  rgb(red = 0.5, green = 0.2, blue = .7, alpha = 1)
colMgK <- rgb(red = 0.71, green = .87, blue = 0.55, alpha = 1)
colMgNa <- rgb(red = 1, green = 1, blue = .61, alpha = 0.8)
colMgLi <- rgb(red = 1, green = .73, blue = 0.42, alpha = 0.8)

# Plot
plotlemaitre(CaCl,ymax=3,ymin=0,defaultcol = colCaCl,mycex=3,mypch=21)
plotlemaitre(CaCl3nM,ymax=3,ymin=0,defaultcol =colCaCl3,add=T,mycex=3,mypch=21)
plotlemaitre(MgK,defaultcol = colMgK,add=T,mycex=3,mypch=22)
plotlemaitre(MgNa,defaultcol = colMgNa,add=T,mycex=3,mypch=23)
plotlemaitre(MgLi,defaultcol =colMgLi ,add=T,mycex=3,mypch=24)
plotlemaitre(Hex,defaultcol = 'white',add=T,mycex=3,mypch=25)
points(1,0,col = 1,bg='blue',cex=3,pch=25)
# Legend
legend('topright',legend=c(expression(paste(Ca^"2+","/",Na^"+")),
                           expression(paste(Ca^"2+","/",Na^"+"," (3nM)")),
                           expression(paste(Mg^"2+","/",K^"+")),
                           expression(paste(Mg^"2+","/",Na^"+")),
                           expression(paste(Mg^"2+","/",Li^"+")),
                           "Finite honeycomb",
                           "Infinite honeycomb"), col=1,
       pt.bg=c(colCaCl,colCaCl3,colMgK,colMgNa,colMgLi,'white','blue'),
       pch=c(21,21,22,23,24,25,25),cex=1.5)

dev.copy2pdf(file='comparison-ions.pdf') # Save to file


##### CLUSTERING #####
## Hex
clsHex <- mycluster(Hex,thresh = 15)
frac <- cluster.stats(clsHex$c)
title(sprintf("Fraction of triangles in largest cluster: %.0f%%",frac))
dev.copy2pdf(file='clusters-Hex.pdf')

# clsCaCl <- mycluster.orig(CaCl$data,thresh = 15)
clsCaCl <- mycluster(CaCl,thresh = 15)
frac <- cluster.stats(clsCaCl$c)
title(sprintf("Fraction of triangles in largest cluster: %.0f%%",frac))
dev.copy2pdf(file='clusters-CaCl.pdf')

# clsMgK <- mycluster.orig(MgK$data,thresh = 15)
clsMgK <- mycluster(MgK,thresh = 15)
frac <- cluster.stats(clsMgK$c)
title(sprintf("Fraction of triangles in largest cluster: %.0f%%",frac))
dev.copy2pdf(file='clusters-MgK.pdf')

# clsMgLi <- mycluster.orig(MgLi$data,thresh = 15)
clsMgLi <- mycluster(MgLi,thresh = 15)
frac <- cluster.stats(clsMgLi$c)
title(sprintf("Fraction of triangles in largest cluster: %.0f%%",frac))
dev.copy2pdf(file='clusters-MgLi.pdf')

# clsMgNa <- mycluster.orig(MgNa$data,thresh = 15)
clsMgNa <- mycluster(MgNa,thresh = 15)
frac <- cluster.stats(clsMgNa$c)
title(sprintf("Fraction of triangles in largest cluster: %.0f%%",frac))
dev.copy2pdf(file='clusters-MgNa.pdf')
#-------------------
clsCaCl3nM <- mycluster(CaCl3nM,thresh = 15)
frac <- cluster.stats(clsCaCl$c)
title(sprintf("Fraction of triangles in largest cluster: %.0f%%",frac))
dev.copy2pdf(file='clusters-CaCl-3nM.pdf')

