library(tidyverse)
library(dplyr)
#library(philentropy)


  #load data
data1 <- read.delim("yeast.tsv")
data <- data1[1:500,]
#head(data)

  #cleaning

#sum(is.na(data))
#data %>% 
#  filter(is.na(data))
library(VIM)
aggr(data)

#change to row mean
for (i in 1:ncol(data)){
  cat("Number of NA's in ", names(data)[i], "=", sum(is.na(data[i])), "\n")
  
  data[i] <- ifelse(is.na(data[[i]]),
                   ave(data[[i]], FUN = function(x) mean(x, na.rm = TRUE)),
                   data[[i]])
}


  # hierarchical clustering



  # kmeans clustering

#initialize k random points as centroids
set.seed(123)
create_centroids <- function(df, k = 3){
  
  x <- data.frame(matrix(ncol = 7, nrow = k))
  colnames(x) <- c(names(df))
  
  x[1:k, 1:ncol(x)] <- replicate(ncol(x), runif(k, min(df), max(df)))
  x$cluster <- 1:k
  return(x)
}

centres <- create_centroids(data)

#for each observation, calculate sum squared error to each k centroid

#for each observation, assign centroid that minimizes error 
#we have calculated as its cluster
e1 <- 0

data$cluster <- NA
dis <- function(a,b){
  a <- a[-8]
  b <- b[-8]
  sq <- (a-b)^2
  return(sum(sq))
}

for (i in 1:nrow(data)){
  a <- c(
    dis(centres[1,], data[i,]),
    dis(centres[2,], data[i,]),
    dis(centres[3,], data[i,])
  )
  
  data$cluster[i] <- which(a == min(a))
  e1 <- e1 + min(a)
}

k <- 3
for (i in 1:nrow(data)){
  a <- c()
  for (j in 1:k){
    a <- append(a, dis(centres[j,], data[i,]))
  }
}



# calculate new centroids
c1 <- data %>%
  group_by(cluster) %>%
  summarize(X0hr = mean(X0hr), 
            X9.5hr = mean(X9.5hr), 
            X11.5hr = mean(X11.5hr),
            X13.5hr = mean(X13.5hr),
            X15.5hr = mean(X15.5hr),
            X18.5hr = mean(X18.5hr),
            X20.5hr = mean(X20.5hr)) %>%
  as.data.frame(.)

centres <- centres %>% 
  rows_upsert(c1, by = "cluster")

# iterate 100 times
tol <- 1e-8
it <- 0

e2 <- 0
while (it < 100){
  
  #e step
  for (i in 1:nrow(data)){
    a <- c(
      dis(centres[1,], data[i,]),
      dis(centres[2,], data[i,]),
      dis(centres[3,], data[i,])
    )
    data$cluster[i] <- which(a == min(a))
    e2 <- e2 + min(a)
  }
  
  #m step
  if (abs(e1-e2)>tol){
    c1 <- data %>%
      group_by(cluster) %>%
      summarize(X0hr = mean(X0hr), 
                X9.5hr = mean(X9.5hr), 
                X11.5hr = mean(X11.5hr),
                X13.5hr = mean(X13.5hr),
                X15.5hr = mean(X15.5hr),
                X18.5hr = mean(X18.5hr),
                X20.5hr = mean(X20.5hr)) %>%
      as.data.frame(.)
    
    centres <- centres %>% 
      rows_upsert(c1, by = "cluster")
    
    it <- it +1
    cat("Iteration = ", it, "Change = ", abs(e1-e2),"\n")
    
    e1 <- e2
    e2 <- 0
  }
  else break
  
}


