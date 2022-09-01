#load data
data <- read.delim("yeast.tsv")
#data <- data1[1:500,]

#cleaning
library(VIM)
aggr(data)

for (i in 1:ncol(data)){
  cat("Number of NA's in ", names(data)[i], "=", sum(is.na(data[i])), "\n")
  
  data[i] <- ifelse(is.na(data[[i]]),
                    ave(data[[i]], FUN = function(x) mean(x, na.rm = TRUE)),
                    data[[i]])
}


set.seed(123)

create_centroids <- function(df, k){
  #Returns k centroids with randomly sampled values
  
  x <- data.frame(matrix(ncol = 7, nrow = k))
  colnames(x) <- c(names(df))
  
  x[1:k, 1:ncol(x)] <- replicate(ncol(x), runif(k, min(df), max(df)))
  x$cluster <- 1:k
  return(x)
}
dis <- function(a,b){
  #Returns Euclidean distance between vectors a,b
  a <- a[-8]
  b <- b[-8]
  sq <- (a-b)^2
  return(sqrt(sum(sq)))
}


k_means <- function(df, k){
  #initalize k random points as centroids
  centres <- create_centroids(data,k)
  #print(centres)
  
  #initial e-step
  e1 <- 0
  for (i in 1:nrow(data)){
    a <- c()
    for (j in 1:k){
      a <- append(a, dis(centres[j,], data[i,]))
    }
    data$cluster[i] <- which(a == min(a))
    e1 <- e1 + min(a)
  }

  #cat("clusters at initial estep: ", data$cluster)
  
  
  
  #initial m-step
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
  
  tol <- 1e-8
  max_iter <- 1e2
  it <- 1
  e2 <- 0
  while (it < max_iter){
    
    #e-step
    for (i in 1:nrow(data)){
      a <- c()
      for (j in 1:k){
        a <- append(a, dis(centres[j,], data[i,]))
      }
      data$cluster[i] <- which(a == min(a))
      e2 <- e2 + min(a)
    }


    
    
    #m-step
    #cat("Iteration = ", it, "Change = ", abs(e1-e2), "e1 = ", e1, "e2 =", e2,"\n")
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
      
      it <- it+1
      e1 <- e2
      e2 <-0
    }
    else break
  }
  
  cat("Iterations: ", it)
  return(centres)
  
}


























