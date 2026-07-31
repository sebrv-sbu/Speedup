library(ggplot2)
library(zoo)
library(dplyr)
library(ggpubr)
library(reshape2)
library(scales)

data_s <- read.table("s_results", header=TRUE)
data_5 <- read.table("5_results", header=TRUE)
data_10 <- read.table("10_results", header=TRUE)
data_20 <- read.table("20_results", header=TRUE)

data_s$P  = rep(1, nrow(data_s))
data_5$P  = rep(5, nrow(data_5))
data_10$P = rep(10,nrow(data_10))
data_20$P = rep(20,nrow(data_20))
data_parallel <- bind_rows(data_5,data_10,data_20)

data_parallel$lambda = rep(0.000005, nrow(data_parallel))

data_univ<-bind_rows(data_s,data_parallel)
E_min=14382.995939

K=100
#changing the "score" value here.
data_univ$log_I  = log(data_univ$I)
data_univ$yax    = log(data_univ$E/E_min)
yax_min          = 0

max_yax   = max(data_univ$yax)
c_max_yax = ceiling(max_yax*K)/K
min_yax   = min(data_univ$yax)
nbins     = ceiling(max_yax-min_yax)*K
print(nbins)
#Manipulations:
# Since we may still change the number of bins and the
# width of each interval, I will use K to be the number such
# that the width of the energy interval is 1/K.
#G0[j]=min(I | c_max-(j-1)/K>y>=c_max - j/K)
#Thus, G0 is the minimum I for runs with final 
#y value of 
# K*c_max - (j-1) > Ky >= K c_max - j
#since K*c_max is an integer, j is an integer, and 1 is an integer
#this is equivalent to
#floor(Ky) = K c_max - j
# we want to find the i such that G0[i] corresponds to the energy
# y. We thus rearrange:
# j = K c_max - floor(Ky)
# Thus, if data_univ[i] has a value of y, then 
# the bin it is in should be this j. 


K_max_yax     = ceiling(K*max(data_univ$yax))
data_univ$bin <- as.integer(K_max_yax - floor(K*data_univ$yax))

data_parallel = data_univ[which(data_univ$P>1),]
data_s        = data_univ[which(data_univ$P == 1),]

c_max_p_yax = ceiling(K*max(data_parallel$yax))/K
K_max_p_yax = ceiling(K*max(data_parallel$yax))

ngraphs = max(data_s$bin)

#Now we begin!
lambdas <- unique(data_s$lambda)

lambda_frame <- matrix(0.0, ncol=ngraphs*2, nrow=length(lambdas))
#Here, we will abuse the fact that each lambda has 100 runs exactly.
for (i in 1:length(lambdas)){
  for (j in 1 : 100){
    curr_row                    = data_s[(i-1)*100+j,]
    bin                         = curr_row$bin
    lambda_frame[i,bin]         = 1/100 + lambda_frame[i,bin]
    lambda_frame[i,ngraphs+bin] = curr_row$I + lambda_frame[i, ngraphs+bin]
  }
}
for (i in 1:length(lambdas)){
  for (j in 1:ngraphs){
    if (lambda_frame[i,ngraphs+j] != 0)
      lambda_frame[i,ngraphs+j] = lambda_frame[i,ngraphs+j]/(100*lambda_frame[i,j])
    else 
      lambda_frame[i,ngraphs+j] = NA
  }
}

lambda_data            <- as.data.frame(lambda_frame)
bin_cols               <- seq(from=c_max_yax, by=-1/K, length.out=max(data_s$bin))
density_cols           <- paste("iterations_", bin_cols, sep='')
colnames(lambda_data)  <- c(bin_cols, density_cols)
lambda_data$log_lambda <- log(lambdas)

melted_lambda_data     <- melt(lambda_data, id=c("log_lambda"))
#This is extremely hard to follow. Basically, I am trying to go from the following format:
#lambda | bin_0.5_density | bin_0_density | .... | bin_-4.6_density | bin_0.5_iterations | ... | bin_-4.6_iterations |
# 0.1          0              0.1                    0                    NA                        NA
# to
# lambda | bin | iterations 
#   0.1     0.5    NA
#   0.1     0      1000
# Something like this. 
# When we first melt the data, we then get
# lambda | var            | value
#  0.1     0.5                0
#  ...      ...              ...
#  0.1     0.5_iterations |   NA 
# ...
# But note that in the actual code, instead of bin_0.5_density, it is 0.5, so we don't
# have to do any string comprehension. 

#Here we cut the melted lambda data in half, the top half is the density stuff,
# the bottom the iteration stuff.
index                 <- seq.int(nrow(melted_lambda_data)/2)
lambda_data           <- melted_lambda_data[index,]
its_m_lambda_data     <- melted_lambda_data[-index,]
colnames(lambda_data) <- c("log_lambda", "bin", "density")
lambda_data$log_I     <- log(melted_lambda_data[-index,]$value)
lambda_data$bin       <- as.numeric(paste(lambda_data$bin))

#Here, we create the function to "transform" the bin into something with the same
#y axis as log(iterations). 
c_bin         <- min(lambda_data[!is.na(lambda_data$log_I),]$log_I)
m_bin         <- max(lambda_data[!is.na(lambda_data$log_I),]$log_I)-c_bin
transform_bin <- function(x){
  x=x-min(lambda_data$bin)
  x=x/(max(lambda_data$bin)-min(lambda_data$bin))
  return (m_bin*x + c_bin)
}
inverse_bin   <- function(x){
  x=(x-c_bin)/m_bin
  x=x*(max(lambda_data$bin)-min(lambda_data$bin))
  x=x+min(lambda_data$bin)
  return(x)
}
lambda_data$trans_bin <- transform_bin(lambda_data$bin)
#Now, we have to "transform" the density into something with the same range as bin.
c_density         <- min(lambda_data$bin)
m_density         <- max(lambda_data$bin)-c_density

plot <- ggplot(data=lambda_data[which(lambda_data$density >0),], aes(x=log_lambda)) 
plot <- plot + geom_point(aes(y=transform_bin(bin), fill=density), shape=22, size=1, stroke=NA)
plot <- plot + geom_point(size=0.5 ,aes(y=log_I,color=bin))
plot <- plot + scale_color_gradientn(colors=c('orange','purple','green'), values=rescale(c(0,0.07,max(lambda_data$bin))))
plot <- plot + scale_fill_gradientn(colors=c("turquoise3","yellow","red4"), values=rescale(c(0,0.2,0.6)))
plot <- plot + scale_y_continuous(name="log iterations", sec.axis=sec_axis( transform=~inverse_bin(.), name="energy score"))
ggsave(filename="lambdas_iterations_serial_logfracscore.pdf", plot=plot)
#Now, lets compare the Intervals.
int_frames <- list()
intervals_l <- list()
k=1
procs=c(5,10,20)
for (proc in procs){
  data <- (data_univ[which(data_univ$P == proc),])
  print(nrow(data))
  intervals <- unique(data$Interval)
  intervals_l[[k]]=intervals
  int_frame <- matrix(0.0, ncol=ngraphs*2, nrow=length(intervals))
  for (i in 1:length(intervals)){
    curr_data <- data[which(data$Interval==intervals[i]),]
    for (j in 1:nrow(curr_data)){
        curr_row <- curr_data[j,]
        bin                      = curr_row$bin
        int_frame[i,bin]         = 1/nrow(curr_data)+int_frame[i,bin]
        int_frame[i,ngraphs+bin] = curr_row$I + int_frame[i,ngraphs+bin]


    }
  }
  for (i in 1 : length(intervals)){
    for (j in 1 : ngraphs){
      if (int_frame[i,ngraphs+j] != 0)
        int_frame[i, ngraphs+j]=int_frame[i,ngraphs+j]/(nrow(data[which(data$Interval==intervals[i]),])*int_frame[i,j])
      else
        int_frame[i,ngraphs+j]=NA
    }
  }
  int_frames[[k]]=int_frame
  k = k+1
}
k=1

#for an explanation of the below, refer to the comment for the serial data construction,
#and replace lambda with interval.
for (proc in procs){
  intervals      <- intervals_l[[k]]
  curr_int_frame <- int_frames[[k]]
  int_data       <- as.data.frame(curr_int_frame)

  bin_cols           <- seq(from=c_max_yax, by=-1/K, length.out=max(data_s$bin))
  density_cols       <- paste("iterations_", bin_cols, sep='')
  colnames(int_data) <- c(bin_cols, density_cols)
  int_data$log_int   <- log(intervals)

  melted_int_data     <- melt(int_data, id=c("log_int"))
  index               <- seq.int(nrow(melted_int_data)/2)
  int_data            <- melted_int_data[index,]
  its_m_int_data      <- melted_int_data[-index,]
  colnames(int_data)  <- c("log_int", "bin", "density")
  int_data$log_I      <- log(melted_int_data[-index,]$value)
  int_data$bin        <- as.numeric(paste(int_data$bin))
  
  
  c_bin         <- min(int_data[!is.na(int_data$log_I),]$log_I)
  m_bin         <- (max(int_data[!is.na(int_data$log_I),]$log_I)-c_bin)/4
  int_data <- int_data[!is.na(int_data$log_I),]
  transform_bin <- function(x){
    x=x-min(int_data$bin)
    x=x/(max(int_data$bin)-min(int_data$bin))
    return (m_bin*x + c_bin)
    }

  inverse_bin   <- function(x){
    x=(x-c_bin)/m_bin
    x=x*(max(int_data$bin)-min(int_data$bin))
    x=x+min(int_data$bin)
    return(x)
    } 

  plot <- ggplot(data=int_data[which(int_data$density >0),], aes(x=log_int)) 
  plot <- plot + geom_point(aes(y=transform_bin(bin), size=density, fill=density), shape=22, stroke=NA)
  plot <- plot + scale_size_continuous(range=c(0.1, 2))
  plot <- plot + scale_fill_gradientn(colors=c("turquoise3","red4"))
  plot <- plot + geom_point(size=0.5 ,aes(y=log_I,color=bin))
  plot <- plot + scale_color_gradientn(colors=c('orange','purple','green'))
  plot <- plot + scale_y_continuous(name="log iterations", sec.axis=sec_axis( transform=~inverse_bin(.), name="energy score"))
  ggsave(filename=paste("int_iterations_", proc, "_logfracscore.pdf", sep=''), plot=plot)

  k=k+1
}

