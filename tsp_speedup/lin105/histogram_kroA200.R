library(ggplot2)
library(zoo)
library(dplyr)
library(ggpubr)

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

E_min=29368
tail(sort(data_univ$E, decreasing=TRUE))
tail(sort(data_univ$E, decreasing=FALSE))
yax_min=log(10^(-3))
K=20

data_univ$disp_E = (data_univ$E-E_min)/E_min
data_univ$log_I  = log(data_univ$I)
data_univ$yax    = log(data_univ$disp_E+10^(-3))

data_s = data_univ[which(data_univ$P == 1),]

max_yax=max(data_univ$yax)
c_max_yax=ceiling(max_yax*K)/K
min_yax=min(data_univ$yax)
nbins=ceiling(max_yax-min_yax)*20

G0s   = rep(0,nbins)
G1s   = rep(0,nbins)
Gavgs = rep(0,nbins)
for (i in 1:nbins){
  G0s[i]=min(data_s[which(
  c_max_yax-(i-1)/K>data_s$yax &
  data_s$yax>=c_max_yax-(i)/K),]$I)
  G1s[i]=max(data_s[which(
  c_max_yax-(i-1)/K>data_s$yax &
  data_s$yax>=c_max_yax-(i)/K),]$I)
  Gavgs[i]=mean(data_s[which(
  c_max_yax-(i-1)/K>data_s$yax &
  data_s$yax>=c_max_yax-(i)/K),]$I)
}


sp0     = rep(0,nrow(data_univ))
sp1     = rep(0,nrow(data_univ))
spavg   = rep(0,nrow(data_univ))

for (i in 1:nrow(data_univ)) 
{
  curr_row=data_univ[i,]
  j=as.integer(K*c_max_yax-floor(K*curr_row$yax))
  if (j > 0){
    sp0[i]   = curr_row$I/(G0s[j]*curr_row$P)
    sp1[i]   = curr_row$I/(G1s[j]*curr_row$P)
    spavg[i] = curr_row$I/(Gavgs[j]*curr_row$P)
  } 
}
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

data_univ$sp0   = sp0
data_univ$sp1   = sp1
data_univ$spavg = spavg

data_parallel = data_univ[which(data_univ$P>1),]


plots <- list()

#ITs TIME!
j=1
for (i in c(5,10,20)){
  plot0 <- ggplot(data_parallel[which(data_parallel$P==i),], aes(y=log10(sp0),x=after_stat(density))) + geom_histogram(binwidth=0.1)
  plot0 <- plot0 + facet_grid(cols=vars(Interval))
  plot0 <- plot0 + theme(plot.title=element_text(hjust=0.5))
  plot0 <- plot0 + ggtitle("Interval")

  plot1 <- ggplot(data_parallel[which(data_parallel$P==i),], aes(y=log10(sp1),x=after_stat(density)))+geom_histogram(binwidth=0.1)
  plot1 <- plot1 + facet_grid(cols=vars(Interval))
  plot1 <- plot1 + theme(strip.background = element_blank(), strip.text=element_blank())

  plotavg <- ggplot(data_parallel[which(data_parallel$P==i),], aes(y=log10(spavg),x=after_stat(density)))+geom_histogram(binwidth=0.1)
  plotavg <- plotavg + facet_grid(cols=vars(Interval))
  plotavg <- plotavg + theme(strip.background = element_blank(), strip.text=element_blank())

  plots[[j]] <- ggarrange(plot0, plot1, plotavg,
                     ncol=1,
                     nrow=3)
  j=j+1
}


ggsave(filename="procs5speedup.pdf", plot=plots[[1]], height=16, width=16)
ggsave(filename="procs10speedup.pdf", plot=plots[[2]], height=16, width=16)
ggsave(filename="procs20speedup.pdf", plot=plots[[3]], height=16, width=16)

c_max_p_yax=ceiling(K*max(data_parallel$yax))/K
K_max_p_yax=ceiling(K*max(data_parallel$yax))
K_max_yax=ceiling(K*max(data_univ$yax))
ngraphs=ceiling(K*(c_max_p_yax-yax_min))


data_s$bin <- as.integer(K_max_yax-floor(K*data_s$yax))

min_bin=K_max_yax-K_max_p_yax
hists <- ggplot(data_s[which(data_s$bin>=min_bin),], aes(x=log_I, y=after_stat(density))) + geom_histogram()
hists <- hists + facet_grid(rows=vars(bin))
ggsave(filename="hists_test.pdf", plot=hists, height=16, width=16)
