library(ggplot2)
library(zoo)
library(dplyr)
library(ggpubr)
library(reshape2)
library(scales)

data_s    <- read.table("../s_results", header=TRUE)
print("s")
data_5_c  <- read.table("../5_c_results", header=TRUE)
print("5_c")
data_10_c <- read.table("../10_c_results", header=TRUE)
print("10_c")
data_20_c <- read.table("../20_c_results", header=TRUE)
print("20_c")
data_20_h <- read.table("../20_h_results", header=TRUE)
print("20_h")

data_s$P      = rep(1, nrow(data_s))
data_5_c$P    = rep(5, nrow(data_5_c))
data_10_c$P   = rep(10,nrow(data_10_c))
data_20_c$P   = rep(20,nrow(data_20_c))
data_20_h$P   = rep(20,nrow(data_20_h))

data_parallel_c <- bind_rows(data_5_c,data_10_c,data_20_c)
names(data_parallel_c)[names(data_parallel_c) == "Interval"] <- "m_0"
data_parallel_c$lambda = rep(0.000005, nrow(data_parallel_c))
data_s$scheme=rep("serial", nrow(data_s))
data_parallel_c$scheme=rep("cdr", nrow(data_parallel_c))
data_parallel_h <- data_20_h #later it will have more!
data_parallel_h$scheme <- rep("hsa",nrow(data_parallel_h))
data_parallel_h$lambda <- rep(0.000005,nrow(data_parallel_h))

data_univ<-bind_rows(data_s,data_parallel_c, data_parallel_h)
E_min=29368

K_E=100
#changing the "score" value here.
data_univ[which(data_univ$scheme != "hsa"),]$I = data_univ[which(data_univ$scheme != "hsa"),]$I-100000 - 100000/data_univ[which(data_univ$scheme != "hsa"),]$P
data_univ[which(data_univ$scheme == "hsa" & data_univ$P==20),]$I = data_univ[which(data_univ$scheme =="hsa" & data_univ$P==20),]$I-100000
-100000/4
data_univ$log_I  = log(data_univ$I)
data_univ$yax    = data_univ$E/E_min

max_yax=max(data_univ$yax)
c_max_yax=ceiling(max_yax*K_E)/K_E
min_yax=min(data_univ$yax)
nbin_Es=ceiling(max_yax-min_yax)*K_E
#Manipulations:
# Since we may still change the number of bin_Es and the
# width of each interval, I will use K_E to be the number such
# that the width of the energy interval is 1/K_E.
#G0[j]=min(I | c_max-(j-1)/K_E>y>=c_max - j/K_E)
#Thus, G0 is the minimum I for runs with final 
#y value of 
# K_E*c_max - (j-1) > K_Ey >= K_E c_max - j
#since K_E*c_max is an integer, j is an integer, and 1 is an integer
#this is equivalent to
#floor(K_Ey) = K_E c_max - j
# we want to find the i such that G0[i] corresponds to the energy
# y. We thus rearrange:
# j = K_E c_max - floor(K_Ey)
# Thus, if data_univ[i] has a value of y, then 
# the bin_E it is in should be this j. 


K_E_max_yax=ceiling(K_E*max(data_univ$yax))
data_univ$bin_E <- as.integer(K_E_max_yax - floor(K_E*data_univ$yax))

data_parallel = data_univ[which(data_univ$P>1),]
data_s = data_univ[which(data_univ$P == 1),]

c_max_p_yax=ceiling(K_E*max(data_parallel$yax))/K_E
K_E_max_p_yax=ceiling(K_E*max(data_parallel$yax))

ngraphs=max(data_s$bin_E)

lambda_int_model_for_print <- lm(data_s[which(data_s$lambda<0.01),]$log_I ~ log(data_s[which(data_s$lambda <0.01),]$lambda ))
summary(lambda_int_model_for_print)
lambda_int_model <- lm(data_s[which(data_s$lambda<0.01),]$log_I ~ log10(data_s[which(data_s$lambda <0.01),]$lambda ))

#Now we begin!
lambdas <- unique(data_s$lambda)

lambda_frame <- matrix(0.0, ncol=ngraphs*2, nrow=length(lambdas))
#Here, we will abuse the fact that each lambda has 100 runs exactly.
for (i in 1:length(lambdas)){
  curr_data <- data_s[which(data_s$lambda == lambdas[i]),]
  for (j in 1 : nrow(curr_data)){
    curr_row                       = curr_data[j,]
    bin_E                          = curr_row$bin_E
    lambda_frame[i,bin_E]          = 1/nrow(curr_data)+ lambda_frame[i,bin_E]
    lambda_frame[i,ngraphs+bin_E]  = curr_row$I + lambda_frame[i, ngraphs+bin_E]
  }
}
for (i in 1:length(lambdas)){
  for (j in 1:ngraphs){
    if (lambda_frame[i,ngraphs+j] != 0)
      lambda_frame[i,ngraphs+j] = lambda_frame[i,ngraphs+j]/(nrow(data_s[which(data_s$lambda==lambdas[i]),])*lambda_frame[i,j])
    else 
      lambda_frame[i,ngraphs+j]=NA
  }
}

lambda_data <- as.data.frame(lambda_frame)
bin_E_cols <- seq(from=c_max_yax, by=-1/K_E, length.out=max(data_s$bin_E))
density_cols <- paste("iterations_", bin_E_cols, sep='')
colnames(lambda_data) <- c(bin_E_cols, density_cols)
lambda_data$log_lambda=log10(lambdas)

melted_lambda_data <- melt(lambda_data, id=c("log_lambda"))
#This is extremely hard to follow. Basically, I am trying to go from the following format:
#lambda | bin_E_0.5_density | bin_E_0_density | .... | bin_E_-4.6_density | bin_E_0.5_iterations | ... | bin_E_-4.6_iterations |
# 0.1          0              0.1                    0                    NA                        NA
# to
# lambda | bin_E | iterations 
#   0.1     0.5    NA
#   0.1     0      1000
# Something like this. 
# When we first melt the data, we then get
# lambda | var            | value
#  0.1     0.5                0
#  ...      ...              ...
#  0.1     0.5_iterations |   NA 
# ...
# But note that in the actual code, instead of bin_E_0.5_density, it is 0.5, so we don't
# have to do any string comprehension. 

#Here we cut the melted lambda data in half, the top half is the density stuff,
# the bottom the iteration stuff.
index <- seq.int(nrow(melted_lambda_data)/2)
lambda_data <- melted_lambda_data[index,]
its_m_lambda_data <- melted_lambda_data[-index,]
colnames(lambda_data) <- c("log_lambda", "bin_E", "density")
lambda_data$log_I <- log(melted_lambda_data[-index,]$value)
lambda_data$bin_E <- as.numeric(paste(lambda_data$bin_E))

#Here, we create the function to "transform" the bin_E into something with the same
#y axis as log(iterations). 
c_bin_E         <- min(lambda_data[!is.na(lambda_data$log_I),]$log_I)
m_bin_E         <- max(lambda_data[!is.na(lambda_data$log_I),]$log_I)-c_bin_E
transform_bin_E <- function(x){
  x=x-min(lambda_data$bin_E)
  x=x/(max(lambda_data$bin_E)-min(lambda_data$bin_E))
  return (m_bin_E*x + c_bin_E)
}
inverse_bin_E <- function(x){
  x=(x-c_bin_E)/m_bin_E
  x=x*(max(lambda_data$bin_E)-min(lambda_data$bin_E))
  x=x+min(lambda_data$bin_E)
  return(x)

}
lambda_data$trans_bin_E <- transform_bin_E(lambda_data$bin_E)
#Now, we have to "transform" the density into something with the same range as bin_E.
c_density         <- min(lambda_data$bin_E)
m_density         <- max(lambda_data$bin_E)-c_density

plot <- ggplot(data=lambda_data[which(lambda_data$density >0),], aes(x=log_lambda)) 
plot <- plot + geom_point(aes(y=transform_bin_E(bin_E), fill=density), shape=22, size=1, stroke=NA)
plot <- plot + geom_point(size=0.5 ,aes(y=log_I,color=bin_E))
plot <- plot + scale_color_gradientn(name="E/E_min",colors=c('orange','purple','green'),values=rescale(c(1,1.07,max(lambda_data$bin_E))))
plot <- plot + scale_fill_gradientn(colors=c("#66cc00","#cc3399", "#43123c"), values=rescale(c(0,0.2,0.6)))
#plot <- plot + scale_fill_gradientn(colors=c("turquoise3","yellow","red4"), values=rescale(c(0,0.2,0.6)))
plot <- plot + scale_y_continuous(name="log iterations", sec.axis=sec_axis( transform=~inverse_bin_E(.), name="energy score"))
plot <- plot + scale_x_continuous(name="log_10(lambda)")
plot <- plot + geom_abline(intercept=coef(lambda_int_model)[1], slope=coef(lambda_int_model)[2])
plot <- plot + theme(element_text(size=20))
ggsave(filename="lambdas_iterations_serial_fracscore.pdf", plot=plot,height=5, width=5*1.5)
#Now, we zoom in on a certain region. We will zoom in to lambda <0.001:

cut_data <- lambda_data[which(lambda_data$log_lambda < log10(0.001)),]
c_bin_E         <- min(cut_data[!is.na(cut_data$log_I),]$log_I)
m_bin_E         <- max(cut_data[!is.na(cut_data$log_I),]$log_I)-c_bin_E
transform_bin_E <- function(x){
  x=x-min(cut_data[!is.na(cut_data$log_I),]$bin_E)
  x=x/(max(cut_data[!is.na(cut_data$log_I),]$bin_E)-min(cut_data[!is.na(cut_data$log_I),]$bin_E))
  return (m_bin_E*x + c_bin_E)
}
inverse_bin_E <- function(x){
  x=(x-c_bin_E)/m_bin_E
  x=x*(max(cut_data[!is.na(cut_data$log_I),]$bin_E)-min(cut_data[!is.na(cut_data$log_I),]$bin_E))
  x=x+min(cut_data[!is.na(cut_data$log_I),]$bin_E)
  return(x)

}

plot <- ggplot(data=cut_data[which(cut_data$density >0),], aes(x=log_lambda)) 
plot <- plot + geom_point(aes(y=transform_bin_E(bin_E), fill=density), shape=22, size=1, stroke=NA)
plot <- plot + geom_point(size=0.5 ,aes(y=log_I,color=bin_E))
plot <- plot + scale_color_gradientn(name="E/E_min",colors=c('orange','purple','green'))#,values=rescale(c(1,1.07,max(lambda_data$bin_E))))
plot <- plot + scale_fill_gradientn(colors=c("#66cc00","#cc3399", "#43123c"), values=rescale(c(0,0.2,0.6)))
plot <- plot + scale_y_continuous(name="log iterations", sec.axis=sec_axis( transform=~inverse_bin_E(.), name="energy score"))
plot <- plot + scale_x_continuous(name="log_10(lambda)")
#plot <- plot + geom_abline(intercept=coef(lambda_int_model)[1], slope=coef(lambda_int_model)[2])
ggsave(filename="zoomed_lambdas_iterations_serial_fracscore.pdf", plot=plot,height=5, width=5*1.5)

#Now, lets compare the Intervals.
#ITS TIME FOR MODIFICATION
#We are going to be redefining the bin_Es here. We are going to introduce 
#Iteration Bins. 
#We are going to use the following function: First, we convert
#f(log_I) : log_I       -> [n_I]
#g(E)     : E           -> [n_E]
#h(i,j)   : [n_I]X[n_E] -> [n_I*n_E]
#We use "pairing functions"

K_I <- 20
inverse_pairing_func <- function(z){
  w <- floor( (sqrt(8*z+1)-1)/2 )
  t <- (w**2 + w)/2
  y <- z - t
  x <- w - y
  return (c(x,y))
}

comp_data <- data_univ[which(data_univ$P > 1 | (data_univ$lambda < 8.5e-6 & data_univ$lambda > 8.4e-7)),]
subtract_bin <- min(comp_data$bin_E)
comp_data$bin_E <- comp_data$bin_E - subtract_bin
print(min(comp_data$E))
print(max(comp_data$E))
print(max(comp_data$bin_E))
print(min(comp_data$bin_E))
print(subtract_bin)
n_E <- max(comp_data$bin_E)
comp_data$log_I <- comp_data$log_I + log(comp_data$P)
max_log_I <- max(comp_data$log_I)
min_log_I <- min(comp_data$log_I)
c_max_log_I <- ceiling(max_log_I*K_I)/K_I
K_I_max_log_I <- ceiling(K_I*max_log_I)
comp_data$bin_I <- as.integer(K_I_max_log_I - floor(K_I*comp_data$log_I))
n_I <- max(comp_data$bin_I)
comp_data$bin <- (((comp_data$bin_I + comp_data$bin_E +1)*(comp_data$bin_I + comp_data$bin_E))/2)+comp_data$bin_I

ngraphs <- ((n_E + n_I + 1)*(n_E + n_I))/2 + n_I


int_frames <- list()
intervals_l <- list()
schemes <- list()
procs_l <- list()
k=1

lambdas <- unique(comp_data$lambda)
data <- comp_data[which(comp_data$P==1),]
int_frame <- matrix(0.0, ncol=ngraphs+3, nrow=length(lambdas))
for (i in 1:length(lambdas)){
  curr_data <- data[which(data$lambda == lambdas[i]),]
  int_frame[i,1] <- lambdas[i]
  int_frame[i,2] <- NA
  int_frame[i,3] <- NA
  for (l in 1:nrow(curr_data)){
    curr_row <- curr_data[l,]
    bin <- curr_row$bin
    int_frame[i,bin+3] = 1/nrow(curr_data) + int_frame[i, bin + 3]
  }
}
int_frames[[k]] <- int_frame
schemes[[k]] <- "serial"
procs_l[[k]] <- 1
k = k+1

procs=c(5,10,20)
for (proc in procs){
  data <- (comp_data[which(comp_data$P == proc & comp_data$scheme == "cdr"),])
  intervals <- unique(data$m_0)
  intervals_l[[k]]=intervals
  int_frame <- matrix(0.0, ncol=ngraphs+3, nrow=length(intervals))
  for (i in 1:length(intervals)){
    curr_data <- data[which(data$m_0==intervals[i]),]
    int_frame[i,1] <- 5e-6
    int_frame[i,2] <- intervals[i]
    int_frame[i,3] <- NA
    for (j in 1:nrow(curr_data)){
        curr_row <- curr_data[j,]
        bin                      = curr_row$bin
        int_frame[i,bin+3]         = 1/nrow(curr_data)+int_frame[i,bin+3]
    }
  }
  int_frames[[k]]=int_frame
  schemes[[k]] = "cdr"
  procs_l[[k]] = proc
  k = k+1
}
procs = c(20)
for (proc in procs){
  data <- comp_data[which(comp_data$P == proc & comp_data$scheme == "hsa"),]
  m_0s <- unique(data$m_0)
  m_1s <- unique(data$m_1)
  l_m_1s <- length(m_1s)
  int_frame <- matrix(0.0, ncol=ngraphs + 3, nrow=length(m_0s)*l_m_1s)
  for (i in 1:length(m_0s)){
    for (l in 1:l_m_1s){
      row <- (i-1)*l_m_1s + l
      m_0 <- m_0s[i]
      m_1 <- m_1s[l]
      curr_data <- data_parallel_h[which(data_parallel_h$m_0==m_0 & data_parallel_h$m_1==m_1),]
      int_frame[row,1] = 5e-6
      int_frame[row,2] = m_0
      int_frame[row,3] = m_1
      for (j in 1:nrow(curr_data)){
        curr_row <- curr_data[j,]
        bin                = curr_row$bin
        int_frame[row,bin+2] = 1/nrow(curr_data)+int_frame[row,bin+2]
        }
     }
  }
  int_frames[[k]]=int_frame
  schemes[[k]] <- "hsa"
  procs_l[[k]] = proc
  k = k+1
}


#We are now completely modifying this. Refer to the explanation in the section on HSA.
columns <- c("lambda", "m_0", "m_1", "scheme", "bin_E", "bin_I", "P")
tot_int_data <- data.frame(matrix(ncol=length(columns), nrow=0))
colnames(tot_int_data) <- c("scheme","log_int", "bin_E", "density", "log_I","P")
for (i in (1:(k-1))){
  curr_int_frame <- int_frames[[i]]
  int_data       <- as.data.frame(curr_int_frame)
  colnames(int_data) <- c("lambda", "m_0", "m_1", seq(1,ngraphs, by=1))
  melted_int_data     <- melt(int_data, id.vars=c("lambda", "m_0", "m_1"))
  colnames(melted_int_data) <- c("lambda", "m_0", "m_1", "bin", "density")
  melted_int_data$bin <- as.integer(paste(melted_int_data$bin))
  melted_int_data[,c("bin_E", "bin_I")] <- inverse_pairing_func(melted_int_data$bin)
  melted_int_data <- melted_int_data[which(melted_int_data$bin_E <= n_E & melted_int_data$bin_I <= n_I),]
  melted_int_data$scheme <- schemes[[i]]
  melted_int_data$P <- procs_l[[i]]

  tot_int_data <- rbind(melted_int_data, tot_int_data)
}
  transform_bin_E <- function(bin){
    orig_bin <- bin + subtract_bin - 1
    return((K_E_max_yax - orig_bin)/K_E)
    }
  transform_bin_I <- function(bin){
    return ((K_I_max_log_I - bin)/K_I)
  }
plots <- list()
k=1
ints <- c(10,20,40,50,80,100,200)
print(K_E_max_yax)
plot_data <-tot_int_data[which(tot_int_data$scheme == "cdr" & tot_int_data$m_0 %in% ints),]
comp_data <- plot_data[which(plot_data$density>0),]
plot_data <- plot_data[which(plot_data$bin_I <= max(comp_data$bin_I) & plot_data$bin_I >= min(comp_data$bin_I) 
                             & plot_data$bin_E >= min(comp_data$bin_E) & plot_data$bin_E <= max(comp_data$bin_E)),]
plot <- ggplot(data=plot_data, aes(x=bin_I, y=bin_E, fill=density))
#  plot <- plot + geom_point(aes(y=transform_bin_E(bin_E), size=density, fill=density), shape=22, stroke=NA)
#  plot <- plot + scale_size_continuous(range=c(0.1, 2))
#  plot <- plot + scale_fill_gradientn(colors=c("turquoise3","red4"))
plot <- plot + geom_tile()
plot <- plot + facet_grid(rows=vars(P), cols=vars(m_0))
#stop here for now...
#plot <- plot + scale_color_gradientn(name="E/E_min",colors=c('orange','purple','green'), values=rescale(c(1,2.125/2, 1.125)), limits=c(1,1.125) )
plot <- plot + scale_x_reverse(name="log(I*P)",breaks=c(20, 34, 48),expand=c(0,0),labels= lapply(c(20,34,48), transform_bin_I))
plot <- plot + scale_y_reverse(name="E/E_min",breaks=c(1,5,9,13), expand=c(0,0), labels=lapply(c(1, 5, 9, 13), transform_bin_E))
plot <- plot + scale_fill_gradientn(colors=c("white","orange","purple"), values=rescale(c(0,0.05,max(plot_data$density))))
plot <- plot + theme(panel.spacing= unit(.01, "lines"),
                         panel.border = element_rect(color="black", fill=NA, linewidth=0.2))
ggsave(filename="cdr_heatmap.pdf", plot=plot, height=5, width=5*2)

data_c <- data_univ[which(data_univ$scheme == "cdr"),]
print(max(tot_int_data[which(tot_int_data$density > 0),]$density))
plot <- ggplot(data=tot_int_data[which(tot_int_data$density > 0),], aes(x=bin_I, y=bin_E, color=density))
plot <- plot + geom_point()
plot <- plot + facet_grid(rows=vars(P))
plot <- plot + scale_color_gradientn(name="density", colors=c("#66cc00","#cc3399", "#43123c"))
ggsave(filename="int_energy_cdr_fracscore.pdf", plot=plot)
#The above graphs look like shit when we put both factors in, so we are only going to put in one.
#Also, going to make the limits for the "bin_E" the same in all cases.

#Now it is time for some of the graphs that Alexandre asked for. The first graph:

min_log_I <- min(data_univ[which(data_univ$P>1),]$log_I+log(data_univ[which(data_univ$P>1),]$P))*0.99
max_log_I <- max(data_univ[which(data_univ$P>1),]$log_I+log(data_univ[which(data_univ$P>1),]$P))*1.01
comp_data <- data_univ[which(data_univ$P > 1 | (data_univ$log_I >= min_log_I & data_univ$log_I<=max_log_I) ),]
print(min(comp_data$lambda))
print(max(comp_data$lambda))
#comp_data$yax <- log(comp_data$yax)
comp_data$log_I <- comp_data$log_I+log(comp_data$P)
plot <- ggplot(data=comp_data, aes(x=log_I, y=yax))#, color=interaction(scheme,P,sep=':')))#, shape=as.factor(P)))
plot <- plot + geom_point(size=1.25)
plot <- plot + scale_x_continuous(name="log(I*P)")
plot <- plot + scale_y_continuous(name="E/E_min")
plot <- plot + facet_grid(rows=vars(interaction(scheme,P,sep=':')))
#plot <- plot + scale_shape_manual(values=c(15,16,17,18))
ggsave(plot=plot, filename="hsa_cdr_serial.pdf",height=7*2, width=7)

m_1_labeller <- function(string){ return (paste("log10(m_1) =", round(as.numeric(paste(string)),2))) }
data_20_h <- data_univ[which(data_univ$P == 20 & data_univ$scheme == "hsa"),]
plot <- ggplot(data=data_20_h,aes(x=log10(m_0), y=log_I,color=yax))
plot <- plot + scale_color_gradientn(name="E/E_min",colors=c('orange','purple','green'), values=rescale(c(1,2.08/2, 1.08)), limits=c(1,1.08) )
plot <- plot + geom_point(size=0.5)
plot <- plot + scale_x_continuous(name="log10(m_0)")
plot <- plot + scale_y_continuous(name="log(I)")
plot <- plot + facet_grid(rows=vars(log10(m_1)),labeller=labeller(.rows=m_1_labeller))
ggsave(plot=plot, filename="m_0m_1vlogI.pdf", height=5, width=5*1.5)

int_frames <- list()
m_0s_l<- list()
m_1s_l<- list()
k=1
procs=c(20)

data_parallel_h <- data_univ[which("hsa"==data_univ$scheme),]
data_parallel_h$bin_E <- data_parallel_h$bin_E - min(data_parallel_h$bin_E)
data_parallel_h$log_I <- data_parallel_h$log_I + log(data_parallel_h$P)
n_E                   <- max(data_parallel_h$bin_E)
max_log_I             <- max(data_parallel_h$log_I)
min_log_I             <- min(data_parallel_h$log_I)
c_max_log_I           <- ceiling(max_log_I*K_I)/K_I
K_I_max_log_I         <- ceiling(K_I*max_log_I)
data_parallel_h$bin_I <- as.integer(K_I_max_log_I-floor(K_I*data_parallel_h$log_I))
n_I                   <- max(data_parallel_h$bin_I)
data_parallel_h$bin   <- (((data_parallel_h$bin_I+data_parallel_h$bin_E+1)*(data_parallel_h$bin_I+data_parallel_h$bin_E))/2)+data_parallel_h$bin_I

ngraphs <- ((n_E + n_I +1)*(n_E + n_I ))/2 + n_I


for (proc in procs){
  m_0s <- unique(data_parallel_h$m_0)
  m_1s <- unique(data_parallel_h$m_1)
  m_0s_l[[k]]=m_0s
  m_1s_l[[k]]=m_1s
  l_m_1s <- length(m_1s)
  int_frame <- matrix(0.0, ncol=ngraphs + 2, nrow=length(m_0s)*l_m_1s)
  for (i in 1:length(m_0s)){
    for (l in 1:l_m_1s){
      row <- (i-1)*l_m_1s + l
      m_0 <- m_0s[i]
      m_1 <- m_1s[l]
      curr_data <- data_parallel_h[which(data_parallel_h$m_0==m_0 & data_parallel_h$m_1==m_1),]
      int_frame[row,1] = m_0
      int_frame[row,2] = m_1
      for (j in 1:nrow(curr_data)){
        curr_row <- curr_data[j,]
        bin                = curr_row$bin
        int_frame[row,bin+2] = 1/nrow(curr_data)+int_frame[row,bin+2]
        }
     }
  }
  int_frames[[k]]=int_frame
  k = k+1
}
k=1

for (proc in procs){
  m_0s <- m_0s_l[[k]]
  m_1s <- m_0s_l[[k]]
  curr_int_frame <- int_frames[[k]]
  int_data       <- as.data.frame(curr_int_frame)
  colnames(int_data) <- c("m_0", "m_1", seq(1,ngraphs,by=1))
  melted_int_data    <- melt(int_data, id.vars=c("m_0","m_1"))
  colnames(melted_int_data) <- c("m_0","m_1","bin","density")
  #melted_int_data$bin_E     <- rep(0,nrow(melted_int_data))
  #melted_int_data$bin_I     <- rep(0,nrow(melted_int_data))
  melted_int_data$bin <- as.integer(paste(melted_int_data$bin))
  melted_int_data[,c("bin_E", "bin_I")] <- inverse_pairing_func(melted_int_data$bin)
  melted_int_data <- melted_int_data[which(melted_int_data$bin_E <= n_E & melted_int_data$bin_I <= n_I),]
  transform_bin_E <- function(bin){
    orig_bin <- bin+min(data_univ[which(data_univ$scheme == "hsa"),]$bin_E) -1
    return ((K_E_max_yax-orig_bin)/K_E)
  }
  transform_bin_I <- function(bin){
    return ((K_I_max_log_I - bin)/K_I)
  }
  for (m_0 in m_0s)
    for (m_1 in m_1s)
      melted_int_data[nrow(melted_int_data)+1,] <- c(m_0, m_1, 0, 0, 0, 0)
#  c_bin_E         <- min(int_data[!is.na(int_data$log_I),]$log_I)
#  m_bin_E         <- (max(int_data[!is.na(int_data$log_I),]$log_I)-c_bin_E)
#  int_data <- int_data[!is.na(int_data$log_I),]
#  transform_bin_E <- function(x){
#  x=x-min(int_data$bin_E)
#  x=x/(max(int_data$bin_E)-min(int_data$bin_E))
#  return (m_bin_E*x + c_bin_E)
#   }

 # inverse_bin_E   <- function(x){
 # x=(x-c_bin_E)/m_bin_E
 # x=x*(max(int_data$bin_E)-min(int_data$bin_E))
 # x=x+min(int_data$bin_E)
 # return(x)
 # } 

    m_0_labeller <- function(string){ return (paste("m_0 =", string)) }
    m_1_labeller <- function(string){ return (paste("m_1 =", string)) }
    plot <- ggplot(data=melted_int_data[which(melted_int_data$bin_I != 0),], aes(x=bin_I, y=bin_E,fill=density)) 
#  plot <- plot + geom_point(aes(y=transform_bin_E(bin_E), size=density, fill=density), shape=22, stroke=NA)
#  plot <- plot + scale_size_continuous(range=c(0.1, 2))
    plot <- plot + scale_fill_gradientn(colors=c("white","orange","purple"), values=rescale(c(0,0.05,max(melted_int_data$density))))
    plot <- plot + geom_tile()
    plot <- plot + coord_fixed()
    plot <- plot + scale_x_reverse(name="log(I*P)",breaks=c(2,4,6,8,10,12,14,16),expand=c(0,0),labels= lapply(c(2,4,6,8,10,12,14,16), transform_bin_I))
    plot <- plot + scale_y_reverse(name="E/E_min",breaks=c(1,3,5,7,9,11), expand=c(0,0), labels=lapply(c(1,3,5,7,9,11), transform_bin_E))
    plot <- plot + facet_grid(rows = vars(m_0), cols = vars(m_1), labeller= labeller(.rows=m_0_labeller, .cols=m_1_labeller))
    plot <- plot + theme(panel.spacing= unit(.01, "lines"),
                         panel.border = element_rect(color="black", fill=NA, linewidth=0.2))
                         #strip.background = element_rect(color="black",size=1))

  #plot <- plot + scale_color_gradientn(colors=c('orange','purple','green'), values=rescale(c(1,2.125/2, 1.125)), limits=c(1,1.125) )
    ggsave(filename="heat_map_fracscore.pdf", plot=plot, width=14, height=14)
}


