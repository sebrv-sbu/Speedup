library(ggplot2)
library(zoo)
library(dplyr)
library(ggpubr)
library(reshape2)
library(scales)

data_s    <- read.table("../s_results", header=TRUE)
data_5_c  <- read.table("../5_c_results", header=TRUE)
data_10_c <- read.table("../10_c_results", header=TRUE)
data_20_c <- read.table("../20_c_results", header=TRUE)
data_20_h <- read.table("../20_h_results", header=TRUE)
data_100_h <- read.table("../100_h_results", header=TRUE)
data_200_h <- read.table("../200_h_results", header=TRUE)

orig_data_s <- data_s

data_s$P      = rep(1, nrow(data_s))
data_5_c$P    = rep(5, nrow(data_5_c))
data_10_c$P   = rep(10,nrow(data_10_c))
data_20_c$P   = rep(20,nrow(data_20_c))
data_20_h$P   = rep(20,nrow(data_20_h))
data_100_h$P  = rep(100,nrow(data_100_h))
data_200_h$P  = rep(200,nrow(data_200_h))

data_parallel_c <- bind_rows(data_5_c,data_10_c,data_20_c)
names(data_parallel_c)[names(data_parallel_c) == "Interval"] <- "m_0"
data_parallel_c$lambda = rep(0.000005, nrow(data_parallel_c))
data_s$scheme=rep("serial", nrow(data_s))
data_parallel_c$scheme=rep("cdr", nrow(data_parallel_c))
data_parallel_h <- bind_rows(data_20_h, data_100_h, data_200_h)
data_parallel_h$scheme <- rep("hsa",nrow(data_parallel_h))
data_parallel_h$lambda <- rep(0.000005,nrow(data_parallel_h))

data_univ<-bind_rows(data_s,data_parallel_c, data_parallel_h)
E_min=14382.995939

K_E=100
#changing the "score" value here.
data_univ[which(data_univ$scheme != "hsa"),]$I = data_univ[which(data_univ$scheme != "hsa"),]$I-100000 - 100000/data_univ[which(data_univ$scheme != "hsa"),]$P
data_univ[which(data_univ$scheme == "hsa" & data_univ$P==20),]$I = data_univ[which(data_univ$scheme =="hsa" & data_univ$P==20),]$I-100000
-100000/4
data_univ[which(data_univ$scheme == "hsa" & data_univ$P==100),]$I = data_univ[which(data_univ$scheme =="hsa" & data_univ$P==100),]$I-100000
-100000/20
data_univ[which(data_univ$scheme == "hsa" & data_univ$P==200),]$I = data_univ[which(data_univ$scheme =="hsa" & data_univ$P==200),]$I-100000
-100000/20


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

lambda_int_model_for_print <- lm(data_s[which(data_s$lambda<0.01),]$log_I ~ log(data_s[which(data_s$lambda <0.01),]$lambda) )#+ I(log(data_s[which(data_s$lambda <0.01),]$lambda**2)))
lambda_int_model_for_print_complete <- lm(data_s$log_I ~ log(data_s$lambda) )
print("Lambda_Int_model_complete:")
summary(lambda_int_model_for_print_complete)
lambda_quad_model_for_print <- lm(data_s$log_I ~ log(data_s$lambda) + I(log(data_s$lambda)**2))
summary(lambda_quad_model_for_print)
lambda_int_model <- lm(data_s[which(data_s$lambda<0.01),]$log_I ~ log10(data_s[which(data_s$lambda <0.01),]$lambda)) #+ I(log10(data_s[which(data_s$lambda <0.01),]$lambda**2)))
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
serial_data_for_later <- lambda_data[,seq.int(ncol(lambda_data)/2)]
colnames(serial_data_for_later) <- seq(1, ngraphs, by =1)
serial_data_for_later$lambda = lambdas
serial_data_for_later <- melt(serial_data_for_later, id=c("lambda"))
colnames(serial_data_for_later) <- c("lambda", "bin", "density")
serial_data_for_later$density <- as.numeric(paste(serial_data_for_later$density))
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
plot <- plot + ggtitle("105 City Test Problem Serial Results")
ggsave(filename="pdfs/lambdas_iterations_serial_fracscore.pdf", plot=plot,height=5, width=5*1.5)
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
ggsave(filename="pdfs/zoomed_lambdas_iterations_serial_fracscore.pdf", plot=plot,height=5, width=5*1.5)

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

comp_data <- data_univ[which(data_univ$P > 1 | (data_univ$lambda < 1e-4)),]
subtract_bin <- min(comp_data$bin_E)
comp_data$bin_E <- comp_data$bin_E - subtract_bin
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
procs = c(20,100, 200)
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
      curr_data <- data[which(data$m_0==m_0 & data$m_1==m_1),]
      int_frame[row,1] = 5e-6
      int_frame[row,2] = m_0
      int_frame[row,3] = m_1
      for (j in 1:nrow(curr_data)){
        curr_row <- curr_data[j,]
        bin                = curr_row$bin
        int_frame[row,bin+3] = 1/nrow(curr_data)+int_frame[row,bin+3]
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
colnames(tot_int_data) <- columns
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
plot_data <-tot_int_data[which(tot_int_data$scheme == "cdr" & tot_int_data$m_0 %in% ints),]
comp_data <- plot_data[which(plot_data$density>0),]
plot_data <- plot_data[which(plot_data$bin_I <= max(comp_data$bin_I) & plot_data$bin_I >= min(comp_data$bin_I) 
                             & plot_data$bin_E >= min(comp_data$bin_E) & plot_data$bin_E <= max(comp_data$bin_E)),]
plot <- ggplot(data=plot_data, aes(x=bin_I, y=bin_E, fill=density))
plot <- plot + geom_tile()
plot <- plot + facet_grid(rows=vars(P), cols=vars(m_0))
plot <- plot + scale_x_reverse(name="log(I*P)",breaks=c(25, 33, 41),expand=c(0,0),labels= lapply(c(25,33, 41), transform_bin_I))
plot <- plot + scale_y_reverse(name="E/E_min",breaks=c(1,5,9,13,17,21), expand=c(0,0), labels=lapply(c(1, 5, 9, 13,17,21), transform_bin_E))
plot <- plot + scale_fill_gradientn(colors=c("white","orange","purple"), values=rescale(c(0,0.05,max(plot_data$density))))
plot <- plot + theme(panel.spacing= unit(.01, "lines"),
                         panel.border = element_rect(color="black", fill=NA, linewidth=0.2))
plot <- plot + ggtitle("105 City Test Run CDR Heatmap")
plot <- plot + annotate("text", x=Inf, y=-Inf, label="Mixing interval (m₀)", hjust=1.1, vjust=-2, angle=0, size=4, fontface="italic") 
plot <- plot +   annotate("text", x=-Inf, y=Inf, label="Number of processors (P)",             hjust=-0.3, vjust=1.3, angle=90, size=4, fontface="italic")
#ggsave(filename="pdfs/cdr_heatmap.pdf", plot=plot, height=5, width=5*2)
#ggsave(filename="pdfs/cdr_heatmap.svg", plot=plot, height=5, width=10)

data_c <- data_univ[which(data_univ$scheme == "cdr"),]
plot <- ggplot(data=tot_int_data[which(tot_int_data$density > 0),], aes(x=bin_I, y=bin_E, color=density))
plot <- plot + geom_point()
plot <- plot + facet_grid(rows=vars(P))
plot <- plot + scale_color_gradientn(name="density", colors=c("#66cc00","#cc3399", "#43123c"))
ggsave(filename="pdfs/int_energy_cdr_fracscore.pdf", plot=plot)
#The above graphs look like shit when we put both factors in, so we are only going to put in one.
#Also, going to make the limits for the "bin_E" the same in all cases.

#Now it is time for some of the graphs that Alexandre asked for. The first graph:

min_log_I <- min(data_univ[which(data_univ$P>1),]$log_I+log(data_univ[which(data_univ$P>1),]$P))*0.99
max_log_I <- max(data_univ[which(data_univ$P>1),]$log_I+log(data_univ[which(data_univ$P>1),]$P))*1.01
comp_data <- data_univ[which(data_univ$P > 1 | (data_univ$log_I >= min_log_I & data_univ$log_I<=max_log_I) ),]
#comp_data$yax <- log(comp_data$yax)
comp_data$log_I <- comp_data$log_I+log(comp_data$P)
plot <- ggplot(data=comp_data, aes(x=log_I, y=yax))#, color=interaction(scheme,P,sep=':')))#, shape=as.factor(P)))
plot <- plot + geom_point(size=1.25)
plot <- plot + scale_x_continuous(name="log(I*P)")
plot <- plot + scale_y_continuous(name="E/E_min")
plot <- plot + facet_grid(rows=vars(interaction(scheme,P,sep=':')))
#plot <- plot + scale_shape_manual(values=c(15,16,17,18))
ggsave(plot=plot, filename="pdfs/hsa_cdr_serial.pdf",height=7*2, width=7)

m_1_labeller <- function(string){ return (paste("m_1 =", string)) }
P_labeller <- function(string) { return (paste("P = ", string)) }
data_20_h <- data_univ[which(data_univ$P == 20 & data_univ$scheme == "hsa"),]
plot <- ggplot(data=data_20_h,aes(x=log10(m_0), y=log_I,color=yax))
plot <- plot + scale_color_gradientn(name="E/E_min",colors=c('orange','purple','green'), values=rescale(c(1,2.08/2, 1.08)), limits=c(1,1.08) )
plot <- plot + geom_point(size=0.5)
plot <- plot + scale_x_continuous(name="log10(m_0)")
plot <- plot + scale_y_continuous(name="log(I)")
plot <- plot + facet_grid(rows=vars(m_1),labeller=labeller(.rows=m_1_labeller))
ggsave(plot=plot, filename="pdfs/m_0m_1vlogI.pdf", height=5, width=5*1.5)

int_frames <- list()
m_0s_l<- list()
m_1s_l<- list()
k=1
procs=c(20,100, 200)

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
  proc_data <- data_parallel_h[which(data_parallel_h$P == proc),]
  m_0s <- unique(proc_data$m_0)
  m_1s <- unique(proc_data$m_1)
  m_0s_l[[k]]=m_0s
  m_1s_l[[k]]=m_1s
  l_m_1s <- length(m_1s)
  int_frame <- matrix(0.0, ncol=ngraphs + 2, nrow=length(m_0s)*l_m_1s)
  for (i in 1:length(m_0s)){
    for (l in 1:l_m_1s){
      row <- (i-1)*l_m_1s + l
      m_0 <- m_0s[i]
      m_1 <- m_1s[l]
      curr_data <- proc_data[which(proc_data$m_0==m_0 & proc_data$m_1==m_1),]
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
    plot <- ggplot(data=melted_int_data, aes(x=bin_I, y=bin_E,fill=density)) 
#  plot <- plot + geom_point(aes(y=transform_bin_E(bin_E), size=density, fill=density), shape=22, stroke=NA)
#  plot <- plot + scale_size_continuous(range=c(0.1, 2))
    plot <- plot + scale_fill_gradientn(colors=c("white","orange","purple"), values=rescale(c(0,0.05,max(melted_int_data$density))))
    plot <- plot + geom_tile()
    plot <- plot + coord_fixed()
    plot <- plot + scale_x_reverse(name="log(I*P)",breaks=seq(1, max(melted_int_data$bin_I), by=4),expand=c(0,0),labels= lapply(seq(1, max(melted_int_data$bin_I), by=4), transform_bin_I))
    plot <- plot + scale_y_reverse(name="E/E_min",breaks=seq(0,max(melted_int_data$bin_E), by=2), expand=c(0,0), labels=lapply(seq(0,max(melted_int_data$bin_E), by=2), transform_bin_E))
    plot <- plot + facet_grid(rows = vars(m_0), cols = vars(m_1), labeller= labeller(.rows=m_0_labeller, .cols=m_1_labeller))
    plot <- plot + theme(panel.spacing= unit(.01, "lines"),
                         panel.border = element_rect(color="black", fill=NA, linewidth=0.2))
                         #strip.background = element_rect(color="black",size=1))

  #plot <- plot + scale_color_gradientn(colors=c('orange','purple','green'), values=rescale(c(1,2.125/2, 1.125)), limits=c(1,1.125) )
  plot<-plot + ggtitle(paste(paste("heat map", proc, sep=" "), "hCDR 105 test run", sep=" "))
    #ggsave(filename=paste(paste("pdfs/heat_map_fracscore_", proc, sep=""),".pdf", sep=""), plot=plot, width=14, height=14)
    #ggsave(filename=paste(paste("pdfs/heat_map_fracscore_", proc, sep=""),".svg", sep=""), plot=plot, width=14, height=14)

    k = k + 1
}

#Now, we want to make some new graphs - CDR versus Serial in particular. In order to do this, we are going to only compare
#Energy Densities.
bin_Es <- unique(tot_int_data$bin_E)
intervals_5  <- unique(tot_int_data[which(tot_int_data$scheme == "cdr" & 5 == tot_int_data$P),]$m_0)
intervals_10 <- unique(tot_int_data[which(tot_int_data$scheme == "cdr" & 10 == tot_int_data$P),]$m_0)
intervals_20 <- unique(tot_int_data[which(tot_int_data$scheme == "cdr" & 20 == tot_int_data$P),]$m_0)
m_0s_20 <- unique(tot_int_data[which("hsa" == tot_int_data$scheme & 20 == tot_int_data$P),]$m_0)
m_1s_20 <- unique(tot_int_data[which("hsa" == tot_int_data$scheme & 20 == tot_int_data$P),]$m_1)



density_sum <- function(bin_E, m_0, m_1, P, scheme){
  if ("hsa" == scheme){
    return ( sum(tot_int_data[which(bin_E == tot_int_data$bin_E & m_0 == tot_int_data$m_0 & m_1 == tot_int_data$m_1 & P == tot_int_data$P & "hsa" == tot_int_data$scheme),]$density))
  }
  else
    return(sum(tot_int_data[which(bin_E == tot_int_data$bin_E & m_0 == tot_int_data$m_0 & P == tot_int_data$P & "cdr"==tot_int_data$scheme),]$density))

}


new_data_frame_5 <- data.frame(expand.grid(intervals_5, bin_Es))
colnames(new_data_frame_5) <- c("m_0", "bin_E")
new_data_frame_5$P <- rep(5, nrow(new_data_frame_5))
new_data_frame_10 <- data.frame(expand.grid(intervals_10, bin_Es))
colnames(new_data_frame_10) <- c("m_0", "bin_E")
new_data_frame_10$P <- rep(10, nrow(new_data_frame_10))
new_data_frame_20 <- data.frame(expand.grid(intervals_20, bin_Es))
colnames(new_data_frame_20) <- c("m_0", "bin_E")
new_data_frame_20$P <- rep(20, nrow(new_data_frame_20))
new_data_frame <-rbind(new_data_frame_10, new_data_frame_20, new_data_frame_5)
new_data_frame$m_1 <- rep(NA, nrow(new_data_frame))
new_data_frame$scheme <- rep("cdr", nrow(new_data_frame))

#temp_200_matrix <- matrix(c(1000, 1000, 100, 1000, 10, 1000, 50, 1000, 10, 10), ncol=2,byrow=TRUE)


m.combs <- c()

expand.grid.2 <- function(lst) {

  if (is.null(m.combs)) {

    m.combs <<- lst

  } else {

    m.current <- m.combs
    n <- nrow(m.combs)

    for (i in 1:nrow(lst)) {

      if(i == 1)  # for first iteration cbind new matrix
        m.combs <<- cbind(m.combs, matrix(rep(lst[i, ], each = n), nrow = n))
      else        # for next iterations rbind new matrix
        m.combs <<- rbind(m.combs, cbind(m.current, matrix(rep(lst[i, ], each = n), nrow = n)))
    }
  }
}

new_data_frame_20_h <- data.frame(expand.grid(m_0s_20, m_1s_20, bin_Es))
colnames(new_data_frame_20_h) <- c("m_0", "m_1", "bin_E")
new_data_frame_20_h$P <- rep(20, nrow(new_data_frame_20_h))
new_data_frame_20_h$scheme <- rep("hsa", nrow(new_data_frame_20_h))
new_data_frame_100_h <- data.frame(expand.grid(m_0s_20, m_1s_20, bin_Es))
colnames(new_data_frame_100_h) <- c("m_0", "m_1", "bin_E")
new_data_frame_100_h$P <- rep(100, nrow(new_data_frame_100_h))
new_data_frame_100_h$scheme <- rep("hsa", nrow(new_data_frame_100_h))
#lst <- list(temp_200_matrix, matrix(bin_Es, ncol=1))
#rapply(lst, expand.grid.2)
#new_data_frame_200_h <- data.frame(m.combs)
new_data_frame_200_h <- data.frame(expand.grid(m_0s_20, m_1s_20, bin_Es))
colnames(new_data_frame_200_h) <- c("m_0", "m_1", "bin_E")
new_data_frame_200_h$P <- rep(200, nrow(new_data_frame_200_h))
new_data_frame_200_h$scheme <- rep("hsa", nrow(new_data_frame_200_h))
new_data_frame <- rbind(new_data_frame, new_data_frame_20_h, new_data_frame_100_h, new_data_frame_200_h)

new_data_frame$bin_E <- as.integer(paste(new_data_frame$bin_E))
new_data_frame$m_0 <- as.integer(paste(new_data_frame$m_0))
new_data_frame$m_1 <- as.integer(paste(new_data_frame$m_1))
new_data_frame$density <- mapply(density_sum, new_data_frame$bin_E, new_data_frame$m_0, new_data_frame$m_1, new_data_frame$P, new_data_frame$scheme)
new_data_frame$bin_E = new_data_frame$bin_E + subtract_bin 
min_bin_E <- min(new_data_frame$bin_E)
max_bin_E <- max(new_data_frame$bin_E)

hellinger_distance_cdr <- function(interval, P, lambda){
  sum=0

  parallel_data<- new_data_frame[which( interval == new_data_frame$m_0 & P == new_data_frame$P & "cdr" == new_data_frame$scheme),]
  serial_data <- serial_data_for_later[which(lambda==serial_data_for_later$lambda),]

  for (bin in (min_bin_E : max_bin_E)){
    row_s <- serial_data[which(bin == serial_data$bin),]
    row_p <- parallel_data[which(bin == parallel_data$bin_E),]
    sum = sum + (sqrt(row_s$density)-sqrt(row_p$density))**2
  }
  for (bin in (1:(min_bin_E -1))){
    row_s <- serial_data[which(bin == serial_data$bin),]
    sum = sum + row_s$density
  }
  return(sum/sqrt(2))

}

hellinger_distance_hsa <- function(m_0, m_1, P, lambda){
  sum = 0

  parallel_data<- new_data_frame[which( m_0 == new_data_frame$m_0& m_1 == new_data_frame$m_1 & P == new_data_frame$P & "hsa" == new_data_frame$scheme ),]
  serial_data <- serial_data_for_later[which(lambda==serial_data_for_later$lambda),]

  for (bin in (min_bin_E : max_bin_E)){
    row_s <- serial_data[which(bin == serial_data$bin),]
    row_p <- parallel_data[which(bin == parallel_data$bin_E),]
    sum = sum + (sqrt(row_s$density)-sqrt(row_p$density))**2
  }
  for (bin in (1:(min_bin_E -1))){
    row_s <- serial_data[which(bin == serial_data$bin),]
    sum = sum + row_s$density
  }
  return(sum/sqrt(2))

}

find_min_distribution <- function(m_0, m_1, P, scheme){
  min_lambda=1
  min_dist=1
  if ("cdr" == scheme){
    for (lambda in lambdas){
      dist = hellinger_distance_cdr(m_0, P, lambda) 
      if (dist <= min_dist){
        min_lambda = lambda 
        min_dist = dist
      }
    }
  }
  else {
    for (lambda in lambdas){
      dist = hellinger_distance_hsa(m_0, m_1, P, lambda) 
      if (dist <= min_dist){
        min_lambda = lambda 
        min_dist = dist
      }
    }
  }
  return (c(min_lambda, min_dist))
}
lambdas <- unique(serial_data_for_later$lambda)
   
dist_data_frame_5 <- data.frame(m_0=intervals_5)
dist_data_frame_5$P <- rep(5, nrow(dist_data_frame_5))
dist_data_frame_10 <- data.frame(m_0=intervals_10)
dist_data_frame_10$P <- rep(10, nrow(dist_data_frame_10))
dist_data_frame_20 <- data.frame(m_0=intervals_20)
dist_data_frame_20$P <- rep(20, nrow(dist_data_frame_20))
dist_data_frame <- rbind(dist_data_frame_5, dist_data_frame_10, dist_data_frame_20)
dist_data_frame$scheme <- rep("cdr", nrow(dist_data_frame))
dist_data_frame$m_1 <- rep(NA, nrow(dist_data_frame))
h_data_20<- tot_int_data[which("hsa" == tot_int_data$scheme & 20 == tot_int_data$P),]
m_0sm_1s20 <- expand.grid(unique(h_data_20$m_0), unique(h_data_20$m_1))
dist_data_frame_20_h <- data.frame(m_0sm_1s20)
colnames(dist_data_frame_20_h)<-c("m_0","m_1")
dist_data_frame_20_h$P <- rep(20, nrow(dist_data_frame_20_h))
dist_data_frame_20_h$scheme <- rep("hsa", nrow(dist_data_frame_20_h))
dist_data_frame_100_h <- data.frame(m_0sm_1s20)
colnames(dist_data_frame_100_h)<-c("m_0","m_1")
dist_data_frame_100_h$P <- rep(100, nrow(dist_data_frame_100_h))
dist_data_frame_100_h$scheme <- rep("hsa", nrow(dist_data_frame_100_h))
dist_data_frame_200_h <- data.frame(m_0sm_1s20)
colnames(dist_data_frame_200_h) <- c("m_0", "m_1")
dist_data_frame_200_h$P <- rep(200, nrow(dist_data_frame_200_h))
dist_data_frame_200_h$scheme <- rep("hsa", nrow(dist_data_frame_200_h))

dist_data_frame <- rbind(dist_data_frame, dist_data_frame_20_h, dist_data_frame_100_h, dist_data_frame_200_h)

dist_data_frame[,c("min_lambda", "min_dist")] <- t(mapply(find_min_distribution, dist_data_frame$m_0, dist_data_frame$m_1, dist_data_frame$P, dist_data_frame$scheme))


avg_I_s <- function(lambda){
  s_lambda_data <- orig_data_s[which(orig_data_s$lambda == lambda),]
  return (sum(s_lambda_data$I)/nrow(s_lambda_data))
}

avg_I_p <- function(m_0, m_1, P, scheme){
  if ("cdr" == scheme){
    p_interval_data <- data_parallel_c[which(data_parallel_c$P == P & data_parallel_c$m_0 == m_0),]
    return(P*sum(p_interval_data$I)/nrow(p_interval_data))
  }
  else{
    p_interval_data <- data_parallel_h[which(data_parallel_h$P == P & data_parallel_h$m_0 == m_0 & data_parallel_h$m_1 == m_1),]
    return(P*sum(p_interval_data$I)/nrow(p_interval_data))
  }
}

dist_data_frame$I_s <- as.numeric(paste(lapply(dist_data_frame$min_lambda, avg_I_s)))
dist_data_frame$PI_p <- as.numeric(paste(mapply(avg_I_p,dist_data_frame$m_0, dist_data_frame$m_1, dist_data_frame$P, dist_data_frame$scheme)))
dist_data_frame$approx_efficiency <- dist_data_frame$I_s/dist_data_frame$PI_p

plot1 <- ggplot(data=dist_data_frame[which(dist_data_frame$scheme == "cdr"),],aes(x=log10(m_0), y=approx_efficiency, color=min_dist))
plot1 <- plot1 + geom_point()
plot1 <- plot1 + facet_grid(rows=vars(P), labeller=labeller(.rows=P_labeller))
plot1 <- plot1 + scale_y_continuous(name="approximate parallel efficiency")
plot1 <- plot1 + scale_x_continuous(name="log_10(m_0)")
plot1 <- plot1 + scale_color_gradientn(name="distance",colors=c('green', 'red','black'), values=rescale(c(0,0.2, 0.3)), limits=c(0,0.3))
plot1 <- plot1 + ggtitle("CDR lin105")
ggsave(plot=plot1, filename="pdfs/cdr_speedup.pdf", width=5, height=5*1.4)

plot <- ggplot(data=dist_data_frame[which("hsa" == dist_data_frame$scheme  & 20 == dist_data_frame$P),],aes(x=log10(m_0), y=approx_efficiency, color=min_dist))
plot <- plot + geom_point()
plot <- plot + facet_grid(rows=vars(m_1),labeller=labeller(.rows=m_1_labeller))
plot <- plot + scale_y_continuous(name = "approximate parallel efficiency", limits=c(0,2))
plot <- plot + scale_x_continuous(name = "log_10(m_0)")
plot <- plot + scale_color_gradientn(name = "distance", colors=c('green', 'red', 'black'), values=rescale(c(0,0.2,0.5)),limits=c(0,0.5))
ggsave(plot=plot, filename="pdfs/hsa_speedup.pdf", width = 5, height = 5*1.4)

plot2 <- ggplot(data=dist_data_frame[which("hsa" == dist_data_frame$scheme  & 20 == dist_data_frame$P & dist_data_frame$min_dist<0.3),],aes(x=log10(m_0), y=approx_efficiency, color=min_dist))
plot2 <- plot2 + geom_point()
plot2 <- plot2 + facet_grid(rows=vars(m_1),labeller=labeller(.rows=m_1_labeller))
plot2 <- plot2 + scale_y_continuous(name = "approximate parallel efficiency", limits=c(0,2))
plot2 <- plot2 + scale_x_continuous(name = "log_10(m_0)")
plot2 <- plot2 + scale_color_gradientn(name = "distance", colors=c('green', 'red', 'black'), values=rescale(c(0,0.2,0.3)),limits=c(0,0.3))
plot2 <- plot2 + ggtitle("HSA lin105 20proc")
ggsave(plot=plot2, filename="pdfs/hsa_speedup_20.pdf", width = 5, height = 5*1.4)
plot3 <- ggplot(data=dist_data_frame[which("hsa" == dist_data_frame$scheme & 100 == dist_data_frame$P & dist_data_frame$min_dist <0.3),],
                aes(x=log10(m_0), y=approx_efficiency, color=min_dist))
plot3 <- plot3 + geom_point()
plot3 <- plot3 + facet_grid(rows=vars(m_1),labeller=labeller(.rows=m_1_labeller))
plot3 <- plot3 + scale_y_continuous(name = "approximate parallel efficiency", limits=c(0,2))
plot3 <- plot3 + scale_x_continuous(name = "log_10(m_0)")
plot3 <- plot3 + scale_color_gradientn(name = "distance", colors=c('green', 'red', 'black'), values=rescale(c(0,0.2,0.3)),limits=c(0,0.3))
plot3 <- plot3 + ggtitle("HSA lin105 100proc")
ggsave(plot=plot3, filename="pdfs/hsa_speedup_100.pdf", width = 5, height = 5*1.4)

plot5 <- ggplot(data=dist_data_frame[which("hsa" == dist_data_frame$scheme & 200 == dist_data_frame$P),],
                aes(x=log10(m_0), y=approx_efficiency, color=min_dist))
plot5 <- plot5 + geom_point()
plot5 <- plot5 + facet_grid(rows=vars(m_1),labeller=labeller(.rows=m_1_labeller))
plot5 <- plot5 + scale_y_continuous(name = "approximate parallel efficiency", limits=c(0,2))
plot5 <- plot5 + scale_x_continuous(name = "log_10(m_0)")
plot5 <- plot5 + scale_color_gradientn(name = "distance", colors=c('green', 'red', 'black'), values=rescale(c(0,0.2,0.3)),limits=c(0,0.3))
plot5 <- plot5 + ggtitle("HSA lin105 200proc")
ggsave(plot=plot5, filename="pdfs/hsa_speedup_200.pdf", width = 5, height = 5*1.4)



plot4 <- ggarrange(plot1,plot2,ncol=2)

ggsave(plot=plot4, filename="pdfs/hsa_cdr_speedup.pdf", width=10, height=7)

data_for_graph <- data.frame(E=orig_data_s[which(orig_data_s$lambda==5e-6),]$E)
data_for_graph$scheme <- rep(1, nrow(data_for_graph))
data_for_graph_temp <- data.frame(E=data_parallel_c[which(20 == data_parallel_c$P & 10 == data_parallel_c$m_0),]$E)
data_for_graph_temp$scheme <- rep(3, nrow(data_for_graph_temp))
data_for_graph <- rbind(data_for_graph, data_for_graph_temp)

plot1 <- ggplot(data=data_for_graph, aes(x=scheme, y=E/E_min))
plot1 <- plot1 + geom_jitter()
plot1 <- plot1 + scale_x_continuous(breaks=c(1,3), labels=c("serial", "cdr"))
plot1 <- plot1 + geom_vline(xintercept=2)
plot1 <- plot1 + scale_y_continuous(name="E/E_min", breaks=seq(1,1.15,by=0.02), limits=c(0.985,1.105), expand=c(0,0))

data_for_graph2 <- data.frame(bin_E=new_data_frame[which( 10 == new_data_frame$m_0 & 20 == new_data_frame$P & "cdr" == new_data_frame$scheme),]$bin_E)
data_for_graph2$bin_E = data_for_graph2$bin_E
data_for_graph2$density <-new_data_frame[which( 10 == new_data_frame$m_0 & 20 == new_data_frame$P & "cdr" == new_data_frame$scheme),]$density
data_for_graph2$scheme <- rep(3, nrow(data_for_graph2))
data_for_graph2_temp <- data.frame(bin_E=serial_data_for_later[which(5e-6 == serial_data_for_later$lambda),]$bin)
data_for_graph2_temp$density <- serial_data_for_later[which(5e-6 == serial_data_for_later$lambda),]$density
data_for_graph2_temp$scheme <- rep(1, nrow(data_for_graph2_temp))
data_for_graph2 <- rbind(data_for_graph2, data_for_graph2_temp)
transform_bin_E <- function(bin){
    return(((K_E_max_yax - bin+1)/K_E) )
    }
data_for_graph2$bin_E <- as.numeric(paste(data_for_graph2$bin_E))
data_for_graph2 <- data_for_graph2[which(transform_bin_E(data_for_graph2$bin_E)<1.11),]
data_for_graph2_spoof <- data.frame(bin_E = seq(max_bin_E,as.integer(K_E_max_yax - floor(K_E*1.101)),by = -1))
data_for_graph2_spoof$scheme <- rep(2, nrow(data_for_graph2_spoof))
data_for_graph2_spoof$density <- rep(0, nrow(data_for_graph2_spoof))
data_for_graph2 <- rbind(data_for_graph2, data_for_graph2_spoof)
#data_for_graph2[nrow(data_for_graph2+1),]=c(163,0,2)
#data_for_graph2[nrow(data_for_graph2+1),]=c(163,0,3)

plot2 <- ggplot(data=data_for_graph2, aes(x=scheme, y=bin_E, fill = density))
plot2 <- plot2 + geom_tile()
plot2 <- plot2 + scale_fill_gradientn(colors=c("white","orange","purple"), values=rescale(c(0,0.05,max(data_for_graph2$density))))
plot2 <- plot2 + scale_y_reverse(name="",expand=c(0,0),breaks=c(),limits=c(163,as.integer(K_E_max_yax - floor(K_E*1.101))))
plot2 <- plot2 + scale_x_continuous(breaks = c(1,3), labels=c("serial", "cdr"), expand=c(0,0))

plot3 <- ggarrange(plot1, plot2, ncol=2)
ggsave(plot=plot3, filename="pdfs/bin_tech.pdf", width=8, height=5)
data_for_graph4 <- tot_int_data[which((tot_int_data$P==200 & tot_int_data$m_0 == 1000 & tot_int_data$m_1 == 1000) | 
                                   (tot_int_data$P == 5 & tot_int_data$m_0 == 100) |
                                   (tot_int_data$P == 1 & tot_int_data$lambda<1.2e-07)) ,]
max_bin_E <- max(data_for_graph4[which(data_for_graph4$density>0),]$bin_E)
min_bin_E <- min(data_for_graph4[which(data_for_graph4$density>0),]$bin_E)
max_bin_I <- max(data_for_graph4[which(data_for_graph4$density>0),]$bin_I)
min_bin_I <- min(data_for_graph4[which(data_for_graph4$density>0),]$bin_I)

data_for_graph4 <- data_for_graph4[which(data_for_graph4$bin_E<= max_bin_E & data_for_graph4$bin_E >= min_bin_E
                                         &data_for_graph4$bin_I<= max_bin_I & data_for_graph4$bin_I >= min_bin_I),]
 

transform_bin_E <- function(bin){
    orig_bin <- bin+min(comp_data$bin_E) -1
    return ((K_E_max_yax-orig_bin)/K_E)
  }
  transform_bin_I <- function(bin){
    return ((K_I_max_log_I - bin)/K_I)
  }

  print("min comp data")
  print(min(comp_data$bin_E))
  print("K_E_max_yax")
  print(K_E_max_yax)
  print("K_I_max_log_I")
  print(K_I_max_log_I)

plot3 <- ggplot(data=data_for_graph4, aes(x=bin_I, y=bin_E, fill=density))
plot3 <- plot3 + geom_tile()
plot3 <- plot3 + facet_grid(cols=vars(scheme))
plot3 <- plot3 + scale_x_reverse(name="log(I*P)",breaks=seq(max_bin_I, min_bin_I, by=-15),expand=c(0,0),labels= lapply(seq(max_bin_I, min_bin_I, by=-15), transform_bin_I))
plot3 <- plot3 + scale_y_reverse(name="E/E_min",breaks=seq(max_bin_E,min_bin_E, by=-2), expand=c(0,0), labels=lapply(seq(max_bin_E, min_bin_E, by=-2), transform_bin_E))
plot3 <- plot3 + scale_fill_gradientn(colors=c("white","orange","purple"), values=rescale(c(0,0.05,max(data_for_graph4$density))))

ggsave(plot=plot3, filename="pdfs/hsa_cdr_serial_heatmaps.pdf", width=10, height=5)
write.csv(tot_int_data[which(tot_int_data$density>0),], "csvs/tot_int_data.csv")
write.csv(orig_data_s, "csvs/orig_data_s.csv")
write.csv(dist_data_frame, "csvs/dist_data_frame.csv")
write.csv(data_univ, "csvs/data_univ.csv")
