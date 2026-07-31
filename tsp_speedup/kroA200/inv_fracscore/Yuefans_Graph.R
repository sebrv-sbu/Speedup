library(ggplot2)
library(zoo)
library(dplyr)
library(reshape2)
library(scales)

dist_data_frame <- read.csv("csvs/dist_data_frame.csv", header=TRUE)
orig_data_s <- read.csv("csvs/orig_data_s.csv", header=TRUE)
tot_int_data <- read.csv("csvs/tot_int_data.csv", header=TRUE)
dist_data_frame$X <- NULL
orig_data_s$X <- NULL
tot_int_data$X <- NULL
K_E_max_yax <- 99 
K_I_max_log_I <- 460
min_comp_data_E <- -1

K_E=100
K_I=20

transform_bin_E <- function(bin){
  orig_bin <- bin+min_comp_data_E +2
  return ((K_E_max_yax - orig_bin)/K_E)
}
transform_bin_I <- function(bin){
  return ((K_I_max_log_I - bin)/K_I)
}


max_lambda <- max(dist_data_frame$min_lambda)

tot_int_data <- tot_int_data[which(tot_int_data$lambda <= max_lambda),]
max_density <- max(tot_int_data$density)

#The following is to pick the intervals we want to represent in the graph.#
graph_data <- tot_int_data[which((((tot_int_data$P == 200 & tot_int_data$m_0 == 1000 & tot_int_data$m_1 == 1000) | 
                                 (tot_int_data$P == 100 & tot_int_data$m_0 == 1000 & tot_int_data$m_1 == 1000)  )|
                                 ((tot_int_data$P == 20  & tot_int_data$m_0 == 1000   & tot_int_data$m_1 == 10) |
                                 (tot_int_data$P == 20  & tot_int_data$m_0 == 150  & tot_int_data$scheme == "cdr"))) |
                                 (((tot_int_data$P == 10  & tot_int_data$m_0 == 200 )|
                                 (tot_int_data$P == 5   & tot_int_data$m_0 == 15))|
                                 (tot_int_data$P == 1   & tot_int_data$lambda %in% c(1e-07, 2.5e-07, 5e-07, 1e-06, 2.5e-06, 5e-06)))),]
max_bin_I <- max(graph_data$bin_I)
min_bin_I <- min(graph_data$bin_I)
max_bin_E <- max(graph_data$bin_E)
min_bin_E <- min(graph_data$bin_E)
graph_data$bin_E <- graph_data$bin_E -0.4
graph_data$bin_I <- graph_data$bin_I -0.5
YuefanPlot <- ggplot(data = graph_data, aes(x=bin_I, y=bin_E, color=density, size=factor(P), shape=scheme))+geom_point(stroke=1)
YuefanPlot <- YuefanPlot + scale_shape_manual(values=c(1,5,0))
YuefanPlot <- YuefanPlot + scale_color_gradientn(colors=c("gray","orange","red", "purple", "deepskyblue"), values=rescale(c(0, 0.05, max_density/3, max_density/2, max_density)))
YuefanPlot <- YuefanPlot + scale_x_reverse(name="log(I*P)",breaks=seq(max_bin_I, min_bin_I, by=-15),expand=c(0.05,0.05),labels= lapply(seq(max_bin_I, min_bin_I, by=-15), transform_bin_I))
YuefanPlot <- YuefanPlot + scale_y_reverse(name="E_min/E",breaks=seq(max_bin_E,min_bin_E, by=-2), expand=c(0.05,0.05), labels=lapply(seq(max_bin_E, min_bin_E, by=-2), transform_bin_E))
YuefanPlot <- YuefanPlot + scale_size_manual(values=c(1,2,4,6,9,12))
#YuefanPlot <- YuefanPlot + geom_text(aes(label=ifelse(P == 1,as.character(lambda),' ')),color="black", hjust=0, vjust=0, show.legend=FALSE)

ggsave(YuefanPlot, filename="pdfs/speedup.svg")
ggsave(YuefanPlot, filename="pdfs/speedup.pdf")

new_density<-rep(0,nrow(graph_data))
for (i in (1:nrow(graph_data))){
  row <- graph_data[i,]

  new_density[i] <- sum(graph_data[which((graph_data$bin_E == row$bin_E &
                                      graph_data$P == row$P) &
                                      (graph_data$lambda==row$lambda & 
                                       graph_data$scheme == row$scheme)),]$density)
}
print(new_density)
graph_data$density <- new_density
graph_data <- graph_data[which((graph_data$P>1 | graph_data$lambda == 1e-07)),]
NewPlot <-ggplot(data=graph_data,aes(x=bin_E, y=density,
                                     color=interaction(P,scheme)))
NewPlot <- NewPlot + scale_shape_manual(values=c(1,5,0))
NewPlot <- NewPlot + geom_point(aes(size=factor(P), shape=scheme))
NewPlot <- NewPlot + scale_size_manual(values=c(1,2,4,6,8,10))
NewPlot <- NewPlot + geom_line()
NewPlot <- NewPlot + scale_x_reverse(name="E_min/E", 
                                     breaks=seq(max_bin_E,min_bin_E, by=-2),
                                     expand=c(0.05, 0.05), 
                                     labels=lapply(seq(max_bin_E,min_bin_E, by = -2),
                                                   transform_bin_E))
ggsave(NewPlot, filename="pdfs/DensityAndEnergy.pdf")
                                     
NewPlot <-ggplot(data=graph_data,aes(x=bin_E, y=density,
                                     color=factor(P), size=factor(P), shape=scheme))
NewPlot <- NewPlot + scale_shape_manual(values=c(1,5,0))
NewPlot <- NewPlot + geom_point()
NewPlot <- NewPlot + scale_x_reverse(name="E_min/E", 
                                     breaks=seq(max_bin_E,min_bin_E, by=-2),
                                     expand=c(0.05, 0.05), 
                                     labels=lapply(seq(max_bin_E,min_bin_E, by = -2),
                                                   transform_bin_E))
ggsave(NewPlot, filename="pdfs/DensityAndEnergyPoint.pdf")
 
#procs=c(5,10,20,100,200)
#schemes = c("cdr", "hsa")
#for (proc in procs){
#  for (scheme in schemes){
#    if(nrow(dist_data_frame[which(dist_data_frame$P == proc & dist_data_frame$scheme == scheme),] )>0){
#      data <- dist_data_frame[which(dist_data_frame$P == proc & dist_data_frame$scheme == scheme & dist_data_frame$min_dist),]
#      print(data)
#    }
#  }
#}  

