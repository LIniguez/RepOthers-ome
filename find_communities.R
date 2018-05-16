#!/usr/bin/env Rscript


args<-commandArgs(trailingOnly=T)
samp<-args[1]
round<-args[2]
maximum<-as.numeric(args[3])

infile<-paste(samp,"vertex_weight.txt",sep="")

library(igraph)
data<-read.table(infile,header=F)
data2<-data[data[,4]>0,]
data2<-data2[data2[,5]>0,]
net<-graph_from_data_frame(data[,c(1,2)],directed = F)
net2<-graph_from_data_frame(data2[,c(1,2)],directed = F)

E(net)$weight<-data[,6]

wt<-cluster_walktrap(net)

tam_WT<-sample(sizes(wt))
sum<-0
file_list<-list()
cont1<-1
cont2<-1
file_list[[cont1]]<-array()
for (i in 1:length(tam_WT)){
	sum<-sum+ tam_WT[i]
	if(sum >= maximum){
		cont1<-cont1+1
		cont2<-1
		file_list[[cont1]]<-array()
		file_list[[cont1]][cont2]<-names(tam_WT[i])
		sum<-tam_WT[i]
	}else{
		file_list[[cont1]][cont2]<-names(tam_WT[i])
		
	}
    cont2<-cont2+1
}
for (i in 1:length(file_list)){
	tempout<-paste(samp,i, "_WT_commuity_",round,".txt",sep="")
	write.table(as.character(V(net)[is.element(membership(wt), file_list[[i]])]$name),
		file=tempout,sep="\t",quote=F)
}

