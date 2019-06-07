#!/usr/bin/env Rscript


args<-commandArgs(trailingOnly=T)
samp<-args[1]
round<-args[2]
maximum<-as.numeric(args[3])
rand_net<-as.numeric(args[4])
sel_seed<-args[5]
nproc<-args[6]


infile<-paste(samp,"vertex_weight.txt",sep="")

library(igraph)
library("BiocParallel")

data<-read.table(infile,header=F)
net<-graph_from_data_frame(data[,c(1,2)],directed = F)

E(net)$weight<-data[,6]

set.seed(sel_seed)
index_rand<-sample(1:length(V(net)))



rntw<-function(I, INDX, RN, NTWRK){
  ini<-((I-1)*RN)+1
  fin<-ini+RN-1
  temp<-INDX[ini:fin]
  induced.subgraph(NTWRK, temp[!is.na(temp)])
}

ranNetw<-lapply(1:(as.integer(length(V(net))/rand_net)+1), rntw,INDX=index_rand, RN=rand_net, NTWRK=net )
rm(net)
wt<-bplapply(ranNetw, cluster_walktrap,  BPPARAM=MulticoreParam(workers=nproc))


file_list<-list()
cont1<-1
file_list[[cont1]]<-character()

y<-sapply(wt, function(x,m=maximum){
  tam_WT<-sample(sizes(x))
  j<-1
  for (i in 1:length(tam_WT)){
    if(sum(tam_WT[j:i]) >=m){
        file_list[[cont1]]<<-c(file_list[[cont1]],as.character(names(membership(x)[is.element(membership(x), names(tam_WT[j:i]))])))
	cont1<<-cont1+1
	file_list[[cont1]]<<-character()
	j<-i
     }
  }
  file_list[[cont1]]<<-c(as.character(names(membership(x)[is.element(membership(x), names(tam_WT[j:length(tam_WT)]))])))
  invisible(0)
})

for (i in 1:length(file_list)){
	tempout<-paste(samp,i, "_WT_commuity_",round,".txt",sep="")
	write.table(file_list[[i]][!is.na(file_list[[i]])],
		file=tempout,sep="\t",quote=F)
}
