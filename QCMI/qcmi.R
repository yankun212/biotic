library("phyloseq")
library("ggClusterNet")
library("tidyverse")
library("WGCNA")
library("igraph")
library("Matrix")
library("vegan")
library("magrittr")
library("reticulate")
library("SpiecEasi")
library("ggstatsplot")
library("qcmi")
library("dplyr")

setwd("C:/Users/maple/Desktop/code/QCMI")

################The ecological association of diffusion limitation is classified based on geographical distance.
####result_dl= assigned_process(link_table_row=edgelist, OTUabd=otu_table, p=0.05 , data= geo,cutoff=0, method="dl")


################Filtered environmental ecological association classification, based on selected environmental differences


env<-read.table("env.txt",header = T)
otu<-read.table("otu.txt",header = T)
edge<-read.table("edge.txt",header = T)


######Detecting Environmental Redundancy
library(Hmisc)
plot(varclus(as.matrix(env) ))

result_DDL= assigned_process(link_table_row=edge, OTUabd=otu, p=0.05 , data= env['DDL'],cutoff=0, method="ef")
result_Na= assigned_process(link_table_row=edge, OTUabd=otu, p=0.05 , data= env['Na'],cutoff=0, method="ef")
result_TN= assigned_process(link_table_row=edge, OTUabd=otu, p=0.05 , data= env['TN'],cutoff=0, method="ef")
result_AP= assigned_process(link_table_row=edge, OTUabd=otu, p=0.05 , data= env['AP'],cutoff=0, method="ef")
result_K= assigned_process(link_table_row=edge, OTUabd=otu, p=0.05 , data= env['K'],cutoff=0, method="ef")
result_SOC= assigned_process(link_table_row=edge, OTUabd=otu, p=0.05 , data= env['SOC'],cutoff=0, method="ef")
result_CN= assigned_process(link_table_row=edge, OTUabd=otu, p=0.05 , data= env['CN'],cutoff=0, method="ef")
result_pH= assigned_process(link_table_row=edge, OTUabd=otu, p=0.05 , data= env['pH'],cutoff=0, method="ef")
result_SWC= assigned_process(link_table_row=edge, OTUabd=otu, p=0.05 , data= env['SWC'],cutoff=0, method="ef")
result_Acid= assigned_process(link_table_row=edge, OTUabd=otu, p=0.05 , data= env['Acid'],cutoff=0, method="ef")
result_DB= assigned_process(link_table_row=edge, OTUabd=otu, p=0.05 , data= env['DB'],cutoff=0, method="ef")
result_DT= assigned_process(link_table_row=edge, OTUabd=otu, p=0.05 , data= env['DT'],cutoff=0, method="ef")




write.csv(result_pH,"result_pH.csv")
write.csv(result_SM,"result_SM.csv")
write.csv(result_DOC,"result_DOC.csv")
write.csv(result_TN,"result_TN.csv")
write.csv(result_dl,"result_dl.csv")


total_link=row.names(edge)
ef_link=union(as.character(row.names(result_DDL)),as.character(row.names(result_Na)))
ef_link=union(ef_link,as.character(row.names(result_TN)))
ef_link=union(ef_link,as.character(row.names(result_AP)))
ef_link=union(ef_link,as.character(row.names(result_SWC)))
ef_link=union(ef_link,as.character(row.names(result_K)))
ef_link=union(ef_link,as.character(row.names(result_SOC)))
ef_link=union(ef_link,as.character(row.names(result_CN)))
ef_link=union(ef_link,as.character(row.names(result_pH)))
ef_link=union(ef_link,as.character(row.names(result_Acid)))
ef_link=union(ef_link,as.character(row.names(result_DT)))
ef_link=union(ef_link,as.character(row.names(result_DB)))

bi_link=setdiff(total_link,ef_link)

efedge=edge[ef_link,]
write.csv(efedge,"rmc_ef_edge.csv")
length(row.names(efedge))

biedge=edge[bi_link,]
write.csv(biedge,"rmc_bi_edge.csv")

ig.bi=graph_from_edgelist(as.matrix(biedge[,1:2]) , directed = FALSE)
ig.bi=set_edge_attr(ig.bi, 'weight', index = E(ig.bi),  as.numeric(biedge[,3]) )

##result_bi= qcmi(igraph= ig.bi, OTU= otu_table, pers.cutoff=0)
result_bi= qcmi(igraph= ig.bi, OTU= otu, pers.cutoff=0)


int<-data.frame(result_bi$interaction.neg,result_bi$interaction.pos)

write.csv(int,"int_rmc.csv")

