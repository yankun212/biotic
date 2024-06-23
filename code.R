######"The different biotic associations among rhizosphere microbial communities affected the tolerance of male and female poplar under salt stress"



rm(list=ls())

####Fig.2i,j  cPCoA
library(ggthemes)
library(ggplot2)
library(factoextra)
library(ggsci)
library(vegan)
library(amplicon)
library(ggConvexHull)

otu<-read.table("ASV_table.txt",header = T)#Import ASV data
group<-read.table("env.txt",header = T)#Import group data

#operation models
beta_cpcoa(otu,group,dis="bray",groupID = "group",ellipse = F)+geom_point(size=5)+
  geom_convexhull(alpha = 0.5, aes(fill = group))+theme_bw()+
  theme(text = element_text(size = 15,face="bold",family = "serif"))+
  scale_fill_npg()+scale_color_npg()

#Adonis
distance <- vegdist(t(otu), method = 'bray')
ub_dis_table <- as.dist(distance, diag = FALSE, upper = FALSE)
adonis(ub_dis_table~group$group, data=otu, permutations = 999)



####networks
library(WGCNA)
library(igraph)
library(ggplot2)

#Correlation matrix calculation
#It is necessary to remove the taxa with more 0 values in the group before analysis
occor<-corAndPvalue(t(otu),method="pearson",use="p")
occor.r = occor$cor 
occor.p = occor$p 
occor.r[occor.p>0.01|abs(occor.r)<0.8] = 0
diag(occor.r)<-0
occor.r[is.na(occor.r)]<-0
#edge
sum(abs(occor.r)>0)/2 
#node
sum(colSums(abs(occor.r))>0)
#Export the data and plot using Gephi
write.csv(occor.r,file="network.csv")
##### What is obtained here is only the co-occurrence network, and two steps of filtering are needed to obtain the microbial biological association network.


#### keystone nodes Fig.3f
library(reshape2)
edge<-read.table("clipboard",header = T) #Import the edge file obtained by the biological association network
#The edge file is converted into an association matrix
occor.r<-dcast(edge,edge$Source~edge$Target,mean)
row.names(occor.r)<-occor.r$`edge$Source`
occor.r[is.na(occor.r)]=0
occor.r<-occor.r[-1]
occor.r[abs(occor.r)>0]=1
adjacency_unweight<-occor.r

# Obtain an undirected network with no weights
igraph <- graph_from_adjacency_matrix(as.matrix(adjacency_unweight), mode = 'undirected', weighted = NULL, diag = FALSE)

#computational node
V(igraph)$degree <- degree(igraph)

set.seed(123)
V(igraph)$modularity <- membership(cluster_fast_greedy(igraph))

#integrated data
nodes_list <- data.frame(
  nodes_id = V(igraph)$name, 
  degree = V(igraph)$degree, 
  modularity = V(igraph)$modularity
)
head(nodes_list) 

#Calculate within-module connectivity (Zi) and among-module connectivity (Pi)
source('hub_node.r')

row.names(nodes_list)<-nodes_list$nodes_id
nodes_list<-nodes_list[-1]

zi_pi <- zi.pi(nodes_list, adjacency_unweight, degree = 'degree', modularity_class = 'modularity')
head(zi_pi)

zi_pi <- na.omit(zi_pi)   # remove NA
zi_pi[which(zi_pi$within_module_connectivities < 2.5 & zi_pi$among_module_connectivities < 0.62),'type'] <- 'Peripherals'
zi_pi[which(zi_pi$within_module_connectivities < 2.5 & zi_pi$among_module_connectivities > 0.62),'type'] <- 'Connectors'
zi_pi[which(zi_pi$within_module_connectivities > 2.5 & zi_pi$among_module_connectivities < 0.62),'type'] <- 'Module hubs'
zi_pi[which(zi_pi$within_module_connectivities > 2.5 & zi_pi$among_module_connectivities > 0.62),'type'] <- 'Network hubs'

ggplot(zi_pi, aes(among_module_connectivities, within_module_connectivities)) +
  geom_point(aes(color = type,shape=type), alpha = 0.5, size = 2) +
  scale_color_manual(values = c('gray','red','blue','purple'), 
                     limits = c('Peripherals', 'Connectors', 'Module hubs', 'Network hubs'))+
  theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'), 
        panel.background = element_blank(), legend.key = element_blank()) +
  labs(x = 'Among-module connectivities', y = 'Within-module connectivities', color = '') +
  geom_vline(xintercept = 0.62) +
  geom_hline(yintercept = 2.5)





####Fig.4C The contribution of environmental factors to microbial community construction was calculated.
library(rdacca.hp)
r<-rdacca.hp(t(otu),env,method="RDA",type="adjR2")
plot(r)




####Fig.5 SEM
library(piecewiseSEM)

#ModuleA
mod <- psem(
  lm(HI~Shannon_B+Acid+CD+Pos_B+Neg_B+Gender,data = a),
  lm(Acid~Na_soil+Gender,data = a),
  lm(Pos_B~Acid+Na_soil+Gender,data = a),
  lm(Neg_B~Acid+Na_soil+Gender+X2f+Pos_xj,data = a),
  lm(X2f~Acid+Na_soil+Gender,data = a),
  lm(Shannon_B~Neg_B+CDf+Na_soil,data = a),
  data = a
)
summary(mod)

#ModuleB
mod2 <- psem(
  lm(HI~Shannon_F+Acid+CD+Neg_F+Gender+Na_soil+Pos_F,data = b),
  lm(Acid~Na_soil+Gender,data = b),
  lm(Pos_F~Acid+Na_soil+Gender+CD,data = b),
  lm(Neg_F~Acid+Na_soil+Gender,data = b),
  lm(CD~Acid+Na_soil+Gender,data = b),
  lm(Shannon_F~+CD+Na_soil+Acid+Gender+Neg_F+Pos_F,data = b),
  data = b
)
summary(mod2)


####Fig.S11 TITAN
library(TITAN2)

#Import the ASV table
#ASVs with an occurrence rate of less than 3 in all samples are deleted.
otu<-read.table("clipboard",header = T)

#Import the environment table
env<-read.table("clipboard",header = T)

#run module
glades.titan2 <- titan(env, t(otu),
                       minSplt = 5, numPerm = 250, boot = TRUE, nBoot = 500, imax = FALSE,
                       ivTot = FALSE, pur.cut = 0.85, rel.cut = 0.85, ncpus = 1, memory = FALSE
)

plot_taxa_ridges(glades.titan2, ytxt.sz=8)

XJ_BM<-glades.titan2$sppmax
write.csv(XJ_BM,"XJ_BM.csv")



#####Fig. S16 Community assembly analysis, BNTI 
library(picante)
# Import ASV tables and trees
comun<-read.table("ASV_table.txt",header = T)#Import ASV data
phy<-read.tree("rooted_tree.tre")#Import tree data
comun<-t(comun)
prune_tree<-prune.sample(comun,phy)
phydist<-cophenetic(prune_tree)
#Calculate the NTI
mntd<-ses.mntd(comun,phydist,null.model="taxa.labels",abundance.weighted=T, runs=999)
#The observed βMNTD
comdist.resultS12<-comdistnt(comun,phydist)
#nullmodel
f<-function(){
  g<-randomizeMatrix(comun,null.model=c("frequency","richness","independentswap","trialswap"),iterations=1000)
  fc<-comdist.result<-comdistnt(g,phydist)
}
mt<-mean(replicate(999, f()))
mt.sd<-sd(replicate(999, f()))
BNTI=(comdist.resultS12-mt)/mt.sd