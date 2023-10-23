
#get directories used in this script
source("~/Documents/Projects/HNSCC_PreCog/data/Pecog_HNSCC_directories.R")
RscriptsPath="~/Documents/scripts/rscripts2/"
source(paste(RscriptsPath,"basic_scripts.R",sep=""))

alldata=readRDS("~/Documents/Projects/HNSCC_PreCog/data/Precog.HNSCC.all.exp.clin.data.rds")
explist=alldata$explist
zscores.list=explist
infolist=alldata$infolist

infolist3=infolist

################################################
#Do survival analysis
##############################################

survlist=list()
for(i in 1:length(infolist3)){
  info=infolist3[[i]]
  timevars=colnames(info)[grep("_Time", colnames(info))]
  status=colnames(info)[grep("_Status", colnames(info))]
  survars=list(status, timevars)
  names(survars)=c("status", "time")
  print(names(infolist3)[i])
  print("OS_Time" %in%  timevars)
  survlist[[i]]=survars
}
names(survlist)=names(infolist3)

#make infolist.surv
survlist=survlist[as.numeric(which(unlist(lapply(survlist, function(x) length(x$status)!=0 & length(x$time)!=0))))]

#Find out which datasets have overall survival. For those that don't, find out which alternative survival variables they have 
survvartype=unlist(lapply(survlist, function(x) ifelse("OS_Time" %in% x$time & "OS_Status" %in% x$status, "OS", gsub("_Status","", x$status))))
survvartype[which(names(survvartype)=="GSE39366")]="RFS"

accessions=names(survvartype)

infolist.surv=infolist3[accessions]

#survvartype=survvartype[accessions]
for(i in 1:length(accessions)){
acc=accessions[i]
print(acc)
survvar=paste(survvartype, "_Status", sep="")[i]
survtimevar=paste(survvartype, "_Time", sep="")[i]
infolist.surv[[which(names(infolist.surv) %in% acc)]][,survvar]=as.numeric(as.character(infolist.surv[[which(names(infolist.surv) %in% acc)]][,survvar]))
infolist.surv[[which(names(infolist.surv) %in% acc)]][,survtimevar]=as.numeric(as.character(infolist.surv[[which(names(infolist.surv) %in% acc)]][,survtimevar]))
}
#

# 
#find out if the number of 0 and 1 survival status for each datset. Use only datasets with 5 or more samples with 0 or 1 in analysis 
surv.enough=list()
for(i in 1:length(accessions)){
acc=accessions[i]
info=infolist.surv[[which(names(infolist.surv) %in% acc)]]
survvar=paste(survvartype, "_Status", sep="")[i]
info[,survvar]=factor(info[,survvar], levels=c("0","1"))
surv.enough[[i]]=length(which(info[,survvar]==levels(info[,survvar])[1]))>=5 & length(which(info[,survvar]==levels(info[,survvar])[2]))>=5
}

infolist.surv=infolist.surv[which(unlist(surv.enough))]
accessions.survival=names(infolist.surv)
explist.survival=explist[accessions.survival]

survvartype=survvartype[accessions.survival]

##check class of survival status and time, both should be numeric
for(i in 1:length(accessions.survival)){
acc=accessions.survival[i]
print(acc)
info=infolist.surv[[which(names(infolist.surv) %in% acc)]]
survvar=paste(survvartype, "_Status", sep="")[i]
print(class(info[,survvar]))
survvar=paste(survvartype, "_Time", sep="")[i]
print(class(info[,survvar]))
}
#All good. All are numeric 

##################################################
#Make boxplot showig survival length for each dataset 
#######################################################

file=paste(Figuresdir, "Precog.meta.analysis.boxplot.survival.data")
pdf(file=paste(file,'.pdf',sep=''), height = 9,  width = 9, family = "Helvetica")
par(mfrow=c(4,5))
for(i in 1:length(accessions.survival)){
acc=accessions.survival[i]
info=infolist.surv[[which(names(infolist.surv) %in% acc)]]
survvar=paste(survvartype, "_Status", sep="")[i]
boxplot(info[,paste(survvartype, "_Time", sep="")[i]]~info[,paste(survvartype, "_Status", sep="")[i]], main=acc, ylab="survival (months)", xlab=survvar, outline=F, ylim=c(0,200))
stripchart(info[,paste(survvartype, "_Time", sep="")[i]]~info[,paste(survvartype, "_Status", sep="")[i]], add=T, ylim=c(0,200), vertical=T, method="jitter")
}
dev.off()

###############################################
#Make survival objects for all datasets 
###############################################

library(survival)

accessions.survival=names(infolist.surv)

Survival_objects=list()
for(i in 1:length(accessions.survival)){
acc=accessions.survival[i]
info=infolist.surv[[which(names(infolist.surv) %in% acc)]]
so=Surv(info[,paste(survvartype, "_Time", sep="")[i]], info[,paste(survvartype, "_Status", sep="")[i]])
names(so)=rownames(info)
so=so[!is.na(so)]
Survival_objects[[acc]]=so
}

lapply(Survival_objects, length)
#save survival objects for each dataset
#saveRDS(Survival_objects, paste(Resultsdir, "Precog.HNSCC.survival.objects.all.HNSCC.rds", sep=""))
#Add survival variable type for easy access 
Survival_times=list()
for(i in 1:length(accessions.survival)){
acc=accessions.survival[i]
info=infolist.surv[[which(names(infolist.surv) %in% acc)]]
info=info[names(Survival_objects[[acc]]),]
sum=summary(info[,paste(survvartype, "_Time", sep="")[i]])
Survival_times[[acc]]=sum
}

###
surv.numbers=list()
for(i in 1:length(accessions.survival)){
acc=accessions.survival[i]
info=infolist.surv[[which(names(infolist.surv) %in% acc)]]
info=info[names(Survival_objects[[acc]]),]
survvar=paste(survvartype, "_Status", sep="")[i]
info[,survvar]=factor(as.character(info[,survvar]), levels=c("0","1"))
surv.numbers[[i]]=c(length(which(info[,survvar]==levels(info[,survvar])[1])), length(which(info[,survvar]==levels(info[,survvar])[2])))
names(surv.numbers[[i]])=c("0","1")
}
names(surv.numbers)=accessions.survival

#record the number of patients in each survival category
saveRDS(surv.numbers, paste(Resultsdir, "Precog.HNSCC.survival.patient.numbers.rds", sep=""))
all(names(surv.numbers) %in% names(Survival_objects))

#names(survvartype)==(names(Survival_objects))
#names(surv.numbers)==(names(Survival_objects))
Survival_Allinfo=list(Survival_objects, survvartype, surv.numbers, Survival_times)
names(Survival_Allinfo)=c("Survival_objects", "survvartype", "surv.numbers","surv.times")
saveRDS(Survival_Allinfo, paste(Resultsdir, "Precog.HNSCC.survival.objects.all.HNSCC.info.rds", sep=""))
Survival_Allinfo=readRDS(paste(Resultsdir, "Precog.HNSCC.survival.objects.all.HNSCC.info.rds", sep=""))

x="TCGA"
!is.na(Survival_Allinfo$Survival_objects[[x]])

#Double-check that expression dataframe colnames match survival object names 
lapply(1:length(explist.survival), function(x) all(colnames(explist.survival[[x]])==names(Survival_objects[[x]])))

########################################################
#Do survival meta-analysis in all HNSCC
########################################################

#
alldata=readRDS("~/Documents/Projects/HNSCC_PreCog/data/Precog.HNSCC.all.exp.clin.data.rds")
explist=alldata$explist
zscores.list=explist
infolist=alldata$infolist

surv.numbers=readRDS(paste(Resultsdir, "Precog.HNSCC.survival.patient.numbers.rds", sep=""))
accessions.survival=names(surv.numbers)

infolist.surv=infolist[which(names(infolist) %in% accessions.survival)]
accessions.survival=names(infolist.surv)

#
RscriptsPath="~/Documents/scripts/rscripts2/"
source(paste(RscriptsPath,"coxph.apply.R",sep=""))

###################################

#Get survival objects. see if I can use these to do the meta analysis 
Survival_objects=readRDS(paste(Datadir, "Precog.HNSCC.survival.objects.all.HNSCC.rds", sep=""))

#Restrict to studies with at least 20 patients with survival data 
keep=names(which(unlist(lapply(Survival_objects, function(x) length(x)>=20))))
#17 studies
Survival_objects=Survival_objects[keep]
accessions.survival=names(Survival_objects)

accessions.survival=intersect(accessions.survival, names(explist))
accessions.survival=setdiff(accessions.survival, "GSE686")
#Down to 16 studies with both gene expression results and survival data 

#need to further restrict to studies that have 20 samples after exclusion of samples that do not have expression data
coxphlengths=list()
for(i in 1:length(accessions.survival)){
acc=accessions.survival[i]
Survival_object=Survival_objects[[which(names(Survival_objects)==acc)]]
exp=explist[[which(names(explist)==acc)]]
#exp=as.data.frame(exp)
OverlapSamples=intersect(names(Survival_object), colnames(exp))
Survival_object=Survival_object[OverlapSamples]
coxphlengths[[acc]]=Survival_object
}

#
accessions.survival=names(which(unlist(lapply(coxphlengths, function(x) length(x)>=20))))
#OK 16 studies with over 20 or more samples with both gene expression data and survival

lapply(Survival_objects, length) %>% unlist() %>% sum

#
source(paste(RscriptsPath, "coxph.apply.R", sep=""))

allcoxph=list()
for(i in 1:length(accessions.survival)){

acc=accessions.survival[i]

Survival_object=Survival_objects[[which(names(Survival_objects)==acc)]]

exp=explist[[which(names(explist)==acc)]]
exp=as.data.frame(exp)

OverlapSamples=intersect(names(Survival_object), colnames(exp))
Survival_object=Survival_object[OverlapSamples]
exp=exp[,OverlapSamples]

survs=apply(exp,1, coxph.apply.provide.survobject.with.pvalue)
survs=survs[which(lapply(survs, length)==5)]
survs=as.data.frame(abind(survs, along=1))
survs=survs[!is.na(survs$Z),]
survs$Z=as.numeric(as.character(survs$Z))
survs$P=as.numeric(as.character(survs$P))
survs$Q=p.adjust(survs$P, method="fdr")
survs$gene=rownames(survs)

allcoxph[[i]]=survs
}
names(allcoxph)=accessions.survival

saveRDS(allcoxph, paste(Resultsdir, "Precog.HNSCC.coxph.survival.allstudies.allHNSCC.updated.06_09_2023.rds", sep=""))
allcoxph=readRDS(paste(Resultsdir, "Precog.HNSCC.coxph.survival.allstudies.allHNSCC.updated.06_09_2023.rds", sep=""))

##############################################################
#Perform meta-analysis of prognostic scores using the Liptak method 
##############################################################

library(abind)
allcoxph.abind=as.data.frame(abind(allcoxph, along=1))
rownames(allcoxph.abind)=NULL
allcoxph.abind$Z=as.numeric(as.character(allcoxph.abind$Z))
allcoxph.abind$N=as.numeric(as.character(allcoxph.abind$N))

#restricting to genes that found within at least two studies, otherwise its not a meta-analysis
tab=table(allcoxph.abind$gene)
meta.genes=names(tab[which(tab>=2)])
length(meta.genes) #23558 genes for which at least two studies have the gene. Running meta analysis on all of these 

allcoxph.abind$gene=as.character(allcoxph.abind$gene)

metaz.coxph=lapply(meta.genes, function(x)  Liptak_combine_z(gene=x, df=allcoxph.abind))
names(metaz.coxph)=meta.genes
metaz.coxph=unlist(metaz.coxph)

#Annotate genes using biomaRt
library(biomaRt)
ensembl=useMart("ensembl", dataset="hsapiens_gene_ensembl")

annot=getBM(attributes=c("hgnc_symbol","entrezgene_id"), filters = "entrezgene_id", values = names(metaz.coxph), mart= ensembl)

metaz.coxph.annot=cbind(metaz.coxph, annot[match(meta.genes, annot$entrezgene_id),])
metaz.coxph.annot$gene=rownames(metaz.coxph.annot)

#saveRDS(metaz.coxph.annot, paste(Resultsdir, "Precog.HNSCC.coxph.Liptak.metaz.allHNSCC.updated.03_30_2020.rds", sep=""))
#metaz.coxph.annot=readRDS(paste(Resultsdir, "Precog.HNSCC.coxph.Liptak.metaz.allHNSCC.updated.03_30_2020.rds", sep=""))
#add the number of studies for which each gene was available
metaz.coxph.annot$n.studies=tab[match(metaz.coxph.annot$gene, names(tab))]

threshold="3.09"
metaz.coxph.annot$signficant=ifelse(abs(metaz.coxph.annot$metaz.coxph)>=threshold,"yes","no")

metaz.coxph.annot.pos=metaz.coxph.annot[metaz.coxph.annot$signficant=="yes" & metaz.coxph.annot$metaz.coxph<0,]
metaz.coxph.annot.pos=metaz.coxph.annot.pos[order(metaz.coxph.annot.pos$metaz.coxph),]
metaz.coxph.annot.pos 
#479 genes positively prognostic genes 

metaz.coxph.annot.neg=metaz.coxph.annot[metaz.coxph.annot$signficant=="yes" & metaz.coxph.annot$metaz.coxph>0,]
metaz.coxph.annot.neg=metaz.coxph.annot.neg[order(metaz.coxph.annot.neg$metaz.coxph, decreasing=T),]
#730 negatively prognostic genes

metaz.coxph.annot$survival.direction=rep(NA, nrow(metaz.coxph.annot))
metaz.coxph.annot$survival.direction[metaz.coxph.annot$metaz.coxph<=(-3.09)]="Pro-survival"
metaz.coxph.annot$survival.direction[metaz.coxph.annot$metaz.coxph>=(3.09)] ="Anti-survival"
metaz.coxph.annot$survival.direction=as.factor(metaz.coxph.annot$survival.direction)
summary(metaz.coxph.annot$survival.direction)

saveRDS(metaz.coxph.annot, paste(Resultsdir, "Precog.HNSCC.coxph.Liptak.metaz.allHNSCC.updated.06_09_2023.rds", sep=""))

#Try remaking phenograph clusters
metaz.coxph.annot=readRDS(paste(Resultsdir, "Precog.HNSCC.coxph.Liptak.metaz.allHNSCC.updated.06_09_2023.rds", sep=""))
pvals.sig=metaz.coxph.annot[!is.na(metaz.coxph.annot$survival.direction),"gene"]

#Make a combined matrix of gene expression data from studies, which will be used to cluster genes based on coexpression 
alldata=readRDS("~/Documents/Projects/HNSCC_PreCog/data/Precog.HNSCC.all.exp.clin.data.rds")
explist=alldata$explist
zscores.list=explist
infolist=alldata$infolist

allgenes=lapply(explist, rownames)
#for each dataset get the intersection with the signficant gene list 
allgenes2=lapply(allgenes, function(x) intersect(pvals.sig, x))

#somwhow make a dataframe showing presence or absence in each dataset 
genes=unique(unlist(allgenes2))
length(genes)
#1209 genes altogether 

genearray=array(NA, c(length(allgenes2), length(genes)))
rownames(genearray)=names(allgenes2)
colnames(genearray)=genes
genearray=as.data.frame(genearray)

gene.present=lapply(allgenes2, function(x)  ifelse(genes %in% x, "yes","no"))
library(abind)
gene.present.df=as.data.frame(abind(gene.present, along=2))
rownames(gene.present.df)=genes

library(plyr)
mat=apply(gene.present.df, 2, function(x) as.numeric(as.character(revalue(as.factor(x), c("yes"="1","no"="0")))))
rownames(mat)=rownames(gene.present.df)

#studies with over 80% of genes 
studiesallgenes=names(which((apply(gene.present.df, 2, function(x) length(which(x=="yes")))/length(genes))>=0.8))
length(studiesallgenes)
mat=mat[,studiesallgenes]

#Now get common genes
gene.present.df.allgenes=gene.present.df[,studiesallgenes]
gene.present.df.allgenes=gene.present.df.allgenes[which(apply(gene.present.df.allgenes, 1, function(x) all(x=="yes"))),]
dim(gene.present.df.allgenes) #958 genes in 20 studies 
958/1209

explist.sig=explist[which(names(explist) %in% colnames(gene.present.df.allgenes))]
ex.commongenes=lapply(explist.sig, function(x) x[match(genes, rownames(x)),])
ex.commongenes.df=abind(ex.commongenes, along=2)
ex.commongenes.df=na.omit(ex.commongenes.df)
ex.commongenes.df=ex.commongenes.df
dim(ex.commongenes.df) #958 1642

library(iCellR)
set.seed(123)
Rphenograph_out <- iCellR::Rphenograph(ex.commongenes.df, k = 45)
#6 clusters, hopefully the same

saveRDS(Rphenograph_out, paste0(Resultsdir,"phenograph.survival.genes.060923.rds"))
Rphenograph_out=readRDS(paste0(Resultsdir,"phenograph.survival.genes.060923.rds"))

df=data.frame(gene=rownames(ex.commongenes.df), pheno.clust=Rphenograph_out[[2]]$membership) %>% 
  mutate(pheno.clust=glue::glue("S{pheno.clust}"))

metaz.coxph.annot$survival.gene.cluster.phenograph=df[match(metaz.coxph.annot$gene, df$gene),"pheno.clust"]

saveRDS(metaz.coxph.annot, paste(Resultsdir, "Precog.HNSCC.coxph.Liptak.metaz.allHNSCC.updated.06_09_2023.rds", sep=""))
metaz.coxph.annot=readRDS(paste0(Resultsdir, "Precog.HNSCC.coxph.Liptak.metaz.allHNSCC.updated.06_09_2023.rds"))


###############################################################
#Make UMAP of survival genes 
###############################################################

Rphenograph_out=readRDS(paste0(Resultsdir,"phenograph.survival.genes.060923.rds"))

library(umap)
set.seed(234)
clust.umap = umap(ex.commongenes.df)
saveRDS(clust.umap, paste0(Resultsdir, "UMAP_PRECOG_HNSCC_survival_gene_clusters.updated.061623.rds"))

surv.cols=gnuplot_colors[c(1:5,7)]

####
fm=readRDS(paste0(Datadir, "Findmarkers_Puram_U54_allcombined.update.061323.rds"))
fm2=fm[match(rownames(clust.umap$layout), fm$gene),]

file=paste(Figuresdir, "UMAP_PRECOG_HNSCC_survival_genes_updated_061723", sep="")
pdf(file=paste0(file,'.pdf',sep=''), height = 5,  width = 5, family = "Helvetica")
plot(clust.umap$layout, ylab="UMAP_2", xlab="UMAP_1", pch=as.numeric(fm2$survival.direction), cex=0.6, col=surv.cols[as.numeric(fm2$survival.gene.cluster.phenograph)])
legend("topleft", legend=c("Pro-survival","Anti-survival"), pch=c(1,2), title="Survival direction", cex=1, bty="n")
legend("bottomright", legend=c("S1","S2","S3","S4","S5","S6"), col=surv.cols, pch=15, title="Survival gene\ncluster", cex=1, bty="n")
dev.off()

###############################
#Boxplot of survival gene clusters 
##################################

metaz.coxph.annot=readRDS(paste0(Resultsdir, "Precog.HNSCC.coxph.Liptak.metaz.allHNSCC.updated.06_09_2023.rds"))

#survival
dat=data.frame(Meta.z=fm$metaz.coxph, Gene.cluster=fm$survival.gene.cluster.phenograph, dir=fm$survival.direction, gene=fm$symbol_AnnotationDbi)
dat=dat[!is.na(fm$survival.direction) & !is.na(dat$Meta.z),]

dat$Gene.cluster2=factor(ifelse(is.na(dat$Gene.cluster),"S7",as.character(dat$Gene.cluster)))
dat$dir=factor(dat$dir, levels=c("Anti-survival", "Pro-survival"))

n=10
dat$Meta.z.rank=rep(1, nrow(dat))
dat$Meta.z.rank[order(dat$Meta.z,decreasing=T)[1:n]]=3
dat$Meta.z.rank[order(dat$Meta.z,decreasing=F)[1:n]]=2
dat$Meta.z.rank=factor(dat$Meta.z.rank)

t1 <- ggplot(transform(dat, xjit=jitter(as.numeric(as.factor(Gene.cluster2)))), aes(x=as.factor(Gene.cluster2),y=Meta.z))
t2 <- geom_boxplot(outlier.shape=NA, aes(fill=Gene.cluster2))
t3 <- geom_point(size=0.8, aes(x=xjit, color = Meta.z.rank))
t5 <- geom_text_repel(aes(x=xjit, label=ifelse(Meta.z.rank!=1,as.character(gene),'')), vjust = 0.5, size=2.5)

p=t1+t2+t3 + facet_grid(dir ~ ., scale='free_y') +t5
p=p + scale_color_manual(values = c("black", "red", "blue"), labels=c("Not top rank","Highest 10","Lowest 10"))+labs(color="Meta-Z rank", x="Gene cluster", y="Meta-z score")
p=p+ scale_x_discrete(labels=c("S1","S2","S3","S4","S5","S6","Not clustered"))

surv.cols=gnuplot_colors[c(1:5,7)]
cols2=c(surv.cols, "lightgrey")

p=p+ scale_fill_manual(values=cols2, name="Survival gene cluster", labels=c("Not.clustered"="No cluster"))

file=paste0(Figuresdir, "jitter_clusters_survival_updated_061723")
pdf(file=paste0(file,'.pdf',sep=''), height = 6,  width = 8, family = "Helvetica")
p
dev.off()

#####################################################################
#Make combined table of survival and lnm-assocaiated genes
#####################################################################

library(HGNChelper)
library("AnnotationDbi")
library("org.Hs.eg.db")

metatab.surv=metaz.coxph.annot %>% filter(!is.na(survival.direction)) 

metatab.lnm=readRDS(paste0(Datadir, "zscores_lnm_all_genes_updated.06.09.23.rds")) %>% 
  dplyr::rename('lnm.gene.cluster.phenograph' = 'cluster', 'metaz.lnm'='zval')

#Make combined table of survival and LNM-associated genes 

surv.genes=metatab.surv %>% dplyr::select("gene","hgnc_symbol","metaz.coxph","survival.direction","survival.gene.cluster.phenograph")

node.genes=metatab.lnm %>% dplyr::select("gene","hgnc_symbol","metaz.lnm","lnm_direction","lnm.gene.cluster.phenograph") %>% 
  dplyr::rename("lnm.direction"="lnm_direction")
  
allgenes=unique(c(surv.genes$gene, node.genes$gene))
length(allgenes) #2001 genes
sv2=surv.genes[match(allgenes, surv.genes$gene),] %>% dplyr::select(-c(gene, hgnc_symbol))
lnm2=node.genes[match(allgenes, node.genes$gene),] %>% dplyr::select(-c(gene, hgnc_symbol))

allgenes2=cbind(allgenes, sv2, lnm2) %>% 
  dplyr::rename("gene"="allgenes")

#Add gene annotation
library(biomaRt)
ensembl=useMart("ensembl", dataset="hsapiens_gene_ensembl")

bm.gene=getBM(attributes=c("hgnc_symbol","entrezgene_id"), filters = "entrezgene_id", values = as.character(allgenes2$gene), mart= ensembl)
bm.gene=bm.gene[match(allgenes2$gene, bm.gene$entrezgene_id),]
allgenes2[,c("hgnc_symbol","entrezgene_id")]=bm.gene[,c("hgnc_symbol","entrezgene_id")]
dim(allgenes2) #2001 gene
length(which(is.na(allgenes2$hgnc_symbol))) #24 missing

##
library("AnnotationDbi")
library("org.Hs.eg.db")

allgenes2$symbol_AnnotationDbi=as.character(mapIds(org.Hs.eg.db, keys=allgenes2$gene, column="SYMBOL", keytype="ENTREZID", multiVals="first"))

length(which(is.na(allgenes2$symbol_AnnotationDbi))) #6 missing

#Match gene symbols to Puram and Stanford scRNA-Seq datasets
int=readRDS("~/Documents/Projects/U54/Puram/data/Integated_Puram_patients200cells_allwithlnm.rds")
ga=as.data.frame(as.matrix(GetAssayData(object = int$RNA)))

genes3=allgenes2$symbol_AnnotationDbi
missinggenes=setdiff(genes3, rownames(int$RNA))
length(missinggenes) #117 genes

library(HGNChelper)
#use hynchelper to update names of missing genes
hg=HGNChelper::checkGeneSymbols(missinggenes, unmapped.as.na=FALSE)
intersect(hg$Suggested.Symbol, rownames(int$RNA))#change "MPP6" "PHB" 

#check gene symbolds of Puram data and search again
hg.int=checkGeneSymbols(rownames(int$RNA), unmapped.as.na=FALSE)
intersect(hg$Suggested.Symbol, hg.int$Suggested.Symbol)
intersect(missinggenes, hg.int$Suggested.Symbol)

change=intersect(missinggenes, hg.int$Suggested.Symbol)
change=hg.int[hg.int$Suggested.Symbol %in% change,]
allgenes2$Matched_to_Puram=as.character(plyr::mapvalues(as.factor(allgenes2$symbol_AnnotationDbi), from=c(change$Suggested.Symbol), to=c(change$x)))

#convert NEMP1 to TMEM194A in gene list so it matched Puram
hg.int[hg.int$Suggested.Symbol=="PHB",]
hg[hg$Suggested.Symbol=="PHB",]
allgenes2$Matched_to_Puram=gsub("PHB1", "PHB", allgenes2$Matched_to_Puram)

hg[hg$Suggested.Symbol=="MPP6",]
hg.int[hg.int$Suggested.Symbol=="MPP6",]
allgenes2$Matched_to_Puram=gsub("PALS2", "MPP6", allgenes2$Matched_to_Puram)

length(which(is.na(allgenes2$Matched_to_Puram))) #Only 6 genes remain missing 

#Match genes to Stanford 
#get Stanford scRNA-Seq dataset and add matching column 
int=readRDS("~/Documents/Projects/U54/CCSB_scRNASeq/data/CCSB_scRNASeq_HNSCC_all_enzymatic_intgrated_50PCs.mito.removed.rds")
DefaultAssay(int)="RNA"

missinggenes=setdiff(genes3, rownames(int$RNA))
length(missinggenes) #143 missing genes

hg=HGNChelper::checkGeneSymbols(missinggenes, unmapped.as.na=FALSE)
intersect(hg$Suggested.Symbol, rownames(int$RNA))#change "C5orf66-AS1" "MPP6"        "PHB" 

#check gene symbolds of Puram data and search again
hg.int=checkGeneSymbols(rownames(int$RNA), unmapped.as.na=FALSE)
intersect(hg$Suggested.Symbol, hg.int$Suggested.Symbol)
intersect(missinggenes, hg.int$Suggested.Symbol)

change=intersect(missinggenes, hg.int$Suggested.Symbol)
change=hg.int[hg.int$Suggested.Symbol %in% change,]

allgenes2$Matched_to_Stanford=as.character(plyr::mapvalues(as.factor(allgenes2$symbol_AnnotationDbi), from=c(change$Suggested.Symbol), to=c(change$x)))

hg.int[hg.int$Suggested.Symbol=="PHB",]
hg[hg$Suggested.Symbol=="PHB",]
allgenes2$Matched_to_Stanford=gsub("PHB1", "PHB", allgenes2$Matched_to_Stanford)

hg.int[hg.int$Suggested.Symbol=="MPP6",]
hg[hg$Suggested.Symbol=="MPP6",]
allgenes2$Matched_to_Stanford=gsub("PALS2", "MPP6", allgenes2$Matched_to_Stanford)

saveRDS(allgenes2, paste0(Datadir, "all.genes.survival.lnm.060923.rds"))

#########################################################
#Adding mean normalized counts for each cell type in the Puram and Stanford scRNA-Seq datasets
#Identifying the cell type with max expression in each 
#########################################################

#Adding mean and max normalized counts for all prognostic genes in Puram scRNA-Seq dataset

library(Seurat)
celltypes=c("B.cell", "Endothelial", "Fibroblast", "Malignant", "Mast", "Myeloid", "T.cell.or.NK.cell")

int=readRDS(paste0(DatadirPuram, "Integated_Puram_prim_patients200cells.rds"))
DefaultAssay(int)="RNA"

keep=colnames(int[,int$cell.type.collapsed!="Unclassified" & int$cell.type.collapsed!="myocyte"])
int=subset(int, cells=keep)

cr=int@assays$RNA@counts #raw counts
cn=int@assays$RNA@data #normalizard counts

cm.cr=colMeans(cr)
cm.cn=colMeans(cn)

#Get marix mean levels of each gene
#fm=readRDS(paste0(Datadir, "all.genes.survival.lnm.060923.rds"))
fm=readRDS(paste0(Datadir, "Findmarkers_Puram_U54_allcombined.update.061323.rds"))

ga=cn
genes=fm$Matched_to_Puram[fm$Matched_to_Puram %in% rownames(ga)]
ga=ga[genes,]

meantab=apply(ga, 1, function(x) aggregate(x~int$cell.type.collapsed, FUN=mean)[,2])
meantab=as.data.frame(t(meantab))
colnames(meantab)=celltypes

meantab=meantab[match(fm$Matched_to_Puram, rownames(meantab)),]
meantab.celltype.Puram=meantab

fm[,paste0("Puram_mean_normalized_counts_", colnames(meantab))]=meantab
cols=paste0("Puram_mean_normalized_counts_", colnames(meantab))
fm$Max_celltype_normalized_counts_Puram=apply(fm[,cols],1, function(x) celltypes[which.max(x)])
fm$Max_celltype_normalized_counts_Puram=as.character(fm$Max_celltype_normalized_counts_Puram)
fm$Max_celltype_normalized_counts_Puram=plyr::revalue(as.factor(fm$Max_celltype_normalized_counts_Puram), c("character(0)"=NA))

#Adding mean normalized counts for all prognostic genes in Stanford scRNA-Seq dataset

int=readRDS(paste0(DatadirCCSBSc, "CCSB_scRNASeq_HNSCC_primary_enzymatic_intgrated_50PCs.mito.removed.rds"))
DefaultAssay(int)="RNA"

keep=colnames(int[,int$celltype!="Unclassified"])
int=subset(int, cells=keep)

#Matching cell type factor order to Puram dataset
int$cell.type.collapsed=factor(int$cell.type.collapsed, levels=c("B cell","Endothelial","Fibroblast","Malignant","Mast", "Myeloid","T cell or NK cell"))

cn=int@assays$RNA@data #normalizard counts
cr=int@assays$RNA@counts #raw counts

cm.cr=colMeans(cr)
cm.cn=colMeans(cn)

ga=cn
genes=fm$Matched_to_Stanford[fm$Matched_to_Stanford %in% rownames(ga)]
ga=ga[genes,]

meantab=apply(ga, 1, function(x) aggregate(x~int$cell.type.collapsed, FUN=mean)[,2])
meantab=as.data.frame(t(meantab))
colnames(meantab)=celltypes
meantab=meantab[match(fm$Matched_to_Stanford, rownames(meantab)),]

meantab.celltype.Stanford=meantab
#save mean normalized counts per cell matrices for Puram and Stanford datasets
meantabs=list(celltype.Puram=meantab.celltype.Puram, celltype.Stanford=meantab.celltype.Stanford)
saveRDS(meantabs, paste0(Datadir, "matrices.mean.normalized.counts.update.061323.rds"))

fm[,paste0("Stanford_mean_normalized_counts_", colnames(meantab))]=meantab
cols=paste0("Stanford_mean_normalized_counts_", colnames(meantab))
fm$Max_celltype_normalized_counts_Stanford=apply(fm[,cols],1, function(x) celltypes[which.max(x)])
fm$Max_celltype_normalized_counts_Stanford=as.character(fm$Max_celltype_normalized_counts_Stanford)
fm$Max_celltype_normalized_counts_Stanford=plyr::revalue(as.factor(fm$Max_celltype_normalized_counts_Stanford), c("character(0)"=NA))

library(tidyverse)
fm=fm %>% 
  mutate(lnm.gene.cluster.phenograph=factor(lnm.gene.cluster.phenograph), survival.gene.cluster.phenograph=factor(survival.gene.cluster.phenograph), lnm.direction=factor(lnm.direction))

#saveRDS(fm, paste0(Datadir, "Findmarkers_Puram_U54_allcombined.update.061323.rds"))
#Adding max normalized counts in Stanford bulk RNA-Seq dataset

cb=readRDS(paste0(DatadirU54Bulk, "U54.normalized.counts.allsamples.gene.names.rds"))
sampleTable=readRDS(paste0(DatadirU54Bulk, "U54_data_clinical_formatted_all.rds"))

sampleTable=sampleTable[sampleTable$sample_type %in% c("E","F","L","T") & sampleTable$cancertype=="HNSCC" & sampleTable$ln=="0",]
all(rownames(sampleTable) %in% colnames(cb)) #TRUE

cb=cb[,match(rownames(sampleTable), colnames(cb))]
cb=cb[rownames(cb) %in% fm$symbol_AnnotationDbi,]
#1967 genes

meantab=apply(cb, 1, function(x) aggregate(x~sampleTable$sample_type, FUN=mean)[,2])

levels(as.factor(sampleTable$sample_type))

celltypes=c("Endothelial","Fibroblast","Leukocyte","Malignant")

meantab=as.data.frame(t(meantab))
colnames(meantab)=celltypes

meantab=meantab[match(fm$symbol_AnnotationDbi, rownames(meantab)),]
meantab.celltype.Bulk.Stanford=meantab
#2001 genes

meantabs$celltype.Stanford.Bulk=meantab.celltype.Bulk.Stanford
saveRDS(meantabs, paste0(Datadir, "matrices.mean.normalized.counts.update.061323.rds"))

fm[,paste0("Stanford_mean_normalized_counts_bulk_", colnames(meantab))]=meantab

celltypes=c("Endothelial","Fibroblast","Leukocyte","Malignant")

cols=paste0("Stanford_mean_normalized_counts_bulk_", colnames(meantab))
fm$Max_celltype_normalized_counts_bulk=apply(fm[,cols],1, function(x) celltypes[which.max(x)])
fm$Max_celltype_normalized_counts_bulk=as.character(fm$Max_celltype_normalized_counts_bulk)
fm$Max_celltype_normalized_counts_bulk=plyr::revalue(as.factor(fm$Max_celltype_normalized_counts_bulk), c("character(0)"=NA))

saveRDS(fm, paste0(Datadir, "Findmarkers_Puram_U54_allcombined.update.061323.rds"))

#####################################################################
#Barplot showing overlap of survival and LNM gene clusters 
#####################################################################

tab=fm %>% dplyr::select(survival.gene.cluster.phenograph, lnm.gene.cluster.phenograph) %>% na.omit()
#surv.cols=gnuplot_colors[c(1:5)]
lnm.cols=gnuplot_colors[c(8:11,13:14)]
surv.cols=gnuplot_colors[c(1:5,7)]

p=ggplot(tab, aes(lnm.gene.cluster.phenograph, ..count..)) + geom_bar(aes(fill = survival.gene.cluster.phenograph), position = "dodge") + 
  scale_fill_manual(values=surv.cols, name="Survival gene\ncluster") +
  labs(x="LNM gene cluster", y="# overlapping genes") + theme_bw()

file=paste0(Figuresdir, "Barplot_overlap_survival_clusters_LNM_clusters_update2023")
pdf(file=paste0(file,'.pdf',sep=''), height = 5,  width = 6, family = "Helvetica")
p
dev.off()

#Hypergeometric test for overlap of clusters L2 and S6 
GenesS6=fm[which(fm$survival.gene.cluster.phenograph=="S6"),"gene"]
GenesL2=fm[which(fm$lnm.gene.cluster.phenograph=="L2"),"gene"]

arraygenes=readRDS(glue::glue(Resultsdir, "Precog.meta.analysis.genes.all.HNSCC.updated.06.09.23.rds"))

phyper(length(intersect(GenesS6,GenesL2))-1,length(intersect(GenesS6, arraygenes)),length(setdiff(arraygenes, GenesS6)),length(GenesL2),lower.tail=F, log.p = FALSE)
#30 overlapping genes 
#p=5.495718e-31

genes=fm[which(fm$survival.gene.cluster.phenograph=="S6" & fm$lnm.gene.cluster.phenograph=="L2"),"symbol_AnnotationDbi"]


###########################################################################
#Expression of survival and LNM gene clusters in bulk U54 data 
###########################################################################

cb=readRDS(paste0(DatadirU54Bulk, "U54.normalized.counts.allsamples.gene.names.rds"))
sampleTable=readRDS(paste0(DatadirU54Bulk, "U54_data_clinical_formatted_all.rds"))

sampleTable=sampleTable[sampleTable$sample_type %in% c("E","F","L","T") & sampleTable$cancertype=="HNSCC" & sampleTable$ln=="0",]
all(rownames(sampleTable) %in% colnames(cb)) #TRUE

cb=cb[,match(rownames(sampleTable), colnames(cb))]
cb=cb[rownames(cb) %in% fm$symbol_AnnotationDbi,]
#1974 genes

fm=readRDS(paste0(Datadir, "Findmarkers_Puram_U54_allcombined.update.061323.rds"))

makesig=function(genetype, group){
  levs=levels(fm[,group])
  sigs=lapply(levs, function(lev) 
  {x=fm[fm[,group]==lev,genetype]
  x=x[!is.na(x)]
  return(x)})
  names(sigs)=levs
  return(sigs)
}

sigs=c(
  makesig(genetype="symbol_AnnotationDbi",group="survival.direction"),
  makesig(genetype="symbol_AnnotationDbi",group="survival.gene.cluster.phenograph"),
  makesig(genetype="symbol_AnnotationDbi",group="lnm.direction"),
  makesig(genetype="symbol_AnnotationDbi",group="lnm.gene.cluster.phenograph"))

sigs=lapply(sigs, function(x) x[x %in% rownames(cb)])

sampleTable$sample_type=plyr::revalue(as.factor(sampleTable$sample_type), c("E"="Endothelial", "F"="Fibroblast", "L"="Leukocyte", "T"="Malignant"))

df=as.data.frame(abind::abind(lapply(1:length(sigs), function(x) as.data.frame(cbind(exp=colMeans(cb[sigs[[x]],]), Cell.type=as.character(sampleTable$sample_type), sig=rep(names(sigs)[x], ncol(cb))))), along=1))
df$exp=as.numeric(as.character(df$exp))
df$Cell.type=as.factor(df$Cell.type)

df$sig=factor(df$sig, levels=c("Anti-survival","Pro-survival","S1","S2","S3","S4","S5","S6","Anti-LNM","Pro-LNM","L1","L2","L3","L4","L5","L6"))
df$Outcome=factor(ifelse(df$sig %in% levels(df$sig)[1:8],"Survival","LNM"), levels=c("Survival","LNM"))

df$exp2=log2(df$exp)

df1=df[df$Outcome=="Survival",]
df2=df[df$Outcome=="LNM",]

p1.2=ggplot(df1, aes(x=Cell.type, y=exp2, fill=Cell.type))+ geom_jitter(position=position_jitterdodge(jitter.width =0.2), size=0.3, alpha=0.8, pch=21) +geom_violin(aes(fill = factor(Cell.type)), alpha=0.6)+ facet_grid(~sig,scales="free", space="free_x") + scale_fill_manual(name = "Cell type", values = cbp[c(2,3,8,4)])+ ylab("Mean log2 expression")+ xlab("Cell type")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) 
p2.2=ggplot(df2, aes(x=Cell.type, y=exp2, fill=Cell.type))+ geom_jitter(position=position_jitterdodge(jitter.width =0.2), size=0.3, alpha=0.8, pch=21) +geom_violin(aes(fill = factor(Cell.type)), alpha=0.6)+ facet_grid(~sig,scales="free", space="free_x") + scale_fill_manual(name = "Cell type", values = cbp[c(2,3,8,4)])+ ylab("Mean log2 expression")+ xlab("Cell type")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) 

file=paste0(Figuresdir, "violinplots.celltype.celltype_updated_061623")
pdf(file=paste0(file,'.pdf',sep=''), height = 8,  width = 8, family = "Helvetica")
ggpubr::ggarrange(p1.2, p2.2, ncol=1, nrow=2, common.legend = TRUE, legend="right")
dev.off()

p1.2=ggplot(df1, aes(x=Cell.type, y=exp2, fill=Cell.type))+ geom_jitter(position=position_jitterdodge(jitter.width =0.2), size=0.3, alpha=0.8, pch=21) +geom_violin(aes(fill = factor(Cell.type)), alpha=0.6)+ facet_grid(~sig,scales="free", space="free_x") + scale_fill_manual(name = "Cell type", values = cbp[c(2,3,8,4)])+ ylab("Mean log2 expression")+ xlab("Cell type")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) +theme(legend.position = "none")
p2.2=ggplot(df2, aes(x=Cell.type, y=exp2, fill=Cell.type))+ geom_jitter(position=position_jitterdodge(jitter.width =0.2), size=0.3, alpha=0.8, pch=21) +geom_violin(aes(fill = factor(Cell.type)), alpha=0.6)+ facet_grid(~sig,scales="free", space="free_x") + scale_fill_manual(name = "Cell type", values = cbp[c(2,3,8,4)])+ ylab("Mean log2 expression")+ xlab("Cell type")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) +theme(legend.position = "none")

file=paste0(Figuresdir, "violinplots.celltype.celltype.no.legend_updated_061623")
pdf(file=paste0(file,'.pdf',sep=''), height = 8,  width = 7, family = "Helvetica")
ggpubr::ggarrange(p1.2, p2.2, ncol=1, nrow=2)
dev.off()

#############################
#Making violin plots in bulk RNA-Seq for PanglaoDB cell type markers 
#############################

pang=read.table(paste0(DatadirU54, "PanglaoDB_markers_27_Mar_2020.tsv"), sep="\t", header=T,  fill=T, quote="", comment.char="")
pang=pang[pang$species=="Hs"|pang$species=="Mm Hs",]
celltypes=levels(as.factor(pang[pang$organ=="Immune system","cell.type"]))

library(HGNChelper)
#use hynchelper to update names of missing genes
hg=HGNChelper::checkGeneSymbols(pang$official.gene.symbol, unmapped.as.na=FALSE)
pang$hgnc=hg$Suggested.Symbol

#match to gene symbols for bulk sorted cells
cb=readRDS(paste0(DatadirU54Bulk, "U54.normalized.counts.allsamples.gene.names.rds"))
sampleTable=readRDS(paste0(DatadirU54Bulk, "U54_data_clinical_formatted_all.rds"))

sampleTable=sampleTable[sampleTable$sample_type %in% c("E","F","L","T") & sampleTable$cancertype=="HNSCC" & sampleTable$ln=="0",]
all(rownames(sampleTable) %in% colnames(cb)) #TRUE

cb=cb[,match(rownames(sampleTable), colnames(cb))]
cb=cb[rownames(cb) %in% fm$symbol_AnnotationDbi,]

hg.cb=checkGeneSymbols(rownames(cb), unmapped.as.na=FALSE)
hg.cb=hg.cb[match(hg$Suggested.Symbol, hg.cb$Suggested.Symbol),]

hg$Bulk_scRNASeq_symbol=hg.cb$x
all(pang$hgnc==hg$Suggested.Symbol)
pang$Bulk_scRNASeq_symbol=hg$Bulk_scRNASeq_symbol

celltypes=c("Endothelial cells", "Epithelial cells", "Fibroblasts", "Mast cells", "B cells", "T cells", "Macrophages", "Dendritic cells")

pang=pang[pang$cell.type %in% celltypes,]

sigs=lapply(celltypes, function(celltype) pang[pang$cell.type==celltype,"Bulk_scRNASeq_symbol"])
names(sigs)=celltypes

sigs=lapply(sigs, function(x) x[x %in% rownames(cb)])

sampleTable$sample_type=plyr::revalue(as.factor(sampleTable$sample_type), c("E"="Endothelial", "F"="Fibroblast", "L"="Leukocyte", "T"="Malignant"))

df=as.data.frame(abind::abind(lapply(1:length(sigs), function(x) as.data.frame(cbind(exp=colMeans(cb[sigs[[x]],]), Cell.type=as.character(sampleTable$sample_type), sig=rep(names(sigs)[x], ncol(cb))))), along=1))

df$exp=as.numeric(as.character(df$exp))
df$Cell.type=as.factor(df$Cell.type)

df$sig
df$exp2=log2(df$exp)

df$sig=factor(df$sig, levels=c("Endothelial cells", "Epithelial cells", "Fibroblasts", "Mast cells", "B cells", "T cells", "Macrophages", "Dendritic cells"))

p=ggplot(df, aes(x=Cell.type, y=exp2, fill=Cell.type))+ geom_jitter(position=position_jitterdodge(jitter.width =0.2), size=0.3, alpha=0.8, pch=21) +geom_violin(aes(fill = factor(Cell.type)),alpha=0.6)+ facet_grid(~sig,scales="free", space="free_x") + scale_fill_manual(name = "Cell type", values = cbp[c(2,3,8,4)])+ ylab("Mean log2 expression")+ xlab("Cell type")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

file=paste0(FiguresdirU54, "violinplots.celltype.markers.Panglaodb_updated_070223")
pdf(file=paste0(file,'.pdf',sep=''), height = 5,  width = 9, family = "Helvetica")
p+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) + theme(legend.position = "none")
dev.off()

############################################################################
#Heatmap of CIBERSORT fractions in U54 bulk data (inferred from Puram)
#############################################################################

sampleTable=readRDS(paste0(DatadirU54Bulk, "U54_data_clinical_formatted_all.rds"))
sampleTable=sampleTable[sampleTable$sample_type %in% c("E","F","L","T") & sampleTable$cancertype=="HNSCC" & sampleTable$ln=="0",] %>% 
  dplyr::mutate(sample_type=plyr::revalue(as.factor(sample_type), c("E"="Endothelial", "F"="Fibroblast", "L"="Leukocyte", "T"="Malignant"))) %>% 
  arrange(sample_type)

cib=read.table(paste0(DatadirU54Bulk, "CIBERSORTx_Adjusted_docker_inferred_Puram.txt"), sep="\t", header=T) %>% 
  filter(!str_detect(Mixture, "MEL")) 
cib$Mixture=gsub("[.]","-", cib$Mixture)

cib=cib %>% 
  dplyr::filter(Mixture %in% sampleTable$filename)

cib=cib[match(sampleTable$filename, cib$Mixture),] %>% 
  filter(P.value<0.05)

#celltypes=c("T.cells.CD8", "T.cells.CD4", "Fibroblast",  "Macrophage",  "B.cell",   "Malignant",  "Mast",  "Dendritic", "Myocyte",  "Endothelial")
celltypes=c("T.cells.CD8", "T.cells.CD4", "Fibroblast",  "Macrophage",  "B.cell",   "Malignant",  "Dendritic", "Endothelial")

tab=cib %>% tibble::column_to_rownames("Mixture") %>% 
  dplyr::select(celltypes) %>% 
  rename_with(~ gsub(".", " ", .x, fixed = TRUE))

df=sampleTable %>% dplyr::select(sample_type) %>% 
  purrr::set_names("Flow-sorted cell type")

cbp=c("#000000", "#E69F00", "#56B4E9", "#009E73","#FD61D1", "#A3A500", "#D55E00", "#CC79A7")
cols.list=list(`Flow-sorted cell type`=cbp[c(2,3,8,4)] %>% purrr::set_names(levels(df$`Flow-sorted cell type`)))

rowannot=ComplexHeatmap::rowAnnotation(df = df, col = cols.list, show_annotation_name = TRUE)

p_cluster=tab %>% 
  as.matrix %>% 
  ComplexHeatmap::Heatmap(col=RColorBrewer::brewer.pal(9,"Reds"), show_row_names =F, name = "CIBERSORTx-inferred\ncell\nfraction", right_annotation = rowannot, clustering_method_rows = 'ward.D2',
                          column_title = "CIBERSORTx-inferred cell type", column_title_side = "bottom")

file=paste0(FiguresdirU54, "Heatmap_CIBERSORT_U54_bulk")
pdf(file=paste0(file,'.pdf',sep=''), height = 6,  width = 6, family = "Helvetica")
p_cluster
dev.off()






##########################################################################
#Plot pan-cancer meta-Z scores for survival and LNM-associated genes 
##########################################################################

fm=readRDS(paste0(Datadir, "Findmarkers_Puram_U54_allcombined.update.061323.rds"))

pc=read.table(paste0(Datadir, "PRECOG_original_Z_score_matrix.txt"), sep="\t", header=T)

pc.adv=as.character(pc[pc$Unweighted.meta.z..all.cancers.>=3.09,1])

pc.overlap=pc[match(fm$symbol_AnnotationDbi, pc$Gene.Symbol),]
fm$meta.z.original.PRECOG=pc.overlap$Unweighted.meta.z..all.cancers.

fm$survival.direction.original.PRECOG=rep(NA, nrow(fm))
fm$survival.direction.original.PRECOG[fm$meta.z.original.PRECOG>=3.09]="Anti-Survival"
fm$survival.direction.original.PRECOG[fm$meta.z.original.PRECOG<=(-3.09)]="Pro-Survival"

n=100
precog.adv=pc[order(pc$Unweighted.meta.z..all.cancers.,decreasing=T)[1:n],1]
precog.fav=pc[order(pc$Unweighted.meta.z..all.cancers.,decreasing=F)[1:n],1]

lnm.cols=gnuplot_colors[c(8:11,13:14)]
cols2=c(lnm.cols, "lightgrey")

options(ggrepel.max.overlaps = Inf) 

surv.cols=gnuplot_colors[c(1:5,7)]
cols2=c(surv.cols, "lightgrey")

dat=data.frame(Meta.z=fm$meta.z.original.PRECOG, Signature=factor(fm$lnm.gene.cluster.phenograph), gene=fm$symbol_AnnotationDbi)
dat=dat[!is.na(fm$lnm.gene.cluster.phenograph) & !is.na(dat$Meta.z),]
dat1=dat

dat=data.frame(Meta.z=fm$meta.z.original.PRECOG, Signature=factor(fm$lnm.direction), gene=fm$symbol_AnnotationDbi)
dat=dat[!is.na(fm$lnm.direction) & !is.na(dat$Meta.z),]
dat2=dat

dat=rbind(dat1, dat2)
dat$Signature=factor(dat$Signature, levels=c("Anti-LNM","Pro-LNM","L1", "L2", "L3", "L4", "L5", "L6"))

dat$Meta.z.rank=rep(1, nrow(dat))
dat$Meta.z.rank[which(dat$gene %in% precog.adv)]=2
dat$Meta.z.rank[which(dat$gene %in% precog.fav)]=3
dat$Meta.z.rank=factor(dat$Meta.z.rank, levels=c("2","3","1"))
dat.LNM=dat
dat.LNM$Outcome=rep("LNM", nrow(dat.LNM))

t1=ggplot(dat, aes(x=as.factor(Signature),y=Meta.z))
t1=t1+ geom_boxplot(outlier.shape=NA) + geom_point(size=0.8, position = position_jitter(w = 0.1, h = 0), aes(color = Meta.z.rank))
t1=t1+ geom_text_repel(aes(x=Signature, label=ifelse(Meta.z.rank!=1,as.character(gene),'')), vjust = 0.5, size=2.5)
p=t1 + scale_color_manual(values = c("red", "blue", "black"), labels=c("Highest 100 (Adverse)","Lowest 100 (Favorable)", "Not top rank"))+labs(color="Pan-cancer meta-z rank", x="Prognostic signature", y="Pan-cancer meta-z score")
p=p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p.LNM=p
p.LNM=p.LNM +ylim(-12,14)
p.LNM=p.LNM+ geom_hline(yintercept=0, linetype="dashed", color = "magenta")
p.LNM=p.LNM + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white", colour = "black",size = 0.5, linetype = "solid"))

#
dat=data.frame(Meta.z=fm$meta.z.original.PRECOG, Signature=factor(fm$survival.gene.cluster.phenograph), gene=fm$symbol_AnnotationDbi)
dat=dat[!is.na(fm$survival.gene.cluster.phenograph) & !is.na(dat$Meta.z),]
dat1=dat

dat=data.frame(Meta.z=fm$meta.z.original.PRECOG, Signature=factor(fm$survival.direction), gene=fm$symbol_AnnotationDbi)
dat=dat[!is.na(fm$survival.direction) & !is.na(dat$Meta.z),]
dat2=dat

dat=rbind(dat1, dat2)
dat$Signature=factor(dat$Signature, levels=c("Pro-survival", "Anti-survival", "S1", "S2", "S3", "S4", "S5", "S6"))

dat$Meta.z.rank=rep(1, nrow(dat))
dat$Meta.z.rank[which(dat$gene %in% precog.adv)]=2
dat$Meta.z.rank[which(dat$gene %in% precog.fav)]=3
dat$Meta.z.rank=factor(dat$Meta.z.rank, levels=c("2","3","1"))
dat.survival=dat
dat.survival$Outcome=rep("Survival", nrow(dat.survival))

t1=ggplot(dat, aes(x=Signature,y=Meta.z))
t1=t1+ geom_boxplot(outlier.shape=NA) + geom_point(size=0.8, position = position_jitter(w = 0.1, h = 0), aes(color = Meta.z.rank))
t1=t1+ geom_text_repel(aes(x=Signature, label=ifelse(Meta.z.rank!=1,as.character(gene),'')), vjust = 0.5, size=2.5)
p=t1 + scale_color_manual(values = c("red", "blue", "black"), labels=c("Highest 100 (Adverse)","Lowest 100 (Favorable)", "Not top rank"))+labs(color="Pan-cancer meta-z rank", x="Prognostic signature", y="Pan-cancer meta-z score")
p=p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p.survival=p
p.survival=p.survival +ylim(-12,14)

p.survival=p.survival + geom_hline(yintercept=0, linetype="dashed", color = "magenta")
p.survival=p.survival + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white", colour = "black",size = 0.5, linetype = "solid"))

library(grid)

p.survival.noleg=p.survival +theme(legend.position = "none")
p.LNM.noleg=p.LNM +theme(legend.position = "none")

file=paste0(Figuresdir, "jitter_clusters_LNM_and_survival_original_precog_genes_updated_061823")
pdf(file=paste0(file,'.pdf',sep=''), height = 8,  width = 12, family = "Helvetica")
grid.newpage()
grid.draw(cbind(ggplotGrob(p.survival.noleg), ggplotGrob(p.LNM.noleg), size = "last"))
dev.off()

#
legend <- cowplot::get_legend(p.survival+theme(legend.position = "bottom"))

file=paste0(Figuresdir, "jitter_clusters_LNM_and_survival_original_precog_genes_legend_updated_061823")
pdf(file=paste0(file,'.pdf',sep=''), height = 8,  width = 12, family = "Helvetica")
grid.newpage()
grid.draw(legend)
dev.off()

#Getting numbers of overlap with PRECOG pan-cancer survival genes 
tab=table(fm$survival.direction.original.PRECOG, fm$lnm.gene.cluster.phenograph)
surv.nums=as.numeric(summary(as.factor(fm$lnm.gene.cluster.phenograph))[1:6])

barplot(sweep(tab, MARGIN=2, surv.nums, "/"), beside=T, col=c(1:2))

barplot(tab, beside=T)

tab=table(fm$survival.direction.original.PRECOG, fm$survival.gene.cluster.phenograph)
surv.nums=as.numeric(summary(as.factor(fm$lnm.gene.cluster.phenograph))[1:6])

barplot(sweep(tab, MARGIN=2, surv.nums, "/"), beside=T, col=c(1:2))

##############################################################
##############################################################
#Association of HPV signatres with survival and LNM adjusted for HPV
##############################################################
##############################################################

RscriptsPath="/Users/kbren/Documents/scripts/rscripts/"
source(paste0(RscriptsPath, "basic_scripts.R"))
source("~/Documents/Projects/HNSCC_PreCog/data/Pecog_HNSCC_directories.R")

fm=readRDS(paste0(Datadir, "Findmarkers_Puram_U54_allcombined.update.061323.rds"))

makesig=function(genetype, group){
  levs=levels(fm[,group])
  sigs=lapply(levs, function(lev) 
  {x=fm[fm[,group]==lev,genetype]
  x=x[!is.na(x)]
  return(x)})
  names(sigs)=levs
  return(sigs)
}

sigs=
  c(makesig(genetype="gene",group="lnm.direction"),
    makesig(genetype="gene",group="lnm.gene.cluster.phenograph"),
    makesig(genetype="gene",group="survival.direction"),
    makesig(genetype="gene",group="survival.gene.cluster.phenograph"))

alldata2=readRDS("~/Documents/Projects/HNSCC_PreCog/data/Precog.HNSCC.all.exp.clin.data.curated.for.PRECOG.rds")

explist=alldata2$explist
infolist=alldata2$infolist

#restrict to datasets with >1000 genes 
gene.n.threshold=1000
accessions.survival=names(which(lapply(explist, nrow)>=gene.n.threshold))
infolist=infolist[accessions.survival]
explist=explist[accessions.survival]

sigslist=list()
#
for(i in 1:length(explist)){
  exp=explist[[i]]
  sigdf=as.data.frame(t(abind::abind(lapply(sigs, function(x) scale(colMeans(exp[rownames(exp) %in% x,], na.rm = T))), along=2)))
  sigslist[[i]]=sigdf
}
names(sigslist)=names(explist)

for(i in 1:length(sigslist)){
  print(all(names(sigslist[[i]])==colnames(explist[[i]])))
}

for(i in 1:length(sigslist)){
  print(all(names(sigslist[[i]])==rownames(infolist[[i]])))
}

alldata3=list(sigslist=sigslist, explist=explist, infolist=infolist)

saveRDS(alldata3, "~/Documents/Projects/HNSCC_PreCog/data/Precog.HNSCC.all.sig.exp.clin.data.updated.180623.rds")

################################
#Run meta-analysis for association of signatures with survival adjusting for HPV status
##############################

explist=sigslist

Survival_objects=readRDS(paste(Resultsdir, "Precog.HNSCC.survival.objects.all.HNSCC.rds", sep=""))
accessions.survival=intersect(accessions.survival, names(Survival_objects))

#19 objects
#Restrict to studies that were part of previous meta-analysis

allcoxph=readRDS(paste(Resultsdir, "Precog.HNSCC.coxph.survival.allstudies.allHNSCC.updated.06_09_2023.rds", sep=""))
accessions.survival=intersect(accessions.survival, names(allcoxph))
#length(accessions.survival) #16 studies. Just Chung removed

explist=explist[accessions.survival]
Survival_objects=Survival_objects[accessions.survival]

explist=lapply(names(Survival_objects), function(x) explist[[x]][,names(Survival_objects[[x]])])
names(explist)=names(Survival_objects)
infolist=lapply(names(Survival_objects), function(x) infolist[[x]][names(Survival_objects[[x]]),])
names(infolist)=names(Survival_objects)

numbers.HPV.surv=abind::abind(lapply(infolist, function(x) data.frame(Negative=length(which(x$COV_HNSCC_HPV_status=="negative")), Positive=length(which(x$COV_HNSCC_HPV_status=="positive")))), along=1)
numbers.HPV.surv=as.data.frame(numbers.HPV.surv)

HPV.studies=rownames(numbers.HPV.surv[which(numbers.HPV.surv[,1]>=10 & numbers.HPV.surv[,2]>=10),])
#Accession #"GSE39366" "GSE65858" "Thurlow"  "TCGA" 

#Have to find studies in which HPV status is mixed 
library(survival)
#Get genes associated with each signature 

genes=names(sigs)

accessions.survival=HPV.studies

genelist=list()
for(j in 1:length(genes)){
  gene=genes[j]
  
  coxres=list()
  for(i in 1:length(accessions.survival)){
    acc=accessions.survival[i]
    exp=as.matrix(explist[[acc]])
    info=infolist[[acc]]
    
    if(gene %in% rownames(exp)){
      cox=coxph(Survival_objects[[acc]]~exp[gene,] + info$COV_HNSCC_HPV_status)
      z.exp=as.list(coef(cox)/sqrt(diag(vcov(cox))))[[1]]
      z.HPV=as.list(coef(cox)/sqrt(diag(vcov(cox))))[[2]]
      se.exp=summary(cox)$coefficients[1,3]
      se.HPV=summary(cox)$coefficients[2,3]
      n = cox$n
      df.exp = data.frame(Z = z.exp, SE = se.exp, N = n, DatasetID = acc)
      df.HPV = data.frame(Z = z.HPV, SE = se.HPV, N = n, DatasetID = acc)
      df=list(exp=df.exp, hpv=df.HPV)
      coxres[[acc]]=df
    }
  }
  genelist[[gene]]=coxres
}

list.genes.metaframes=lapply(genelist, function(x) as.data.frame(abind::abind(lapply(x, function(y) y$exp), along=1)))
for(i in 1:length(list.genes.metaframes)){
  list.genes.metaframes[[i]]$Signature=rep(names(list.genes.metaframes)[i], nrow(list.genes.metaframes[[i]]))
}

#now apply liptak test 
metaz.list.exp.adj.HPV=unlist(lapply(list.genes.metaframes,function(mdf)sum(as.numeric(mdf$N)*as.numeric(mdf$Z))/sqrt(sum(as.numeric(mdf$N)^2))))
#This is just the association of HPV with survival adjusted for the gene expression signatures, so pretty useless
list.hpv.metaframes=lapply(genelist, function(x) as.data.frame(abind::abind(lapply(x, function(y) y$hpv), along=1)))

allresults=list(metaframes=list.genes.metaframes, meta.zs=metaz.list.exp.adj.HPV)

saveRDS(allresults, paste0(Resultsdir, "coxph.survival.adj.HPV.updated.LNM.z.score.061823.rds"))
allresults=readRDS(paste0(Resultsdir, "coxph.survival.adj.HPV.updated.LNM.z.score.061823.rds"))
list.genes.metaframes=allresults$metaframes
metaz.list.exp.adj.HPV=allresults$meta.zs

tab=as.data.frame(t(abind::abind(lapply(list.genes.metaframes, function(x) as.numeric(x$Z)), along=2)))
colnames(tab)=rownames(list.genes.metaframes[[1]])
tab=tab[9:nrow(tab),]

rownames(tab)=paste0(rownames(tab)," | Meta-Z score: ", round(metaz.list.exp.adj.HPV,2)[9:length(metaz.list.exp.adj.HPV)])

studymeta=readRDS(paste0(Datadir, "table_1_adding_authors_OPH_HPV_info.rds"))
studymeta=studymeta[match(colnames(tab), studymeta$Study.accession),]
colnames(tab)=studymeta$Study_label_author
tab=as.matrix(tab)

#
p=ComplexHeatmap::Heatmap(tab, cluster_rows = FALSE, cluster_columns = FALSE, name="Z-score", row_title = "Survival gene signatures", column_title = "Primary HNC patient studies", column_title_side = "bottom")
library(ComplexHeatmap)

file=paste(Figuresdir, "heatmap_meta_Z_survival_adj_HPV_061823")
pdf(file=paste0(file,'.pdf',sep=''), height = 4,  width = 6, family = "Helvetica")
ComplexHeatmap::draw(p, heatmap_legend_side="left", padding = unit(c(2, 2, 2, 20), "mm"))
dev.off()

saveRDS(tab, paste0(Figuresdir.grobs, "heatmap_meta_Z_survival.sourcedata_061823.rds"))

#############################
#Test association of LNM-associated signatures with LNM adjusted for HPV
#############################

alldata2=readRDS("~/Documents/Projects/HNSCC_PreCog/data/Precog.HNSCC.all.exp.clin.data.curated.for.PRECOG.rds")

explist=alldata2$explist
infolist=alldata2$infolist

#restrict to datasets with >1000 genes 
gene.n.threshold=1000
accessions.survival=names(which(lapply(explist, nrow)>=gene.n.threshold))
explist=explist[accessions.survival]

sigslist=list()
#
for(i in 1:length(explist)){
  exp=explist[[i]]
  sigdf=as.data.frame(t(abind::abind(lapply(sigs, function(x) scale(colMeans(exp[rownames(exp) %in% x,], na.rm = T))), along=2)))
  sigslist[[i]]=sigdf
}
names(sigslist)=names(explist)
explist=sigslist

#LNM
node.list=readRDS(paste(Resultsdir,"Precog.meta.analysis.node.status.all.HNSCC.updated.06.09.23.rds", sep=""))
accessions.node=intersect(names(explist), names(node.list$effect.sizes))

numbers.HPV.surv=abind::abind(lapply(infolist, function(x) data.frame(Negative=length(which(x$COV_HNSCC_HPV_status=="negative")), Positive=length(which(x$COV_HNSCC_HPV_status=="positive")))), along=1)
numbers.HPV.surv=as.data.frame(numbers.HPV.surv)
HPV.studies=rownames(numbers.HPV.surv[which(numbers.HPV.surv[,1]>=10 & numbers.HPV.surv[,2]>=10),])

accessions.node=intersect(accessions.node, HPV.studies)
accessions.node #"GSE33205" "GSE39366" "GSE65858" "Thurlow"  "TCGA"    

genelist=list()
for(j in 1:length(genes)){
  gene=genes[j]
  
  glmres=list()
  for(i in 1:length(accessions.node)){
    acc=accessions.node[i]
    exp=as.matrix(explist[[acc]])
    info=infolist[[acc]]
    
    if(gene %in% rownames(exp)){
      mod=glm(infolist[[acc]]$COVAR_N_status~as.matrix(explist[[acc]])[gene,] + infolist[[acc]]$COV_HNSCC_HPV_status, family = "binomial")
      m=summary(mod)
      z.exp=coef(m)[2,3]
      z.HPV=coef(m)[3,3]
      se.exp=coef(m)[2,2]
      se.HPV=coef(m)[3,2]
      n=stats::nobs(mod)
      df.exp = data.frame(Z = z.exp, SE = se.exp, N = n, DatasetID = acc)
      df.HPV = data.frame(Z = z.HPV, SE = se.HPV, N = n, DatasetID = acc)
      df=list(exp=df.exp, hpv=df.HPV)
      glmres[[acc]]=df
    }
  }
  genelist[[gene]]=glmres
}

list.genes.metaframes=lapply(genelist, function(x) as.data.frame(abind::abind(lapply(x, function(y) y$exp), along=1)))
for(i in 1:length(list.genes.metaframes)){
  list.genes.metaframes[[i]]$Signature=rep(names(list.genes.metaframes)[i], nrow(list.genes.metaframes[[i]]))
}

#now apply liptak test 
metaz.list.exp.node.adj.HPV=unlist(lapply(list.genes.metaframes,function(mdf)sum(as.numeric(mdf$N)*as.numeric(mdf$Z))/sqrt(sum(as.numeric(mdf$N)^2))))
#This is just the association of HPV with survival adjusted for the gene expression signatures, so pretty useless
list.hpv.node.metaframes=lapply(genelist, function(x) as.data.frame(abind::abind(lapply(x, function(y) y$hpv), along=1)))
#unlist(lapply(list.hpv.metaframes,function(mdf)sum(as.numeric(mdf$N)*as.numeric(mdf$Z))/sqrt(sum(as.numeric(mdf$N)^2))))

allresults=list(metaframes=list.genes.metaframes, meta.zs=metaz.list.exp.node.adj.HPV)

saveRDS(allresults, paste0(Resultsdir, "glm.node.adj.HPV.updated.06_09_2023.rds"))
allresults=readRDS(paste0(Resultsdir, "glm.node.adj.HPV.updated.06_09_2023.rds"))
list.genes.metaframes=allresults$metaframes
metaz.list.exp.node.adj.HPV=allresults$meta.zs

tab=as.data.frame(t(abind::abind(lapply(list.genes.metaframes, function(x) as.numeric(x$Z)), along=2)))
colnames(tab)=rownames(list.genes.metaframes[[1]])
rownames(tab)=paste0(rownames(tab)," | Meta-Z score: ", round(metaz.list.exp.node.adj.HPV,2))

saveRDS(tab, paste0(Figuresdir.grobs, "heatmap_meta_Z_node_LNM_adj_HPV_with_surv_updated.06_09_2023.sourcedata.rds"))
tab=readRDS(paste0(Figuresdir.grobs, "heatmap_meta_Z_node_LNM_adj_HPV_with_surv_updated.06_09_2023.sourcedata.rds"))


tab=as.matrix(tab)
p=ComplexHeatmap::Heatmap(tab, cluster_rows = FALSE, cluster_columns = FALSE, name="Z-score", row_title = "LNM signatures", column_title = "Primary HNC patient studies", column_title_side = "bottom")

tab=tab[1:8,]

#changing column names to study names
studymeta=readRDS(paste0(Datadir, "table_1_adding_authors_OPH_HPV_info.rds"))
studymeta=studymeta[match(colnames(tab), studymeta$Study.accession),]
colnames(tab)=studymeta$Study_label_author

p=ComplexHeatmap::Heatmap(tab, cluster_rows = FALSE, cluster_columns = FALSE, name="Z-score", row_title = "LNM signatures", column_title = "Primary HNC patient studies", column_title_side = "bottom")

file=paste(Figuresdir, "heatmap_meta_Z_node_LNM_adj_HPV_updated_061823")
pdf(file=paste0(file,'.pdf',sep=''), height = 3,  width = 6, family = "Helvetica")
ComplexHeatmap::draw(p, heatmap_legend_side="left", padding = unit(c(2, 2, 2, 20), "mm"))
dev.off()


##################################################################################
##################################################################################
#Making table indicating cell types in which each gene is highest expressed (Highest normalized counts)
##################################################################################
##################################################################################

fm=readRDS(paste0(Datadir, "Findmarkers_Puram_U54_allcombined.update.061323.rds"))

#Make table of genes
colnames=c(
  "Max_celltype_normalized_counts_Puram",
  "Max_celltype_normalized_counts_Stanford",
  "Max_celltype_normalized_counts_bulk")

list.nums=lapply(colnames, function(x) table(fm[,x], fm$lnm.gene.cluster.phenograph))
list.percentages=lapply(colnames, function(x) round(prop.table(table(fm[,x], fm$lnm.gene.cluster.phenograph),2)*100,1))
list.percentages.lnm.clusters=list.percentages

mp.l=matrix(paste(list.nums[[1]], " (", list.percentages[[1]],")", sep=""),7,6)
ms.l=matrix(paste(list.nums[[2]], " (", list.percentages[[2]],")", sep=""),7,6)
msb.l=matrix(paste(list.nums[[3]], " (", list.percentages[[3]],")", sep=""),4,6)

dimnames(mp.l)=dimnames(list.nums[[1]])
dimnames(ms.l)=dimnames(list.nums[[2]])
dimnames(msb.l)=dimnames(list.nums[[3]])

mp.l=as.data.frame(cbind("Study"=rep("Puram scRNA-Seq",nrow(mp.l)), "Cell type"=rownames(mp.l),mp.l))
ms.l=as.data.frame(cbind("Study"=rep("Stanford scRNA-Seq",nrow(ms.l)), "Cell type"=rownames(ms.l),ms.l))
msb.l=as.data.frame(cbind("Study"=rep("Stanford bulk RNA-Seq",nrow(msb.l)), "Cell type"=rownames(msb.l),msb.l))

list.nums=lapply(colnames, function(x) table(fm[,x], fm$lnm.direction))
list.percentages=lapply(colnames, function(x) round(prop.table(table(fm[,x], fm$lnm.direction),2)*100,1))
list.percentages.lnm.dir=list.percentages

mp.ld=matrix(paste(list.nums[[1]], " (", list.percentages[[1]],")", sep=""),7,2)
ms.ld=matrix(paste(list.nums[[2]], " (", list.percentages[[2]],")", sep=""),7,2)
msb.ld=matrix(paste(list.nums[[3]], " (", list.percentages[[3]],")", sep=""),4,2)

dimnames(mp.ld)=dimnames(list.nums[[1]])
dimnames(ms.ld)=dimnames(list.nums[[2]])
dimnames(msb.ld)=dimnames(list.nums[[3]])

mp.l=cbind(mp.l[,1:2], mp.ld, mp.l[,3:8])
ms.l=cbind(ms.l[,1:2], ms.ld, ms.l[,3:8])
msb.l=cbind(msb.l[,1:2], msb.ld, msb.l[,3:8])

list.percentages.lnm.all=list(cbind(list.percentages.lnm.dir[[1]],list.percentages.lnm.clusters[[1]]),
                              cbind(list.percentages.lnm.dir[[2]],list.percentages.lnm.clusters[[2]]),
                              cbind(list.percentages.lnm.dir[[3]],list.percentages.lnm.clusters[[3]]))


#survival genes
list.nums=lapply(colnames, function(x) table(fm[,x], fm$survival.gene.cluster.phenograph))
list.percentages=lapply(colnames, function(x) round(prop.table(table(fm[,x], fm$survival.gene.cluster.phenograph),2)*100,1))
list.percentages.survival.clusters=list.percentages

mp.s=matrix(paste(list.nums[[1]], " (", list.percentages[[1]],")", sep=""),7,6)
ms.s=matrix(paste(list.nums[[2]], " (", list.percentages[[2]],")", sep=""),7,6)
msb.s=matrix(paste(list.nums[[3]], " (", list.percentages[[3]],")", sep=""),4,6)

dimnames(mp.s)=dimnames(list.nums[[1]])
dimnames(ms.s)=dimnames(list.nums[[2]])
dimnames(msb.s)=dimnames(list.nums[[3]])

list.nums=lapply(colnames, function(x) table(fm[,x], fm$survival.direction))
list.percentages=lapply(colnames, function(x) round(prop.table(table(fm[,x], fm$survival.direction),2)*100,1))
list.percentages.survival.dir=list.percentages

list.percentages.survival.all=list(cbind(list.percentages.survival.dir[[1]],list.percentages.survival.clusters[[1]]),
                                   cbind(list.percentages.survival.dir[[2]],list.percentages.survival.clusters[[2]]),
                                   cbind(list.percentages.survival.dir[[3]],list.percentages.survival.clusters[[3]]))

mp.sd=matrix(paste(list.nums[[1]], " (", list.percentages[[1]],")", sep=""),7,2)
ms.sd=matrix(paste(list.nums[[2]], " (", list.percentages[[2]],")", sep=""),7,2)
msb.sd=matrix(paste(list.nums[[3]], " (", list.percentages[[3]],")", sep=""),4,2)

dimnames(mp.sd)=dimnames(list.nums[[1]])
dimnames(ms.sd)=dimnames(list.nums[[2]])
dimnames(msb.sd)=dimnames(list.nums[[3]])

mp.s=cbind(mp.sd, mp.s)
ms.s=cbind(ms.sd, ms.s)
msb.s=cbind(msb.sd, msb.s)

mp=cbind(mp.l[,1:2], mp.s, mp.l[,3:ncol(mp.l)])
ms=cbind(ms.l[,1:2], ms.s, ms.l[,3:ncol(ms.l)])
msb=cbind(msb.l[,1:2], msb.s, msb.l[,3:ncol(msb.l)])

ms.all=rbind(mp, ms, msb)

#Need to addbinary level survival and lnm associations
ms.all$`Cell type`=gsub("B ","B/Plasma ", gsub(" cell or ","/", gsub("[.]"," ", ms.all$`Cell type`)))

saveRDS(ms.all, paste0(Resultsdir, "numbers.percentage.genes.highest.expressed.celltypes.allstudies.updated_062223.rds"))

########

#Combining survival and LNM genes 
tp=as.data.frame(cbind(list.percentages.survival.all[[1]], list.percentages.lnm.all[[1]]))
tp$Cell.type=rownames(tp)
tp$Study=rep("Puram scRNA-Seq", nrow(tp))
ts=as.data.frame(cbind(list.percentages.survival.all[[2]], list.percentages.lnm.all[[2]]))
ts$Cell.type=rownames(ts)
ts$Study=rep("Stanford scRNA-Seq", nrow(ts))

tb=as.data.frame(cbind(list.percentages.survival.all[[3]], list.percentages.lnm.all[[3]]))
tb$Cell.type=rownames(tb)
tb$Study=rep("Stanford bulk RNA-Seq", nrow(tb))

table.percentages.all=rbind(tp, ts, tb)

saveRDS(table.percentages.all, paste0(Resultsdir, "percentage.genes.highest.expressed.celltypes.allstudies.updated_062223.rds"))
table.percentages.all=readRDS(paste0(Resultsdir, "percentage.genes.highest.expressed.celltypes.allstudies.updated_062223.rds"))

table.percentages.all$Cell.type=gsub("B ","B/Plasma ", gsub(" cell or ","/", gsub("[.]"," ", table.percentages.all$Cell.type)))

#Make image of table, which can be used in excel to make table
library(gt)
library(paletteer)
library(webshot)
library(dplyr)
#webshot::install_phantomjs()
library(forcats)

tb=table.percentages.all %>% 
  mutate(`Cell.type` = fct_relevel(Cell.type, c("Malignant","Endothelial","Fibroblast","B/Plasma cell", "Mast", "Myeloid", "T/NK cell"))) %>% 
  arrange(Cell.type) %>% group_by(Study) %>% gt(rowname_col = "Cell.type")%>% 
  data_color(columns = c(`Anti-survival`:`L6`),
             colors = scales::col_numeric(
               paletteer::paletteer_d(
                 palette = "ggsci::red_material") %>% as.character(),
               domain = c(0:100)))
tb=tb %>% 
  tab_spanner(label = "Survival signatures",
              columns = c(`Anti-survival`:`S6`))%>% 
  tab_spanner(label = "LNM signatures",
              columns = c(`Anti-LNM`:`L6`))

tb %>%
  gtsave(
    "table.highest.expression.celltype.table.gt.updated.062223.png", expand = 10,
    path = Resultsdir)

tb %>%
  gtsave(
    "table.highest.expression.celltype.table.gt.updated.062223.rtf", expand = 10,
    path = Resultsdir)

saveRDS(tb, paste0(Resultsdir, "table.highest.expression.celltype.table.gt.updated.062223.rds"))


######################################
######################################
#Heatmaps of survival and LNM-associated genes, figure 1
#########################################
######################################
######################################
######################################
#making elaborate survival gene heatmap with p-values 
#########################################
######################################

fm=readRDS(paste0(Datadir, "Findmarkers_Puram_U54_allcombined.update.061323.rds"))

#Identify non-HPV related studies
HPV.stat.studies=readRDS(paste0(Datadir, "HPV.status.studies.meta.analysis.rds"))

allcoxph=readRDS(paste0(Resultsdir, "Precog.HNSCC.progostic.genes.Liptak.metaz.updated.06_09_2023.rds"))

#This 
for(i in 1:length(allcoxph)){
  df=as.data.frame(allcoxph[[i]])
  df$gene=rownames(df)
  allcoxph[[i]]=df
}

library(abind)
allcoxph.abind=as.data.frame(abind(allcoxph, along=1))
rownames(allcoxph.abind)=NULL

allcoxph.abind$Z=as.numeric(as.character(allcoxph.abind$Z))
allcoxph.abind$N=as.numeric(as.character(allcoxph.abind$N))
MetaFrame=allcoxph.abind

#PRECOG.HNSCC.All=readRDS(paste0(Resultsdir,"Precog.HNSCC.progostic.genes.Liptak.metaz.updated.03_30_2020.rds"))
#MetaFrame=allcoxph.abind
metaz.coxph=readRDS(paste(Resultsdir, "Precog.HNSCC.coxph.Liptak.metaz.allHNSCC.updated.06_09_2023.rds", sep=""))

threshold=3.09
metaz.coxph.neg=metaz.coxph[as.numeric(metaz.coxph$metaz.coxph)>3.09,]
metaz.coxph.neg=metaz.coxph.neg[order(metaz.coxph.neg$metaz.coxph, decreasing = T),]

metaz.coxph.pos=metaz.coxph[metaz.coxph$metaz.coxph<(-3.09),]
metaz.coxph.pos=metaz.coxph.pos[order(metaz.coxph.pos$metaz.coxph),]

metaz.coxph.sig=rbind(metaz.coxph.pos,metaz.coxph.neg)
metaz.coxph.sig=metaz.coxph.sig[order(metaz.coxph.sig$metaz.coxph),]

genes=as.character(metaz.coxph.sig$gene)
studies=names(allcoxph)

zarray=array(NA, c(length(genes), length(studies)))
rownames(zarray)=genes
colnames(zarray)=studies

for(i in 1:length(studies)){
  study=studies[i]
  MF=MetaFrame[MetaFrame$DatasetID==study,]
  zarray[,i]=as.numeric(as.character(MF[match(rownames(zarray), MF$gene),"Z"]))
}

#Get first author names instead of accessions
studymeta=readRDS(paste0(Datadir, "table_1_adding_authors_OPH_HPV_info.rds"))
studymeta=studymeta[match(colnames(zarray), studymeta$Study.accession),]

colnames(zarray)=studymeta$Study_label_author
studymeta$HPV.stat.studies

studymeta=studymeta[order(studymeta$HPV.stat.studies),]
zarray=zarray[,studymeta$Study_label_author]

#Order columns by HPV status
#HPV.stat.studies.zarray=HPV.stat.studies[names(HPV.stat.studies) %in% colnames(zarray)]
library(ComplexHeatmap)

df1=data.frame(c(rep("Pro-survival", nrow(metaz.coxph.pos)), rep("Anti-survival", nrow(metaz.coxph.neg))))
rownames(df1)=rownames(zarray)
colnames(df1)="Survival direction"

genecols=c("red","cyan") #othe color combo on wiki
names(genecols)=c("Pro-survival", "Anti-survival")
genecols=list("Survival direction"=genecols)

rowannot=rowAnnotation(df = df1, col = genecols, show_annotation_name = TRUE)

#rownames(zarray)=fm[match(rownames(zarray), fm$gene),"hgnc_symbol"]
rownames(zarray)=ifelse(!is.na(fm[match(rownames(zarray), fm$gene),"hgnc_symbol"]),fm[match(rownames(zarray), fm$gene),"hgnc_symbol"],rownames(fm))

genes=c("TP53INP1", "TP73", "BTG3", "CD247", "CD19", "IL2RG", "CD2","CD3D","CD3E","CD3G","CD5","CD6","CD7","SNAI2", "ITGB6", "FN1","P4HA1", "P4HA2", "GAPDH", "ENO1","VEGFC","CXCL1","CCND1") 
genes=genes[genes %in% rownames(zarray)]     

geneannot = rowAnnotation(foo = anno_mark(at = match(genes, rownames(zarray)), labels = genes))

#make barplots showing the proportion of HPV as well as the subanatomic locations of the samples withn the studies 
Survival_objects=readRDS(paste(Resultsdir, "Precog.HNSCC.survival.objects.all.HNSCC.rds", sep=""))

alldata2=readRDS("~/Documents/Projects/HNSCC_PreCog/data/Precog.HNSCC.all.exp.clin.data.curated.for.PRECOG.rds")
infolist=alldata2$infolist

infolist=infolist[names(Survival_objects)]
infolist=lapply(names(Survival_objects), function(x) infolist[[x]][names(Survival_objects[[x]]),])
names(infolist)=names(Survival_objects)

tab=as.data.frame(t(abind::abind(lapply(infolist, function(x) table(x$COV_HNSCC_HPV_status)), along=2)))
tab$NAs=unlist(lapply(infolist, function(x) length(which(is.na(x$COV_HNSCC_HPV_status)))))
tab=prop.table(as.matrix(tab),1)
tab=tab[studymeta$Study.accession,]
tab=tab*100
#tab=as.data.frame(prop.table(as.matrix(tab),1))

tab2=as.data.frame(t(abind::abind(lapply(infolist, function(x) table(x$COV_HNSCC_major_anatomical_subdivision)), along=2)))
tab2$NAs=unlist(lapply(infolist, function(x) length(which(is.na(x$COV_HNSCC_major_anatomical_subdivision)))))
tab2=prop.table(as.matrix(tab2),1)
tab2=tab2[studymeta$Study.accession,]
tab2=tab2*100

all(studymeta$Study_label_author==colnames(zarray))

cbp.site=c(gnuplot_colors[3:6],"white")
cbp.HPV=c("grey","black","white")

ha=HeatmapAnnotation('HPV status'=anno_barplot(as.data.frame(tab), gp = gpar(fill = cbp.HPV)),
                     'Site'=anno_barplot(as.data.frame(tab2), gp = gpar(fill = cbp.site)),
                     height = unit(30,"mm"), gap = unit(2.3, "mm"))

#add sample numbers. More clearly identify HPV negative cases

cbp.site=c(gnuplot_colors[3:6],"white")
cbp.HPV=c("grey","black","white")

studymeta$HPV.stat.studies2=plyr::revalue(studymeta$HPV.stat.studies, c("All.negative"="No", "Includes.positive.OPH"="Yes", "Possibly.includes.positive.OPH"="Incomplete data"))
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7")


ha.test=HeatmapAnnotation('HPV status'=anno_barplot(as.data.frame(tab), gp = gpar(fill = cbp.HPV)),
                          'Site'=anno_barplot(as.data.frame(tab2), gp = gpar(fill = cbp.site)),
                          'Includes HPV+ve OPC'= studymeta$HPV.stat.studies2,
                          height = unit(40,"mm"), gap = unit(2.3, "mm"),
                          col=list('Includes HPV+ve OPC' = c("No" = cbp1[1], "Yes" = cbp1[2], "Incomplete data" = cbp1[3])), show_legend = c(TRUE, TRUE, FALSE))

col_fun3 = circlize::colorRamp2(c(min(metaz.coxph.sig$metaz.coxph), 0, max(metaz.coxph.sig$metaz.coxph)), c("green","black","red"))
rowannot.list.surv=list(`Meta-z-score`=col_fun3)

#meta-Z annotation with no legend
geneannot2=rowAnnotation(`Meta-z-score`=metaz.coxph.sig$metaz.coxph, foo = anno_mark(at = match(genes, rownames(zarray)), labels = genes), col=rowannot.list.surv, show_legend=FALSE)
p.surv=Heatmap(zarray, show_heatmap_legend = FALSE, cluster_rows = F, cluster_columns = F, show_row_names = F, left_annotation = rowannot, right_annotation=geneannot2, name="Cox regression\nz-score", row_title = "Survival-associated genes", column_title = "Primary HNC patient studies", column_title_side = "bottom", heatmap_legend_param = list(legend_direction="horizontal"), top_annotation = ha.test)

#Then make false legends and add to top 

zs=c(-10,-5,0,5,10)
ps=2*pnorm(-abs(zs))
ps=signif(ps, digits=3)
labs=paste0(zs," [",ps,"]")
val=zs
col_fun2 = circlize::colorRamp2(c(min(val), 0, max(val)), c("green","black","red"))
lgd = Legend(col_fun = col_fun2, title = "Meta-Z score\n[Corresponding p-value]",labels = labs, at=zs)

zs2=c(-4,-2,0,2,4)
ps2=2*pnorm(-abs(zs2))
ps2=signif(ps2, digits=3)
labs2=paste0(zs2," [",ps2,"]")

val2=zs2
col_fun3 = circlize::colorRamp2(c(min(val2), 0, max(val2)), c("blue","white","red"))
lgd2 = Legend(col_fun = col_fun3, title = "Cox regression Z-score\n[Corresponding p-value]",labels = labs2, at=zs2)
lgd_list=list(lgd2, lgd)

file=paste0(Figuresdir, "heatmap_survival_genes_Liptak_meta_Z_complexheatmap_with_studyannotation_ordered_HPV_allgenes_seperate_legends_with_metaz_with_pvals_updated_062623")
pdf(file=paste0(file,'.pdf',sep=''), height = 10,  width = 6, family = "Helvetica")
draw(p.surv,  annotation_legend_list = lgd_list, heatmap_legend_side='top', annotation_legend_side='top')
decorate_annotation("HPV status", {
  grid.text("% HNCs", unit(-7, "mm"), just = "bottom", rot = 90)
})
decorate_annotation("Site", {
  grid.text("% HNCs", unit(-7, "mm"), just = "bottom", rot = 90)
})
dev.off()

#################################
#Extra annotation legends needed for both survival and LNM 
##################################

lgd.hpv = Legend(labels = c("Negative","Positive"), title = "HPV status", legend_gp = gpar(fill = cbp.HPV),title_position = "topcenter", direction = "horizontal")
lgd.site = Legend(labels = c("Hypopharynx","Larynx","Oral cavity","Oropharynx"), title = "Subanatomic site", legend_gp = gpar(fill = cbp.site),title_position = "topcenter", direction = "horizontal")
lgd.hpv.oph = Legend(labels = levels(studymeta$HPV.stat.studies2), title = "Includes\nHPV+ve OPC", legend_gp = gpar(fill = cbp1[1:3]),title_position = "topcenter", direction = "horizontal")

lgd_list=list(lgd.hpv, lgd.site, lgd.hpv.oph)

#Make blank heatmap with shared legends
p.blank=Heatmap(matrix(nr=0, nc=ncol(zarray)))

file=paste0(Figuresdir, "heatmap_effectsizes_legends_updated_062623")
pdf(file=paste0(file,'.pdf',sep=''), height = 6,  width = 6, family = "Helvetica")
draw(p.blank, annotation_legend_list = lgd_list, annotation_legend_side="left", padding = unit(c(2, 10, 2, 2), "mm"))
dev.off()

######################################
######################################
#LNM heatmap with p-values corresponding to meta-z scores 
#########################################
######################################

allgenes2=readRDS(paste0(Datadir, "all.genes.survival.lnm.060923.rds"))

threshold=3.09
metaz.pro.lnm=allgenes2[which(as.numeric(allgenes2$metaz.lnm)>=3.09),]
metaz.pro.lnm=metaz.pro.lnm[order(metaz.pro.lnm$metaz.lnm, decreasing = T),]

metaz.anti.lnm=allgenes2[which(allgenes2$metaz.lnm<=(-3.09)),]
metaz.anti.lnm=metaz.anti.lnm[order(metaz.anti.lnm$metaz.lnm),]

metaz.lnm.sig=rbind(metaz.anti.lnm, metaz.pro.lnm)
metaz.lnm.sig=metaz.lnm.sig[order(metaz.lnm.sig$metaz.lnm),]
genes=as.character(metaz.lnm.sig$gene)

allgenes2[match(genes, allgenes2$gene),"symbol_AnnotationDbi"]

node.list=readRDS(paste0(Resultsdir,"Precog.meta.analysis.node.status.all.HNSCC.updated.06.09.23.rds"))

meta1=node.list$meta
metatab=cbind(
  unlist(lapply(meta1, function(x) x$pval)),
  unlist(lapply(meta1, function(x) x$b)))

colnames(metatab)=c("pval","RE model b")
metatab=as.data.frame(metatab)
metatab=metatab[order(metatab$pval),]
metatab$qval=p.adjust(metatab$pval)

#metatab=metatab[genes,]
#restrict to genes with absolute meta-Z >=3.09

pvals=metatab
colnames(pvals)[2]="beta"

siggenes=genes

#There are no z-scores for the mean-difference-based method. Need to just get a table of effect-sizes that are used to calculate random effects model p-value
effect.sizes=node.list$effect.sizes
effect.sizes.s=lapply(effect.sizes, function(x) x[names(x) %in% siggenes])

studies=names(effect.sizes)

zarray=array(NA, c(length(siggenes), length(studies)))
rownames(zarray)=siggenes
colnames(zarray)=studies

for(i in 1:length(studies)){
  study=studies[i]
  x=effect.sizes[[study]]
  es=unlist(lapply(x, function(y) y$es))
  es=es[match(siggenes, names(es))]
  print(all(names(es)==rownames(zarray), na.rm=T))
  zarray[,study]=es
}

#Removing GSE686, which doesn't have data for any signficant gene 
zarray=zarray[,names(which(apply(zarray, 2, function(x) !all(is.na(x)))))]

fm=readRDS(paste0(Datadir, "Findmarkers_Puram_U54_allcombined.update.061323.rds"))
metazs=fm[match(rownames(zarray), fm$gene),"metaz.lnm"]

#order zarray by HPV status
studymeta=readRDS(paste0(Datadir, "table_1_adding_authors_OPH_HPV_info.rds"))
studymeta=studymeta[match(colnames(zarray), studymeta$Study.accession),]

colnames(zarray)=studymeta$Study_label_author
studymeta$HPV.stat.studies

studymeta=studymeta[order(studymeta$HPV.stat.studies),]
zarray=zarray[,studymeta$Study_label_author]

#
val=ifelse(unlist(lapply(rownames(zarray), function(x) pvals[x,"beta"]))<0,"Anti-LNM","Pro-LNM")
df1=data.frame(val)
rownames(df1)=rownames(zarray)
colnames(df1)="LNM direction"

genecols=c("green","magenta")
names(genecols)=c("Anti-LNM", "Pro-LNM")
genecols=list("LNM direction"=genecols)

rowannot=rowAnnotation(df = df1, col = genecols, show_annotation_name = TRUE)

#make barplots showing the proportion of HPV as well as the subanatomic locations of the samples withn the studies 
LNM_objects=readRDS(paste0(Datadir, "Precog.HNSCC.lymph.node.metastasis.objects.all.HNSCC.rds"))
LNM_objects$`E-MTAB-1328`$COVAR_N_status

alldata2=readRDS("~/Documents/Projects/HNSCC_PreCog/data/Precog.HNSCC.all.exp.clin.data.curated.for.PRECOG.rds")
infolist=alldata2$infolist

infolist=infolist[names(LNM_objects)]
infolist=lapply(names(LNM_objects), function(x) infolist[[x]][rownames(LNM_objects[[x]]),])
names(infolist)=names(LNM_objects)

tab=as.data.frame(t(abind::abind(lapply(infolist, function(x) table(x$COV_HNSCC_HPV_status)), along=2)))
tab$NAs=unlist(lapply(infolist, function(x) length(which(is.na(x$COV_HNSCC_HPV_status)))))
tab=prop.table(as.matrix(tab),1)
tab=tab[studymeta$Study.accession,]
tab=tab*100
#tab=as.data.frame(prop.table(as.matrix(tab),1))

tab2=as.data.frame(t(abind::abind(lapply(infolist, function(x) table(x$COV_HNSCC_major_anatomical_subdivision)), along=2)))
tab2$NAs=unlist(lapply(infolist, function(x) length(which(is.na(x$COV_HNSCC_major_anatomical_subdivision)))))
tab2=prop.table(as.matrix(tab2),1)
tab2=tab2[studymeta$Study.accession,]
tab2=tab2*100

cbp.site=c(gnuplot_colors[3:6],"white")
cbp.HPV=c("grey","black","white")

studymeta$HPV.stat.studies2=plyr::revalue(studymeta$HPV.stat.studies, c("All.negative"="No", "Includes.positive.OPH"="Yes", "Possibly.includes.positive.OPH"="Incomplete data"))

#Need to add examples of important genes 
rownames(zarray)=allgenes2[match(rownames(zarray), allgenes2$gene),"hgnc_symbol"]

#Get original PRECOG and select most significant genes
pc=read.table(paste0(Datadir, "PRECOG_original_Z_score_matrix.txt"), sep="\t", header=T)

pc.adv=pc[order(pc$Unweighted.meta.z..all.cancers., decreasing=T),"Gene.Symbol"]

pc.fav=pc[pc$Unweighted.meta.z..all.cancers.<(-3.09),"Gene.Symbol"]
intersect(pc.fav, metaz.anti.lnm[1:50,"symbol_AnnotationDbi"])
#Only "SH3D19" overlapping

ha.test.lnm=HeatmapAnnotation('HPV status'=anno_barplot(as.data.frame(tab), gp = gpar(fill = cbp.HPV)),
                              'Site'=anno_barplot(as.data.frame(tab2), gp = gpar(fill = cbp.site)),
                              'Includes HPV+ve OPC'= studymeta$HPV.stat.studies2,
                              height = unit(40,"mm"), gap = unit(2.3, "mm"),
                              col=list('Includes HPV+ve OPC' = c("No" = cbp1[1], "Yes" = cbp1[2], "Incomplete data" = cbp1[3])), show_legend = c(TRUE, TRUE, FALSE))


#p.lnm=Heatmap(zarray, cluster_rows = F, cluster_columns = F,  name="Scaled mean\nexpression difference", row_title = "LNM-associated genes", column_title = "Primary HNC patient studies", show_row_names = F, left_annotation = rowannot, column_title_side = "bottom", top_annotation = ha,  heatmap_legend_param = list(legend_direction="horizontal"))
lgd.hpv = Legend(labels = c("Negative","Positive"), title = "HPV status", legend_gp = gpar(fill = cbp.HPV),title_position = "topcenter", direction = "horizontal")
lgd.site = Legend(labels = c("Hypopharynx","Larynx","Oral cavity","Oropharynx"), title = "Subanatomic site", legend_gp = gpar(fill = cbp.site),title_position = "topcenter", direction = "horizontal")
lgd.hpv.oph = Legend(labels = levels(studymeta$HPV.stat.studies2), title = "Includes HPV+ve OPH", legend_gp = gpar(fill = cbp1[1:3]),title_position = "topcenter", direction = "horizontal")
#lgd.hpv.oph = Legend(labels = levels(studymeta$HPV.stat.studies2), title = "Includes HPV+ve OPH", legend_gp = gpar(fill = cbp1[1:3]),title_position = "topcenter", direction = "horizontal")

lgd_list=list(lgd.hpv, lgd.site, lgd.hpv.oph)

##
col_fun3 = circlize::colorRamp2(c(min(metazs), 0, max(metazs)), c("green","black","red"))
rowannot.list.lnm=list(`Meta-z-score`=col_fun3)

genes=c("CDK1", "CCNB2", "CHEK2", "AURKA", "AURKB", "MKI67","MCM3", "MCM7", "MCM8", "MCM6", "POLA1", "RAD51", "MSH6","IVL", "KRT2", "KRT6B", "KRT10", "KRT16", "KRT23", "KRT75", "KRT78", "KRT80", "KRTDAP", "CDKN1A") 
geneannot2=rowAnnotation(`Meta-z-score`=metazs, foo = anno_mark(at = match(genes, rownames(zarray)), labels = genes), col=rowannot.list.lnm, show_legend=FALSE)

p.lnm=Heatmap(zarray, show_heatmap_legend = FALSE, cluster_rows = F, cluster_columns = F,  name="Scaled mean\nexpression difference", row_title = "LNM-associated genes", column_title = "Primary HNC patient studies", show_row_names = F, left_annotation = rowannot, column_title_side = "bottom", top_annotation = ha.test.lnm,  heatmap_legend_param = list(legend_direction="vertical"), right_annotation = geneannot2)
p.lnm

###########
zs=c(-10,-5,0,5,10)
ps=2*pnorm(-abs(zs))
ps=signif(ps, digits=3)
labs=paste0(zs," [",ps,"]")
val=zs
col_fun2 = circlize::colorRamp2(c(min(val), 0, max(val)), c("green","black","red"))
lgd = Legend(col_fun = col_fun2, title = "Meta-Z score\n[Corresponding p-value]",labels = labs, at=zs)

#
diffs=zs=c(-2,-1,0,1,2)
val2=diffs
col_fun3 = circlize::colorRamp2(c(min(val2), 0, max(val2)), c("blue","white","red"))
lgd2 = Legend(col_fun = col_fun3, title = "Scaled mean\nexpression difference",labels = val2, at=val2)
lgd_list=list(lgd2, lgd)

#
file=paste0(Figuresdir, "heatmap_effectsizes_node_complexheatmap_with_studyannotation_updated_101621_ordered_HPV_allgenes_seperate_legends_with_metaz_with_pvals_updated_062623")
pdf(file=paste0(file,'.pdf',sep=''), height = 10,  width = 6, family = "Helvetica")
draw(p.lnm,  annotation_legend_list = lgd_list, heatmap_legend_side='top', annotation_legend_side='top')
decorate_annotation("HPV status", {
  grid.text("% HNCs", unit(-7, "mm"), just = "bottom", rot = 90)
})
decorate_annotation("Site", {
  grid.text("% HNCs", unit(-7, "mm"), just = "bottom", rot = 90)
})
dev.off()





##############################################################################
##############################################################################
#Making supplementary table 2. Table of prognostic genes 
##############################################################################
##############################################################################

fm=readRDS(paste0(Datadir, "Findmarkers_Puram_U54_allcombined.update.061323.rds"))

#Add meta-z scores for associations of genes with grade
metaz.glm.annot=readRDS(paste(Resultsdir, "Precog.HNSCC.glm.grade.Liptak.metaz.allHNSCC.updated.061623.rds", sep=""))
metaz.glm.annot=metaz.glm.annot[match(fm$gene, metaz.glm.annot$gene),]
metaz.glm.annot$grade.direction=rep(NA, nrow(metaz.glm.annot))

metaz.glm.annot$grade.direction[metaz.glm.annot$signficant=="yes" & metaz.glm.annot$metaz.glm >0]="Pro-grade"
metaz.glm.annot$grade.direction[metaz.glm.annot$signficant=="yes" & metaz.glm.annot$metaz.glm <0]="Anti-grade"

fm$grade.meta.z=metaz.glm.annot$metaz.glm
fm$grade.direction=metaz.glm.annot$grade.direction

#Add original PRECOG data 
pc=read.table(paste0(Datadir, "PRECOG_original_Z_score_matrix.txt"), sep="\t", header=T)

pc.overlap=pc[match(fm$symbol_AnnotationDbi, pc$Gene.Symbol),]
fm$meta.z.original.PRECOG=pc.overlap$Unweighted.meta.z..all.cancers.

fm$survival.direction.original.PRECOG=rep(NA, nrow(fm))
fm$survival.direction.original.PRECOG[fm$meta.z.original.PRECOG>=3.09]="Anti-Survival"
fm$survival.direction.original.PRECOG[fm$meta.z.original.PRECOG<=(-3.09)]="Pro-Survival"

#############################################
#Add cell clusters with highest expression of prognostic genes in Puram and Stanford
###########################################

int=readRDS(paste0(DatadirPuram, "Integated_Puram_prim_patients200cells_update_160623.rds"))
DefaultAssay(int)="RNA"

keep=colnames(int[,int$cell.type.collapsed!="Unclassified" & int$cell.type.collapsed!="myocyte"])
int=subset(int, cells=keep)

#cr=int@assays$RNA@counts #raw counts
cn=int@assays$RNA@data #normalizard counts

fm=readRDS(paste0(Datadir, "Findmarkers_Puram_U54_allcombined.update.061323.rds"))

ga=cn
genes=fm$Matched_to_Puram[fm$Matched_to_Puram %in% rownames(ga)]
ga=ga[genes,]

cellclusts=paste0("cluster_", levels(int$FindClusters.res0.8))

meantab=apply(ga, 1, function(x) aggregate(x~int$FindClusters.res0.8, FUN=mean)[,2])
meantab=as.data.frame(t(meantab))
colnames(meantab)=cellclusts

colnames(meantab)=paste0("Puram_mean_normalized_counts_", colnames(meantab))
fmcl=meantab
cols=colnames(fmcl)

fmcl$Max_cellclusts_normalized_counts_Puram=apply(fmcl[,cols],1, function(x) cellclusts[which.max(x)])
fmcl$Max_cellclusts_normalized_counts_Puram=as.character(fmcl$Max_cellclusts_normalized_counts_Puram)

fmcl$Max_cellclusts_normalized_counts_Puram=factor(plyr::revalue(as.factor(fmcl$Max_cellclusts_normalized_counts_Puram), c("character(0)"=NA)), levels=paste0("cluster_",seq(0,14)))

library(dplyr)
fmcl=fmcl %>% tibble::rownames_to_column("Matched_to_Puram")
fmcl=fmcl[match(fm$Matched_to_Puram, fmcl$Matched_to_Puram),]
fm$`Cell cluster highest mean expression Puram scRNA-Seq`=as.character(gsub("cluster_", "", fmcl$Max_cellclusts_normalized_counts_Puram))

#Add highest cell cluster for Stanford scRNA-Seq dataset
int=readRDS(paste0(DatadirCCSBSc, "CCSB_scRNASeq_HNSCC_primary_enzymatic_intgrated_50PCs.mito.removed_update_160623.rds"))
DefaultAssay(int)="RNA"

keep=colnames(int[,int$cell.type.collapsed!="Unclassified" & int$cell.type.collapsed!="myocyte"])
int=subset(int, cells=keep)

#cr=int@assays$RNA@counts #raw counts
cn=int@assays$RNA@data #normalizard counts

ga=cn
genes=fm$Matched_to_Stanford[fm$Matched_to_Stanford %in% rownames(ga)]
ga=ga[genes,]

cellclusts=paste0("cluster_", levels(int$FindClusters.res0.8))

meantab=apply(ga, 1, function(x) aggregate(x~int$FindClusters.res0.8, FUN=mean)[,2])
meantab=as.data.frame(t(meantab))
colnames(meantab)=cellclusts

colnames(meantab)=paste0("Stanford_mean_normalized_counts_", colnames(meantab))
fmcl=meantab
cols=colnames(fmcl)

fmcl$Max_cellclusts_normalized_counts_Stanford=apply(fmcl[,cols],1, function(x) cellclusts[which.max(x)])
fmcl$Max_cellclusts_normalized_counts_Stanford=as.character(fmcl$Max_cellclusts_normalized_counts_Stanford)

fmcl=fmcl %>% tibble::rownames_to_column("Matched_to_Stanford")
fmcl=fmcl[match(fm$Matched_to_Stanford, fmcl$Matched_to_Stanford),]
fm$`Cell cluster highest mean expression Stanford scRNA-Seq`=as.character(gsub("cluster_", "", fmcl$Max_cellclusts_normalized_counts_Stanford))

#saveRDS(fm, paste0(Datadir, "Findmarkers_Puram_U54_allcombined.update.061323.rds"))

####################################
#Making supp table 2
####################################

fm$p.value.coxph=2*pnorm(-abs(fm$metaz.coxph))
fm$p.value.lnm=2*pnorm(-abs(fm$metaz.lnm))

colnames1=c(
  "gene",
  "symbol_AnnotationDbi", 
  "metaz.coxph",
  "p.value.coxph",
  "survival.direction",
  "survival.gene.cluster.phenograph",
  "metaz.lnm",
  "p.value.lnm",
  "lnm.direction",
  "lnm.gene.cluster.phenograph",
  "grade.meta.z",
  "grade.direction",
  "meta.z.original.PRECOG",
  "survival.direction.original.PRECOG",
  "Max_celltype_normalized_counts_Stanford",
  "Max_celltype_normalized_counts_Puram",
  "Max_celltype_normalized_counts_bulk",
  "Cell cluster highest mean expression Puram scRNA-Seq",  
  "Cell cluster highest mean expression Stanford scRNA-Seq")

all(colnames1 %in% colnames(fm)) #TRUE

tab=fm[,colnames1]

colnames=c("entrezgene ID",	"Gene symbol",	
           "Survival meta-z score",	"Survival p-value", "Survival directon of association", "Survival gene cluster",
           "LNM meta-z score",	"LNM p-value", "LNM directon of association", "LNM gene cluster",
           "Grade meta-z score", "Grade directon of association", "Survival pan-cancer (PRECOG) meta-z score", "Survival pan-cancer (PRECOG) directon of association",
           "Cell type highest mean expression Stanford scRNA-Seq", "Cell type highest mean expression Puram scRNA-Seq", "Cell type highest mean expression Stanford bulk RNA-Seq",
           "Cell cluster highest mean expression Puram scRNA-Seq",  "Cell cluster highest mean expression Stanford scRNA-Seq")
cbind(colnames, colnames(tab))

colnames(tab)
colnames(tab)=colnames

tab$`Survival p-value`=signif(tab$`Survival p-value`, digits=3)
tab$`LNM p-value`=signif(tab$`LNM p-value`, digits=3)
tab$`Survival meta-z score`=round(tab$`Survival meta-z score`,3)
tab$`LNM meta-z score`=round(tab$`LNM meta-z score`,3)
tab$`Grade meta-z score`=round(tab$`Grade meta-z score`,3)
tab$`Survival pan-cancer (PRECOG) meta-z score`=round(tab$`Survival pan-cancer (PRECOG) meta-z score`,3)

tab$`Cell type highest mean expression Stanford scRNA-Seq`=as.character(plyr::revalue(tab$`Cell type highest mean expression Stanford scRNA-Seq`, c("B.cell"="B/Plasma cell", "T.cell.or.NK.cell"="T/NK cell")))
tab$`Cell type highest mean expression Puram scRNA-Seq`=as.character(plyr::revalue(tab$`Cell type highest mean expression Puram scRNA-Seq`, c("B.cell"="B/Plasma cell", "T.cell.or.NK.cell"="T/NK cell")))
tab$`Cell type highest mean expression Stanford bulk RNA-Seq`=as.character(plyr::revalue(tab$`Cell type highest mean expression Stanford bulk RNA-Seq`, c("B.cell"="B/Plasma cell", "T.cell.or.NK.cell"="T/NK cell")))
tab$`Survival pan-cancer (PRECOG) directon of association`=gsub("Survival","survival",tab$`Survival pan-cancer (PRECOG) directon of association`)
#tab=readRDS(paste0(Resultsdir, "Genes.survival.lnm.genes.v2_updated_080823.rds"))

tab$`Survival meta-z score HPV-ve`<-fm$metaz.coxph.HPV.negative
tab$`LNM meta-z score HPV-ve`<-fm$metaz.lnm.HPV.negative

saveRDS(tab, paste0(Resultsdir, "Genes.survival.lnm.genes.v2_updated_080823.rds"))
write.table(tab, file=paste0(Resultsdir, "Genes.survival.lnm.genes.v2_updated_080823.txt"), sep="\t", row.names = FALSE, na="", quote=F)

#################################################################################
#################################################################################
#Meta analysis of genes associated with LNM status
######################################################################
#####################################################################

#node objects
#get gene expression and clinical data for all studies
alldata=readRDS(paste(Datadir, "Precog.HNSCC.all.exp.clin.data.rds", sep=""))
infolist=alldata$infolist
explist=alldata$explist

#
accessions=names(explist)
excluded_studies=c("GSE30788_discovery", "GSE30788_validation", "GSE686")
accessions=setdiff(accessions, excluded_studies)

#
#get lymph node objects
node_objects=readRDS(paste(Datadir, "Precog.HNSCC.lymph.node.metastasis.objects.all.HNSCC.rds", sep=""))

accessions.expnode=intersect(accessions, names(node_objects))
#20 datasets

#####restrictong further to datasets with at least five pos and five neg patients with gene exprssion data 
nodelengths=list()
for(i in 1:length(accessions.expnode)){
  acc=accessions.expnode[i]
  node_object=node_objects[[acc]]
  exp=as.data.frame(explist[[acc]])
  OverlapSamples=intersect(rownames(node_object), colnames(exp))
  node_object=node_object[OverlapSamples,]
  nodelengths[[acc]]=length(which(node_object$COVAR_N_status=="0"))>=5 & length(which(node_object$COVAR_N_status=="1"))>=5
}
accessions.expnode=names(which(unlist(nodelengths)))
#20 datasets left

################################################
#make matching lists of expression and node objects
###################################################

varlist=node_objects
varfactor="COVAR_N_status"

explist.node=list()
varlist_objects=list()

accessions.var=accessions.expnode

for(i in 1:length(accessions.var)){
  acc=accessions.var[i]
  exp.var=explist[[acc]]
  var.exp=varlist[[acc]]
  
  OverlapSamples=intersect(colnames(exp.var), rownames(var.exp))
  exp.var=exp.var[,OverlapSamples]
  var.exp=var.exp[OverlapSamples,]
  
  explist.node[[acc]]=exp.var
  varlist_objects[[acc]]=var.exp
}

#restrict to studies with gene in half or more studies. Otherwise top genes are those that are only in one ot two studies 
allgenes=unique(unlist(lapply(explist.node, rownames)))
#For each gene, find the number of studies the gene is in 
allgeneinstances=unlist(lapply(explist.node, rownames))

#Find the number of studies for each gene is available 
nstudies=lapply(allgenes, function(x) length(grep(paste("\\b",x, "\\b", sep=""), allgeneinstances)))
names(nstudies)=allgenes
#check a few randomly and they're all fine
saveRDS(nstudies, paste0(Datadir, "nstudies.ge.updated.06.09.23.rds"))
nstudies=readRDS(paste0(Datadir, "nstudies.ge.updated.06.09.23.rds"))

num=round(length(accessions.expnode)*.5) #10
genes.available.all=names(unlist(nstudies)[unlist(nstudies)>=num])
length(genes.available.all) #17264 genes in at least half (n=10) studies 
saveRDS(genes.available.all, paste(Resultsdir,"Precog.meta.analysis.genes.all.HNSCC.updated.06.09.23.rds", sep=""))
genes.available.all=readRDS(paste(Resultsdir,"Precog.meta.analysis.genes.all.HNSCC.updated.06.09.23.rds", sep=""))

#genes.available.all=readRDS(paste(Resultsdir,"Precog.meta.analysis.genes.all.HNSCC.rds", sep=""))

source(paste(RscriptsPath, "meta.analysis.mean.difference.functions.R", sep=""))
node.list=make.meta.node_multistudies(expressionlist=explist.node, studyaccessions=accessions.var, clinlist=varlist_objects, genes=genes.available.all, varfactor="COVAR_N_status")

lapply(explist.node, ncol) %>% unlist() %>% sum()


saveRDS(node.list, paste(Resultsdir,"Precog.meta.analysis.node.status.all.HNSCC.updated.06.09.23.rds", sep=""))
lnm.list=readRDS(paste(Resultsdir,"Precog.meta.analysis.node.status.all.HNSCC.updated.06.09.23.rds", sep=""))

#
meta2=lnm.list$meta

metatab=cbind(
  unlist(lapply(meta2, function(x) x$pval)),
  unlist(lapply(meta2, function(x) x$b)),
  unlist(lapply(meta2, function(x) x$zval)))

metatab=as.data.frame(metatab)
dim(metatab[metatab$V3>3.09,])
dim(metatab[metatab$V3<(-3.09),])
colnames(metatab)=c("pval","beta","zval")
metatab$gene=as.character(rownames(metatab))

####
#metatab=metatab[abs(metatab$zval)>=3.09,]

library(biomaRt)
ensembl=useMart("ensembl", dataset="hsapiens_gene_ensembl")

bm.gene=getBM(attributes=c("hgnc_symbol","entrezgene_id"), filters = "entrezgene_id", values = as.character(rownames(metatab)), mart= ensembl)
bm.gene=bm.gene[match(rownames(metatab), bm.gene$entrezgene_id),]

metatab=cbind(metatab, bm.gene)
metatab=metatab[order(abs(metatab$zval), decreasing = T),]
metatab$gene=rownames(metatab)
metatab$lnm_direction=rep(NA, nrow(metatab))
metatab$lnm_direction[metatab$zval<0]="Anti-LNM"
metatab$lnm_direction[metatab$zval>0]="Pro-LNM"
##

saveRDS(metatab, paste0(Datadir, "zscores_lnm_all_genes_updated.06.09.23.rds"))
metatab=readRDS(paste0(Datadir, "zscores_lnm_all_genes_updated.06.09.23.rds"))

##########################################################################
#Get lnm gene Phenograph clusters
##########################################################################

alldata=readRDS("~/Documents/Projects/HNSCC_PreCog/data/Precog.HNSCC.all.exp.clin.data.rds")
explist=alldata$explist
zscores.list=explist
infolist=alldata$infolist

#
allgenes=lapply(explist, rownames)
#for each dataset get the intersection with the signficant gene list 
allgenes2=lapply(allgenes, function(x) intersect(metatab$gene, x))

#somwhow make a dataframe showing presence or absence in each dataset 
genes=unique(unlist(allgenes2))
length(genes) #877 genes altogether 

genearray=array(NA, c(length(allgenes2), length(genes)))
rownames(genearray)=names(allgenes2)
colnames(genearray)=genes
genearray=as.data.frame(genearray)

gene.present=lapply(allgenes2, function(x)  ifelse(genes %in% x, "yes","no"))
library(abind)
gene.present.df=as.data.frame(abind(gene.present, along=2))
rownames(gene.present.df)=genes

library(plyr)
mat=apply(gene.present.df, 2, function(x) as.numeric(as.character(revalue(as.factor(x), c("yes"="1","no"="0")))))
rownames(mat)=rownames(gene.present.df)

studiesallgenes=names(which((apply(gene.present.df, 2, function(x) length(which(x=="yes")))/length(genes))>=0.8))
length(studiesallgenes)
mat=mat[,studiesallgenes]

#Now get common genes
gene.present.df.allgenes=gene.present.df[,studiesallgenes]
gene.present.df.allgenes=gene.present.df.allgenes[which(apply(gene.present.df.allgenes, 1, function(x) all(x=="yes"))),]
dim(gene.present.df.allgenes) #742 genes in 20 studies 

explist.sig=explist[which(names(explist) %in% colnames(gene.present.df.allgenes))]
ex.commongenes=lapply(explist.sig, function(x) x[match(genes, rownames(x)),])
ex.commongenes.df=abind(ex.commongenes, along=2)
ex.commongenes.df=na.omit(ex.commongenes.df)
ex.commongenes.df=ex.commongenes.df

#pheno=MakePhenographClusters(t(ex.commongenes.df))
dim(ex.commongenes.df)
#742/877  genes in 20 studies (1642 primary HNCs)

library(iCellR)
set.seed(123)
#Rphenograph_out <- iCellR::Rphenograph(ex.commongenes.df, k = 45)
#6 clusters
saveRDS(Rphenograph_out, paste0(Resultsdir,"phenograph.lnm.genes.060923.rds"))

df=data.frame(gene=rownames(ex.commongenes.df), pheno.clust=Rphenograph_out[[2]]$membership) %>% 
  mutate(pheno.clust=glue::glue("L{pheno.clust}"))

metatab$cluster=df[match(metatab$gene, df$gene),"pheno.clust"]

saveRDS(metatab, paste0(Datadir, "zscores_lnm_all_genes_updated.06.09.23.rds"))
metatab=readRDS(paste0(Datadir, "zscores_lnm_all_genes_updated.06.09.23.rds"))


#########################################
#Make umap of LNM gene clusters
#########################################

Rphenograph_out=readRDS(paste0(Resultsdir,"phenograph.lnm.genes.060923.rds"))

library(umap)
set.seed(234)
clust.umap = umap(ex.commongenes.df)
saveRDS(clust.umap, paste0(Resultsdir, "UMAP_PRECOG_HNSCC_LNM_gene_clusters.updated.061623.rds"))

lnm.cols=gnuplot_colors[c(8:11,13:14)]

fm=readRDS(paste0(Datadir, "Findmarkers_Puram_U54_allcombined.update.061323.rds"))
fm2=fm[match(rownames(clust.umap$layout), fm$gene),]

file=paste(Figuresdir, "UMAP_PRECOG_HNSCC_LNM_genes_updated_061723", sep="")
pdf(file=paste0(file,'.pdf',sep=''), height = 5,  width = 5, family = "Helvetica")
plot(clust.umap$layout, ylab="UMAP_2", xlab="UMAP_1",  cex=0.6, pch=as.numeric(fm2$lnm.direction), col=lnm.cols[as.numeric(fm2$lnm.gene.cluster.phenograph)])
legend("right", legend=c("Anti-LNM","Pro-LNM"), pch=c(1,2), title="LNM direction", cex=1, bty="n")
legend("bottomleft", legend=c("L1","L2","L3","L4","L5","L6"), col=lnm.cols, pch=15, title="LNM gene\ncluster", cex=1, bty="n")
dev.off()

################################################
#Boxplot of LNM-associated genes
################################################

library(dplyr)
library(ggrepel)
options(ggrepel.max.overlaps = Inf) 

metatab=readRDS(paste0(Datadir, "zscores_lnm_all_genes_updated.06.09.23.rds"))

metatab=metatab %>% dplyr::rename('lnm.gene.cluster.phenograph' = 'cluster', 'metaz.lnm'='zval')

library("AnnotationDbi")
library("org.Hs.eg.db")

metatab$symbol_AnnotationDbi=as.character(mapIds(org.Hs.eg.db, keys=metatab$gene, column="SYMBOL", keytype="ENTREZID", multiVals="first"))

dat=data.frame(Meta.z=metatab$metaz.lnm, Gene.cluster=factor(metatab$lnm.gene.cluster.phenograph), gene=metatab$symbol_AnnotationDbi)
dat=dat[!is.na(metatab$lnm_direction) & !is.na(dat$gene),]

dat$Gene.cluster2=factor(ifelse(is.na(dat$Gene.cluster),"Not.clustered",as.character(dat$Gene.cluster)))

dat$dir=ifelse(dat$Meta.z<(-3.09),"Anti-LNM","Pro-LNM")
dat$dir=factor(ifelse(dat$Meta.z>(3.09),"Pro-LNM","Anti-LNM"), levels=c("Pro-LNM","Anti-LNM"))

n=10
dat$Meta.z.rank=rep(1, nrow(dat))
dat$Meta.z.rank[order(dat$Meta.z,decreasing=T)[1:n]]=2
dat$Meta.z.rank[order(dat$Meta.z,decreasing=F)[1:n]]=3
dat$Meta.z.rank=factor(dat$Meta.z.rank)

t1 <- ggplot(transform(dat, xjit=jitter(as.numeric(as.factor(Gene.cluster2)))), aes(x=as.factor(Gene.cluster2),y=Meta.z))
t2 <- geom_boxplot(outlier.shape=NA, aes(fill=Gene.cluster2))
t3 <- geom_point(size=0.8, aes(x=xjit, color = Meta.z.rank))
t5 <- geom_text_repel(aes(x=xjit, label=ifelse(Meta.z.rank!=1,as.character(gene),'')), vjust = 0.5, size = 2.5)

lnm.cols=gnuplot_colors[c(8:11,13:14)]
cols2=c(lnm.cols, "lightgrey")

p=t1+t2+t3 + facet_grid(dir ~ ., scale='free_y') +t5
p=p + scale_color_manual(values = c("black", "red", "blue"), labels=c("Not top rank","Highest 10","Lowest 10"))+labs(color="Meta-Z rank", x="Gene cluster", y="Meta-z score")
p=p+ scale_fill_manual(values=cols2, name="LNM gene cluster", labels=c("Not.clustered"="No cluster"))
p

file=paste0(Figuresdir, "jitter_clusters_lnm_updated_061723")
pdf(file=paste0(file,'.pdf',sep=''), height = 6,  width = 8, family = "Helvetica")
p
dev.off()

#################################################################################
#################################################################################
#Meta analysis of genes associated with tumor grade 
######################################################################
#####################################################################

#script: Precog_meta_analysis_genes_associated_with_differentiation.R
glm.apply.provide.genesig=function(genesig, varfactor, genesetname){
  library(survival)
  library(abind)
  
  out <- tryCatch(
    {
      glm=summary(glm(colMeans(exp[rownames(exp) %in% genesig,])~as.numeric(info[,varfactor])))
      
      z = coef(glm)[6]
      se = coef(glm)[4]
      n = length(na.omit(info[,varfactor]))
      p=coef(glm)[8]
      gs=genesetname
      df = data.frame(Z = z, SE = se, N = n, DatasetID = gs, P=p)
    },
    error=function(cond) {
      message("Here's the original error message:")
      message(cond)
      return(NA)
    }, 
    warning=function(cond) {
      message("Here's the original warning message:")
      message(cond)
      return(df)
    }
  )
  return(out)  
}

library(abind)
library(metafor)

make.meta.glm=function(glmlist, glmlistname){
  allglm.abind=abind(glmlist, along=1)
  
  allglm.abind=as.data.frame(abind(glmlist, along=1))
  rownames(allglm.abind)=NULL
  allglm.abind$Z=as.numeric(as.character(allglm.abind$Z))
  allglm.abind$N=as.numeric(as.character(allglm.abind$N))
  
  metaz=sum(allglm.abind$N*allglm.abind$Z)/sqrt(sum(allglm.abind$N^2))
  
  #apply metfor and make forest plot
  d=as.numeric(allglm.abind$Z)
  se=as.numeric(as.character(allglm.abind$SE))
  names1=as.character(allglm.abind$DatasetID)
  g=metafor::rma(d, sei=se, data=allglm.abind)
  
  forestplot=forest(g, slab=allglm.abind$DatasetID, main=glmlistname)
  
  res=list(allglm.abind, metaz, g, forestplot)
  names(res)=c("glm.table","meta.Z","rma_results","forestplot")
  
  return(res)
}

#############################################
#########################################
#Find genes associated with grade
#########################################
############################################

source("~/Documents/Projects/HNSCC_PreCog/data/Pecog_HNSCC_directories.R")
source(paste0(RscriptsPath, "basic_scripts.R"))

#Run meta-analysis of genes associated with level of differentiation in PRECOG

#have formatted differentiation level
alldata=readRDS("~/Documents/Projects/HNSCC_PreCog/data/Precog.HNSCC.all.exp.clin.data.rds")
explist=alldata$explist
infolist=alldata$infolist

#get study grade/level of differentiation
varfactor="COV_HNSCC_study_grade"
lapply(infolist, function(x) summary(as.factor(x[,varfactor])))

accessions.grade=names(which(unlist(lapply(infolist, function(x) !all(is.na(x[,varfactor]))))))
accessions.grade=setdiff(accessions.grade, "GSE686")

infolist.grade=infolist[accessions.grade]
explist.grade=explist[accessions.grade]

#double checking that colnames match between clinical and gene expression data
lapply(1:length(infolist.grade), function(x) all(infolist.grade[[x]][,"SAMPLID"]==colnames(explist.grade[[x]])))

genes=rownames(explist.grade[[1]])[1:10]

#Number of patients with grade data
lapply(explist.grade, ncol) %>% unlist() %>% sum()

#################################
#
#################################

###Below is signature level
#Could run marker genes for each cell Seurat cell type from Puram 

var="COV_HNSCC_study_grade"

glmlists=list()

for(i in 1:length(accessions.grade)){
  
  acc=accessions.grade[i]
  info=infolist.grade[[acc]]
  exp=explist.grade[[acc]]
  exp=as.matrix(exp)
  
  glms=lapply(rownames(exp), function(x) summary(glm(exp[x,]~as.numeric(info[,varfactor]))))
  names(glms)=rownames(exp)
  glm.dfs=lapply(glms, function(glm) data.frame(Z = coef(glm)[6], SE = coef(glm)[4], N = length(na.omit(info[,varfactor])), DatasetID = acc, P=coef(glm)[8]))
  glm.df=as.data.frame(abind::abind(glm.dfs, along=1))
  glm.df$Z=as.numeric(as.character(glm.df$Z))
  glm.df$SE=as.numeric(as.character(glm.df$SE))
  glm.df$N=as.numeric(as.character(glm.df$N))
  glm.df$P=as.numeric(as.character(glm.df$P))
  glm.df$gene=names(glms)
  rownames(glm.df)=names(glms)
  
  glmlists[[acc]]=glm.df
}

saveRDS(glmlists, paste0(Resultsdir, "glms.all.genes.all.studies.histological.grade.updated.061623.rds"))
#13 studies
glmlists=readRDS(paste0(Resultsdir, "glms.all.genes.all.studies.histological.grade.updated.061623.rds"))

sum(unlist(lapply(glmlists, function(x) unique(x$N)))) #N=1075 primary tumors

###################
#Combine p-values using Liptak's test
###############

library(abind)
allglm.abind=as.data.frame(abind(glmlists, along=1))
rownames(allglm.abind)=NULL
allglm.abind$Z=as.numeric(as.character(allglm.abind$Z))
allglm.abind$N=as.numeric(as.character(allglm.abind$N))

#restricting to genes that found within at least two studies, otherwise its not a meta-analysis
tab=table(allglm.abind$gene)
meta.genes=names(tab[which(tab>=2)])
length(meta.genes) #25058 genes for which at least two studies have the gene. Running meta analysis on all of these 

Liptak_combine_z
function(gene, df){
  mdf=df[df$gene==gene,]
  metaz=sum(mdf$N*mdf$Z)/sqrt(sum(mdf$N^2))
  return(metaz)
}
metaz.glm=lapply(meta.genes, function(x)  Liptak_combine_z(gene=x, df=allglm.abind))
names(metaz.glm)=meta.genes
metaz.glm=unlist(metaz.glm)

#Annotate genes using biomaRt
library(biomaRt)
ensembl=useMart("ensembl", dataset="hsapiens_gene_ensembl")

annot=getBM(attributes=c("hgnc_symbol","entrezgene_id"), filters = "entrezgene_id", values = names(metaz.glm), mart= ensembl)

metaz.glm.annot=cbind(metaz.glm, annot[match(meta.genes, annot$entrezgene_id),])
metaz.glm.annot$gene=rownames(metaz.glm.annot)

library("AnnotationDbi")
library("org.Hs.eg.db")

metaz.glm.annot$AnnotationDbi=mapIds(org.Hs.eg.db, keys=metaz.glm.annot$gene, column="SYMBOL", keytype="ENTREZID", multiVals="first")

#add the number of studies for which each gene was available
metaz.glm.annot$n.studies=tab[match(metaz.glm.annot$gene, names(tab))]

metaz.glm.annot=readRDS(paste(Resultsdir, "Precog.HNSCC.glm.grade.Liptak.metaz.allHNSCC.updated.061623.rds", sep=""))

threshold="3.09"
metaz.glm.annot$signficant=ifelse(abs(metaz.glm.annot$metaz.glm)>=threshold,"yes","no")
summary(as.factor(metaz.glm.annot$signficant))
#5045/25058 genes

metaz.glm.annot$direction=rep(NA, nrow(metaz.glm.annot))
metaz.glm.annot$direction[metaz.glm.annot$metaz.glm<=-3.09]="Anti-grade"
metaz.glm.annot$direction[metaz.glm.annot$metaz.glm>=3.09]="Pro-grade"
summary(as.factor(metaz.glm.annot$direction))
1960+3108 #5068 genes

saveRDS(metaz.glm.annot, paste(Resultsdir, "Precog.HNSCC.glm.grade.Liptak.metaz.allHNSCC.updated.061623.rds", sep=""))

metaz.glm.annot.down=metaz.glm.annot[metaz.glm.annot$signficant=="yes" & metaz.glm.annot$metaz.glm<0,]
metaz.glm.annot.down=metaz.glm.annot.down[order(metaz.glm.annot.down$metaz.glm),]
#1937 genes negatively associated with grade

metaz.glm.annot.up=metaz.glm.annot[metaz.glm.annot$signficant=="yes" & metaz.glm.annot$metaz.glm>0,]
metaz.glm.annot.up=metaz.glm.annot.up[order(metaz.glm.annot.up$metaz.glm, decreasing = T),]
#3108 genes positively associated with grade

#####################################################
#Write table of grade-associated genes 
#####################################################

metaz.glm.annot=readRDS(paste(Resultsdir, "Precog.HNSCC.glm.grade.Liptak.metaz.allHNSCC.updated.061623.rds", sep=""))

colnames(metaz.glm.annot)
colnames=c("entrezgene_id", "AnnotationDbi", "metaz.glm","p.value","q.value" ,"direction")

tab=metaz.glm.annot[,colnames]
colnames(tab)=c("entrezgene ID",	"Gene symbol",	"Meta-z score",	"P-value","Q-value (FDR)", "Grade directon of association")

tab=tab[abs(tab$`Meta-z score`)>=3.09,]
dim(tab)
tab=tab[order(abs(tab$`Meta-z score`), decreasing = T),]

tab$`P-value`=signif(tab$`P-value`, digits=3)
tab$`Q-value (FDR)`=signif(tab$`Q-value (FDR)`, digits=3)
tab$`Meta-z score`=round(tab$`Meta-z score`,3)
tab$`entrezgene ID`=as.character(tab$`entrezgene ID`)

#
saveRDS(tab, paste0(Resultsdir, "Genes.grade.all.updated.061623.rds"))
write.table(tab, file=paste0(Resultsdir, "Genes.grade.all.updated.061623.txt"), sep="\t", row.names = FALSE)

####




acc="TCGA"
exp=as.matrix(explist[[acc]])
info=infolist[[acc]]

#Use BIRC5 as an example of an stem-like malignant cell expressed gene that is associated with metastasis and is also positively associated with grade 
symbol="FAM83G"
symbol="CHPT1"
symbol="BIRC5"
symbol="NEMP1"
symbol="TOP2A"

gene=metaz.glm.annot[which(metaz.glm.annot$AnnotationDbi==symbol),"gene"]
gene

file=paste0(Figuresdir, "NEMP1_expression_grade.updated.061623")
pdf(file=paste0(file,'.pdf',sep=''), height = 6,  width = 6, family = "Helvetica")
boxplot(exp[gene,]~info[,"COV_HNSCC_study_grade"], main=symbol, ylab="Gene expression", xlab="Histological grade")
dev.off()

#
head(metaz.glm.annot.up, 50)

################################################################
#Make heatmap of z-scores for genes associated with grade 
##############################################################

glmlists=readRDS(paste0(Resultsdir, "glms.all.genes.all.studies.histological.grade.updated.061623.rds"))
library(abind)
allglm.abind=as.data.frame(abind(glmlists, along=1))
rownames(allglm.abind)=NULL
allglm.abind$Z=as.numeric(as.character(allglm.abind$Z))
allglm.abind$N=as.numeric(as.character(allglm.abind$N))

MetaFrame=allglm.abind

genes=as.character(metaz.glm.annot[metaz.glm.annot$signficant=="yes","gene"])
studies=levels(as.factor(MetaFrame$DatasetID))

genes=genes[order(metaz.glm.annot[genes,"metaz.glm"])]

#zarray=array(NA, c(length(genes), length(studies)))
#rownames(zarray)=genes
#colnames(zarray)=studies

zarray=as.data.frame(abind(lapply(studies, function(x) glmlists[[x]][match(genes, rownames(glmlists[[x]])),"Z"]), along=2))
colnames(zarray)=studies
rownames(zarray)=genes

saveRDS(zarray, paste0(Resultsdir,"z.score.matrix.meta.analysis.grade.updated.061623.rds"))
zarray=readRDS(paste0(Resultsdir,"z.score.matrix.meta.analysis.grade.updated.061623.rds"))

library(ComplexHeatmap)

df=data.frame(Direction=ifelse(metaz.glm.annot[genes,"metaz.glm"]<0,"Pro-grade","Anti-grade"))

cols=c("magenta","green")
names(cols)=c("Pro-grade", "Anti-grade")
rowcols=list(Direction=cols)

rowannot=ComplexHeatmap::rowAnnotation(df = df, col = rowcols, show_annotation_name = TRUE)

studymeta=readRDS(paste0(Datadir, "table_1_adding_authors_OPH_HPV_info.rds"))
studymeta=studymeta[match(colnames(zarray), studymeta$Study.accession),]
colnames(zarray)=studymeta$Study_label_author

zarray=as.matrix(zarray)

p=ComplexHeatmap::Heatmap(zarray, cluster_rows = FALSE, name = "Z-score", cluster_columns = FALSE, show_row_names =F, show_column_names =T, right_annotation=rowannot, row_title = "Pathological grade-associated genes", column_title_side = "bottom", column_title = "Primary HNC patient studies", raster_device = 'png', use_raster=FALSE)

file=paste0(Figuresdir, "Heatmap_z_scores_differentiation.updated.061623")
pdf(file=paste0(file,'.pdf',sep=''), height = 6,  width = 6, family = "Helvetica")
print(p)
dev.off()

#######################
#Heatmap z-scores grade
########################

metazs=metaz.glm.annot[match(rownames(zarray), metaz.glm.annot$gene),"metaz.glm"]

col_fun3 = circlize::colorRamp2(c(min(metazs), 0, max(metazs)), c("green","black","red"))
rowannot.list.lnm=list(`Meta-z-score`=col_fun3)

geneannot2=rowAnnotation(`Meta-z-score`=metazs, col=rowannot.list.lnm, annotation_legend_param = list(`Meta-z-score` = list(direction = "horizontal")))

p=ComplexHeatmap::Heatmap(zarray, cluster_rows = FALSE, name = "Z-score", cluster_columns = FALSE, show_row_names =F, show_column_names =T, left_annotation=rowannot,right_annotation=geneannot2, heatmap_legend_param = list(legend_direction="horizontal"), row_title = "Pathological grade-associated genes", column_title_side = "bottom", column_title = "Primary HNC patient studies", raster_device = 'png', use_raster=FALSE)

file=paste0(Figuresdir, "Heatmap_z_scores_differentiation_with_metaz.updated.061623")
pdf(file=paste0(file,'.pdf',sep=''), height = 6,  width = 6, family = "Helvetica")
print(p)
dev.off()

geneannot2=rowAnnotation(`Meta-z-score`=metazs, col=rowannot.list.lnm, show_legend=FALSE)
#geneannot2=rowAnnotation(`Meta-z-score`=metazs, foo = anno_mark(at = match(genes, rownames(zarray)), labels = genes), col=rowannot.list.lnm, show_legend=FALSE)

zs=c(-10,-5,0,5,10)
ps=2*pnorm(-abs(zs))
ps=ps=signif(ps, digits=3)
ps=format(ps, scientific=TRUE)
labs=paste0(zs," [",ps,"]")
val=zs
col_fun2 = circlize::colorRamp2(c(min(val), 0, max(val)), c("green","black","red"))
lgd = Legend(col_fun = col_fun2, title = "Meta-z score [Corresponding p-value]",labels = labs, at=zs)

zs2=c(-5,0,5)
ps2=2*pnorm(-abs(zs2))
ps2=signif(ps2, digits=3)
ps2=format(ps2, scientific=TRUE)
labs2=paste0(zs2," [",ps2,"]")

val2=zs2
col_fun3 = circlize::colorRamp2(c(min(val2), 0, max(val2)), c("blue","white","red"))
lgd2 = Legend(col_fun = col_fun3, title = "Linear regression z-score [Corresponding p-value]",labels = labs2, at=zs2)
lgd_list=list(lgd2, lgd)

p=ComplexHeatmap::Heatmap(zarray, show_heatmap_legend = FALSE, cluster_rows = FALSE, name = "z-score", cluster_columns = FALSE, show_row_names =F, show_column_names =T, left_annotation=rowannot,right_annotation=geneannot2, heatmap_legend_param = list(legend_direction="horizontal"), row_title = "Pathological grade-associated genes", column_title_side = "bottom", column_title = "Primary HNC patient studies", raster_device = 'png', use_raster=FALSE)

file=paste0(Figuresdir, "Heatmap_z_scores_differentiation_with_metaz_with_pvals.updated.061623")
pdf(file=paste0(file,'.pdf',sep=''), height = 6,  width = 10, family = "Helvetica")
draw(p,  annotation_legend_list = lgd_list, heatmap_legend_side='left', annotation_legend_side='right')
dev.off()

#saving data
datalist=list(zarray=zarray, metazs=metazs)
saveRDS(datalist, paste0(Figuresdir.grobs, "data.Heatmap_z_scores_differentiation_with_metaz_with_pvals.updated.061623.rds"))


####################################################
####################################################
#Venn diagram overlap of LNM-associated genes with grade genes 
####################################################
####################################################

metaz.glm.annot=readRDS(paste(Resultsdir, "Precog.HNSCC.glm.grade.Liptak.metaz.allHNSCC.updated.061623.rds", sep=""))

head(metaz.glm.annot)

fm=readRDS(paste0(Datadir, "Findmarkers_Puram_U54_allcombined.update.061323.rds"))

pro.lnm=fm[which(fm$lnm.direction=="Pro-LNM"),"gene"]
anti.lnm=fm[which(fm$lnm.direction=="Anti-LNM"),"gene"]

anti.grade=metaz.glm.annot[metaz.glm.annot$metaz.glm<(-3.09),"gene"]
pro.grade=metaz.glm.annot[metaz.glm.annot$metaz.glm>3.09,"gene"]

length(intersect(pro.lnm, pro.grade))#324/3108 pro-LNM genes
length(intersect(pro.lnm, anti.grade))#0 genes

length(intersect(anti.lnm, anti.grade))#374/1960 LNM genes
intersect(anti.lnm, pro.grade)#0 genes

#install.packages("VennDiagram")
library(VennDiagram)

plot1=draw.pairwise.venn(length(pro.lnm), length(pro.grade), length(intersect(pro.lnm, pro.grade)), category = c("Pro-LNM", "Pro-grade"), lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2))

file=paste(Figuresdir, "Venn_overlap_proLNM_proGrade.updated.061623")
pdf(file=paste(file,'.pdf',sep=''), height = 3.5,  width = 3.5, family = "Helvetica")
par(oma=c(1,1,1,5))
draw.pairwise.venn(length(pro.lnm), length(pro.grade), length(intersect(pro.lnm, pro.grade)), category = c("Pro-LNM", "Pro-grade"), lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2))
dev.off()

plot2=draw.pairwise.venn(length(anti.lnm), length(anti.grade), length(intersect(anti.lnm, anti.grade)), category = c("Anti-LNM", "Anti-grade"), lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2))

file=paste(Figuresdir, "Venn_overlap_antiLNM_antiGrade.updated.061623")
pdf(file=paste(file,'.pdf',sep=''), height = 3.5,  width = 3.5, family = "Helvetica")
par(oma=c(1,1,1,5))
draw.pairwise.venn(length(anti.lnm), length(anti.grade), length(intersect(anti.lnm, anti.grade)), category = c("Anti-LNM", "Anti-grade"), lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2))
dev.off()

#
threshold="3.09"
metaz.glm.annot$signficant=ifelse(abs(metaz.glm.annot$metaz.glm)>=threshold,"yes","no")
metaz.glm.annot$p.value=2*pnorm(-abs(metaz.glm.annot$metaz.glm))
metaz.glm.annot$q.value=p.adjust(metaz.glm.annot$p.value, method="fdr")

saveRDS(metaz.glm.annot, paste(Resultsdir, "Precog.HNSCC.glm.grade.Liptak.metaz.allHNSCC.updated.061623.rds", sep=""))
metaz.glm.annot=readRDS(paste0(Resultsdir, "Precog.HNSCC.glm.grade.Liptak.metaz.allHNSCC.updated.061623.rds"))

#################################
#Plot meta-z scores lnm versus grade
#################################

lnm.list=readRDS(paste0(Resultsdir, "Precog.meta.analysis.node.status.all.HNSCC.updated.06.09.23.rds"))
meta2=lnm.list$meta

metatab=cbind(
  unlist(lapply(meta2, function(x) x$pval)),
  unlist(lapply(meta2, function(x) x$b)),
  unlist(lapply(meta2, function(x) x$zval)))

metatab=as.data.frame(metatab)
colnames(metatab)=c("pval","beta","zval")
metatab$gene=rownames(metatab)

metaz.glm.annot$rank_positively_associated_grade=rank(metaz.glm.annot$metaz.glm*-1)
metaz.glm.annot$rank_negatively_associated_grade=rank(metaz.glm.annot$metaz.glm)

metatab[,c("metaz.glm","rank_positively_associated_grade","rank_negatively_associated_grade","AnnotationDbi")]=metaz.glm.annot[match(metatab$gene, metaz.glm.annot$gene),c("metaz.glm","rank_positively_associated_grade","rank_negatively_associated_grade","AnnotationDbi")]

metatab$rank_positively_associated_lnm=rank(metatab$zval*-1)
metatab$rank_negatively_associated_lnm=rank(metatab$zval)

###Make scatter plot of meta-z scores for LNM and grade
library(ggrepel)

remotes::install_github("JosephCrispell/basicPlotteR")

coords=subset(metatab, rank_positively_associated_grade<=50 & rank_positively_associated_lnm<=50|rank_negatively_associated_grade<=50 & rank_negatively_associated_lnm<=50)

cc=round(cor(metatab$zval, metatab$metaz.glm, use="pairwise.complete.obs"),2)

file=paste(Figuresdir, "plot_metazs_lnm_versus_grade_smoothscatter_with_labels.updated.060923", sep="")
pdf(file=paste0(file,'.pdf',sep=''), height = 7,  width = 7, family = "Helvetica")
with(metatab, smoothScatter(metaz.glm, zval, xlab="Meta-z score tumor grade", ylab="Meta-z score LNM", xlim=c(-12.5,10.5)))
#with(subset(metatab, abs(zval)>=3.09), points(metaz.glm, zval, pch=20, col="red"))
abline(v=-3.09, lty=2, col="magenta")
abline(v=3.09, lty=2, col="magenta")
abline(h=-3.09, lty=2, col="green")
abline(h=3.09, lty=2, col="green")
legend("bottomright", legend=paste0("r=", cc), col="red", bty="n")
abline(lm(metatab$zval~metatab$metaz.glm), lty=2, col="red")
basicPlotteR::addTextLabels(coords$metaz.glm, coords$zval, coords$AnnotationDbi, cex.label=1, col.label="black", col.line="grey")
dev.off()


####################################################
###############################################
#Heatmap of association of signatures with grade
###############################################
####################################################

#Get forest plot figures to test for linear association with grade 

#have formatted differentiation level
alldata=readRDS("~/Documents/Projects/HNSCC_PreCog/data/Precog.HNSCC.all.exp.clin.data.curated.for.PRECOG.rds")

explist=alldata$explist
infolist=alldata$infolist

#restrict to datasets with >1000 genes 
gene.n.threshold=1000
accessions=names(which(lapply(explist, nrow)>=gene.n.threshold))
explist=explist[accessions]

#Get gene signatures 
fm=readRDS(paste0(Datadir, "Findmarkers_Puram_U54_allcombined.update.061323.rds"))

makesig=function(genetype, group){
  levs=levels(fm[,group])
  sigs=lapply(levs, function(lev) 
  {x=fm[fm[,group]==lev,genetype]
  x=x[!is.na(x)]
  return(x)})
  names(sigs)=levs
  return(sigs)
}

#sigs=c(makesig(genetype="gene",group="lnm.direction"),
#       makesig(genetype="gene",group="lnm.gene.cluster.phenograph"))

sigs=c(makesig(genetype="gene",group="lnm.direction"),
       makesig(genetype="gene",group="lnm.gene.cluster.phenograph"),
       makesig(genetype="gene",group="survival.direction"),
       makesig(genetype="gene",group="survival.gene.cluster.phenograph"))

siglist=list()
#
for(i in 1:length(explist)){
  exp=explist[[i]]
  sigdf=as.data.frame(t(abind::abind(lapply(sigs, function(x) scale(colMeans(exp[rownames(exp) %in% x,], na.rm = T))), along=2)))
  siglist[[i]]=sigdf
}
names(siglist)=names(explist)
explist=siglist

#get study grade/level of differentiation
varfactor="COV_HNSCC_study_grade"
lapply(infolist, function(x) summary(as.factor(x[,varfactor])))

accessions.grade=names(which(unlist(lapply(infolist, function(x) !all(is.na(x[,varfactor]))))))
accessions.grade=intersect(accessions.grade, names(explist))

infolist.grade=infolist[accessions.grade]
explist.grade=siglist[accessions.grade]

#double checking that colnames match between clinical and gene expression data
lapply(1:length(infolist.grade), function(x) all(infolist.grade[[x]][,"SAMPLID"]==colnames(explist.grade[[x]])))

###
glm.apply.provide.genesig=function(signame, varfactor, genesetname){
  library(survival)
  library(abind)
  
  out <- tryCatch(
    {
      glm=summary(glm(exp[signame,]~as.numeric(info[,varfactor])))
      
      z = coef(glm)[6]
      se = coef(glm)[4]
      n = length(na.omit(info[,varfactor]))
      p=coef(glm)[8]
      gs=genesetname
      df = data.frame(Z = z, SE = se, N = n, DatasetID = gs, P=p)
    },
    error=function(cond) {
      message("Here's the original error message:")
      message(cond)
      return(NA)
    }, 
    warning=function(cond) {
      message("Here's the original warning message:")
      message(cond)
      return(df)
    }
  )
  return(out)  
}


library(abind)
library(metafor)

make.meta.glm=function(glmlist, glmlistname){
  allglm.abind=abind(glmlist, along=1)
  
  allglm.abind=as.data.frame(abind(glmlist, along=1))
  rownames(allglm.abind)=NULL
  allglm.abind$Z=as.numeric(as.character(allglm.abind$Z))
  allglm.abind$N=as.numeric(as.character(allglm.abind$N))
  
  metaz=sum(allglm.abind$N*allglm.abind$Z)/sqrt(sum(allglm.abind$N^2))
  
  #apply metfor and make forest plot
  d=as.numeric(allglm.abind$Z)
  se=as.numeric(as.character(allglm.abind$SE))
  names1=as.character(allglm.abind$DatasetID)
  g=metafor::rma(d, sei=se, data=allglm.abind)
  
  forestplot=forest(g, slab=allglm.abind$DatasetID, main=glmlistname)
  
  res=list(allglm.abind, metaz, g, forestplot)
  names(res)=c("glm.table","meta.Z","rma_results","forestplot")
  
  return(res)
}

#Run meta-analysis
sigs=rownames(explist.grade[[1]])

var="COV_HNSCC_study_grade"

glmlists=list()

for(j in 1:length(sigs)){
  
  sig=sigs[[j]]
  glmlist=list()
  
  for(i in 1:length(accessions.grade)){
    acc=accessions.grade[i]
    info=infolist.grade[[acc]]
    exp=as.matrix(explist.grade[[acc]])
    
    glmlist[[acc]]=glm.apply.provide.genesig(sig, var, acc)
  }
  
  glmlists[[j]]=glmlist
}
names(glmlists)=names(sigs)

#
meta.analysis.signatures.grade=lapply(1:length(glmlists), function(x) make.meta.glm(glmlists[[x]], names(glmlists)[x]))
names(meta.analysis.signatures.grade)=names(glmlists)

g=meta.analysis.signatures.grade$node_status_down$rma_results
allglm.abind=meta.analysis.signatures.grade$node_status_down$glm.table

metasigs=meta.analysis.signatures.grade

pvals=unlist(lapply(metasigs, function(x) x$rma_results$pval))
pvals2=pvalr(unlist(lapply(metasigs, function(x) x$rma_results$pval)))
metazs=unlist(lapply(1:length(metasigs), function(x) round(metasigs[[x]]$meta.Z,2)))
names(metazs)=rownames(explist.grade$TCGA)
#make zarray to make heatmap

#
genes=names(pvals2)
studies=accessions.grade

zarray=abind::abind(lapply(glmlists, function(x) as.numeric(as.character(as.data.frame(abind(x, along=1))[,"Z"]))), along=2)
rownames(zarray)=studies
colnames(zarray)=sigs
#names(glmlists$Survival_all_positive)==studies
zarray=t(zarray)[1:8,]

results=list(meta=meta.analysis.signatures.grade, zarray=zarray, glmlists=glmlists)
saveRDS(results, paste0(Resultsdir, "meta.analysis.signatures.grade.updated.061723.rds"))
results=readRDS(paste0(Resultsdir, "meta.analysis.signatures.grade.updated.061723.rds"))

#Make version that is just LNM signatures 
#Need to formally test the association of signatures with HPV if I'm doing it for tumor size
#rownames(zarray)=gsub("_"," ",rownames(zarray))
rownames(zarray)=paste0(rownames(zarray)," | Meta-Z score: ", metazs[1:8])

library(ComplexHeatmap)

#Make version for just LNM
p=ComplexHeatmap::Heatmap(zarray[1:8,], cluster_rows = F, cluster_columns = F,  name="Z-score", row_title = "Prognostic signatures", column_title = "Primary HNC patient studies", column_title_side = "bottom", show_row_names = T,  heatmap_legend_param = list(legend_direction="vertical"))

file=paste(Figuresdir, "Heatmap_PRECOG_signatures_grade_just_LNM_updated_061723")
pdf(file=paste(file,'.pdf',sep=''), height = 3.5,  width = 8, family = "Helvetica")
ComplexHeatmap::draw(p, heatmap_legend_side="left", padding = unit(c(2, 2, 2, 20), "mm"))
dev.off()

##############################################
#Heatmap dot plot of correlations between EMT and progostic sigs 
###############################################

signature="EMT.score"

int=readRDS("~/Documents/Projects/U54/Puram/data/Integated_Puram_prim_patients200cells_update_160623.rds")
DefaultAssay(int)="RNA"

keep=colnames(int[,int$cell.type.collapsed!="Unclassified" & int$cell.type.collapsed!="myocyte"])
int=subset(int, cells=keep)

clustnames=c("Anti-survival","Pro-survival", paste0("S", c(1:6)),"Anti-LNM","Pro-LNM", paste0("L", c(1:6)))

stats=lapply(1:length(clustnames), function(i) {
  surv.clust=clustnames[i]
  clustname=clustnames[i]
  cc=round(cor(scale(int@meta.data[,surv.clust][int$cell.type.collapsed=="Malignant"]), scale(int@meta.data[,signature][int$cell.type.collapsed=="Malignant"]), use="pairwise.complete.obs"),2)
  p=cor.test(scale(int@meta.data[,surv.clust][int$cell.type.collapsed=="Malignant"]), scale(int@meta.data[,signature][int$cell.type.collapsed=="Malignant"]), use="pairwise.complete.obs")$p.value
  stats=data.frame(sig=clustname, cc=cc, p=p)
  return(stats)
})

dat=as.data.frame(abind::abind(stats, along=1))
dat.Puram=dat
dat.Puram$Study=rep("Puram", nrow(dat.Puram))

int=readRDS(paste0(DatadirCCSBSc, "CCSB_scRNASeq_HNSCC_primary_enzymatic_intgrated_50PCs.mito.removed.rds"))
DefaultAssay(int)="RNA"

int$cell.type.collapsed=factor(int$cell.type.collapsed, levels=c("B cell","Endothelial","Fibroblast","Malignant","Mast", "Myeloid","T cell or NK cell"))

keep=colnames(int[,int$celltype!="Unclassified"])
int=subset(int, cells=keep)

stats=lapply(1:length(clustnames), function(i) {
  surv.clust=clustnames[i]
  clustname=clustnames[i]
  cc=round(cor(scale(int@meta.data[,surv.clust][int$cell.type.collapsed=="Malignant"]), scale(int@meta.data[,signature][int$cell.type.collapsed=="Malignant"]), use="pairwise.complete.obs"),2)
  p=cor.test(scale(int@meta.data[,surv.clust][int$cell.type.collapsed=="Malignant"]), scale(int@meta.data[,signature][int$cell.type.collapsed=="Malignant"]), use="pairwise.complete.obs")$p.value
  stats=data.frame(sig=clustname, cc=cc, p=p)
  return(stats)
})

dat=as.data.frame(abind::abind(stats, along=1))
dat.Stanford=dat
dat.Stanford$Study=rep("Stanford", nrow(dat.Stanford))

dat=rbind(dat.Puram, dat.Stanford)

#
dat$cc=as.numeric(dat$cc)
dat$p=as.numeric(dat$p)
dat$Minlog10p=-log10(dat$p)
dat$sig=factor(dat$sig, levels=clustnames)
dat$Study=factor(dat$Study)

require(ggplot2)

dat$fill2=rep(0, nrow(dat))

dev.off()

p=ggplot(dat, aes(y = sig, x = Study)) +  geom_tile(aes(fill=fill2), show.legend = FALSE)
p=p+ geom_point(aes(colour = cc, size =Minlog10p))
p=p +scale_color_gradient2(name = "Pearson r",low = "green", mid="black",  high = "magenta")+scale_fill_gradient2(low = "lightgrey", high = "lightgrey")
p=p+ scale_y_discrete(limits = rev(unique(sort(dat$sig))))
p=p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white", colour = "black",size = 0.5, linetype = "solid"))
p=p+ labs(y="Prognostc signature")
p=p+scale_size_continuous(name="-log10 p-value")

file=paste(FiguresdirPuram, "Puram_heatmap_EMT_versus_prognostic_signatures_updated_180623", sep="")
pdf(file=paste0(file,'.pdf',sep=''), height = 5,  width = 4, family = "Helvetica")
p
dev.off()


##########################################################
#Plotting grade in TCGA
##########################################################

dat2=readRDS("~/Documents/Projects/HNSCC_PreCog/data/TCGA_HNSC_Gencode32_tximport_counts_Combat_entrezid.rds")

alldata=readRDS("~/Documents/Projects/HNSCC_PreCog/data/Precog.HNSCC.all.exp.clin.data.rds")

explist=alldata$explist
infolist=alldata$infolist

#exp=explist[[which(names(explist) %in% "TCGA")]]
info=infolist[[which(names(infolist) %in% "TCGA")]]
info$barcode=substr(gsub("[.]","-",info$EXPR),1,15)

exp=dat2
exp=exp[,substr(colnames(exp),14,15)=="01"|substr(colnames(exp),14,15)=="11"]
colnames(exp)=substr(colnames(exp),1,15)

ann=data.frame(barcode=colnames(exp))
ann$COV_HNSCC_study_grade=as.character(info[match(ann$barcode, info$barcode),"COV_HNSCC_study_grade"])
ann[substr(ann$barcode,14,15)=="11","COV_HNSCC_study_grade"]="Normal"

ann$COV_HNSCC_study_grade=factor(ann$COV_HNSCC_study_grade, levels=c("Normal","G1","G2","G3","G4"))
rownames(ann)=ann$barcode
all(rownames(ann)==colnames(exp))

fm=readRDS(paste0(Datadir, "Findmarkers_Puram_U54_allcombined.update.061323.rds"))

siglist=list(
  Adverse=fm[which(fm$survival.direction=="Anti-survival"),"gene"],
  Favorable=fm[which(fm$survival.direction=="Pro-survival"),"gene"],
  L1=fm[which(fm$lnm.gene.cluster.phenograph=="L1"),"gene"],
  L2=fm[which(fm$lnm.gene.cluster.phenograph=="L2"),"gene"],
  L3=fm[which(fm$lnm.gene.cluster.phenograph=="L3"),"gene"],
  L4=fm[which(fm$lnm.gene.cluster.phenograph=="L4"),"gene"],
  L5=fm[which(fm$lnm.gene.cluster.phenograph=="L5"),"gene"],
  L6=fm[which(fm$lnm.gene.cluster.phenograph=="L6"),"gene"],
  "Pro-LNM"=fm[which(fm$lnm.direction=="Pro-LNM"),"gene"],
  "Anti-LNM"=fm[which(fm$lnm.direction=="Anti-LNM"),"gene"]
)

siglist=lapply(siglist, function(x) x[x %in% rownames(exp)])

sig=siglist$`Anti-LNM`
tab.L1=data.frame(Signature=scale(colMeans(exp[sig,])), Grade=ann$COV_HNSCC_study_grade)
tab.L1$Direction=rep("Anti-LNM", nrow(tab.L1))

sig=siglist$`Pro-LNM`
tab.L4=data.frame(Signature=scale(colMeans(exp[sig,])), Grade=ann$COV_HNSCC_study_grade)
tab.L4$Direction=rep("Pro-LNM", nrow(tab.L4))

tab=rbind(tab.L1, tab.L4)

tab$Direction=factor(tab$Direction)
tab=na.omit(tab)

saveRDS(tab, paste0(Figuresdir.grobs, "data.boxplots.Anti_LNM.Pro_LNM.grade_TCGA_htSEQ_counts_ggplot2_updated_062023.rds"))
tab=readRDS(paste0(Figuresdir.grobs, "data.boxplots.Anti_LNM.Pro_LNM.grade_TCGA_htSEQ_counts_ggplot2_updated_062023.rds"))

cols=RColorBrewer::brewer.pal(5,"YlOrRd")

p=ggplot(tab, aes(x=Grade, y=Signature)) + geom_boxplot(outlier.shape=NA, aes(fill=Grade), lwd=0.2)+ geom_jitter(shape=16, size=0.2, position=position_jitter(0.2)) + facet_grid(~Direction)+ labs(x = "Sample type [Grade]", y="Scaled mean gene expression")
p=p+ scale_fill_manual(name = "Sample type [Grade]", values =cols, labels = c("Normal" = "Tumor-adjacent normal", "G1" = "Malignant [G1]", "G2" = "Malignant [G2]", "G3" = "Malignant [G3]","G4" = "Malignant [G4]"))+ theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
p=p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white", colour = "black",size = 0.5, linetype = "solid"))
p

file=paste0(Figuresdir, "boxplots.Anti_LNM.Pro_LNM.grade_TCGA_htSEQ_counts_ggplot2_updated_062023")
pdf(file=paste0(file,'.pdf',sep=''), height = 2.8,  width = 5, family = "Helvetica")
p
dev.off()

#######################################################################################################
#Meta-analysis of survival-associated genes adjusted for age and sex
#######################################################################################################

alldata2=readRDS("~/Documents/Projects/HNSCC_PreCog/data/Precog.HNSCC.all.exp.clin.data.curated.for.PRECOG.rds")

explist=alldata2$explist
infolist=alldata2$infolist

Survival_objects=readRDS(paste(Resultsdir, "Precog.HNSCC.survival.objects.all.HNSCC.rds", sep=""))
infolist=infolist[names(Survival_objects)]

accessions=names(infolist)

for(acc in accessions){
  infolist[[acc]]=infolist[[acc]][names(Survival_objects[[acc]]),] 
}
#restrict to studies in which age is numeric, to remove the one in which it's categorical

infolist=infolist[lapply(infolist, function(x) is.numeric(x$COVAR_age)) %>% unlist %>% which() %>% names()]

tab.vars=lapply(infolist, function(x) table(x$COVAR_sex)) %>% abind::abind(along=2) %>% t %>% as.data.frame()
tab.vars$n.age=lapply(infolist, function(x) length(which(!is.na(x)))) %>% unlist

accessions=tab.vars %>% filter(female>5 & male>5) %>% rownames()
accessions=setdiff(accessions, "GSE686")
#"GSE23558"    "GSE3292"     "GSE41116"    "GSE65858"      "GSE75538"    "GSE85195"    "GSE95805"    "TCGA"        "E-MTAB-1328" "GSE39366"    "Thurlow"

infolist=infolist[accessions]
explist=explist[accessions]

Survival_objects=readRDS(paste(Resultsdir, "Precog.HNSCC.survival.objects.all.HNSCC.rds", sep=""))
Survival_objects=Survival_objects[accessions]

Survival_objects.surv=list()
infolist.surv=list()
explist.surv=list()

#Go through and restict to studies with at least 20 samples with survival data 
for(acc in accessions){
  Survival_object=Survival_objects[[acc]]
  info=infolist[[acc]]
  exp=explist[[acc]]
  OverlapSamples=intersect(rownames(info), names(Survival_object))
  
  Survival_objects.surv[[acc]]=Survival_object[OverlapSamples]
  infolist.surv[[acc]]=info[OverlapSamples,]
  explist.surv[[acc]]=exp[,OverlapSamples]
}

keep=names(which(unlist(lapply(Survival_objects.surv, function(x) length(x)>=20))))
keep=setdiff(keep, "GSE686") #Manually excluding study with only 1830 genes as it always causes trouble and is the Chung repeat sample dataset
#"GSE23558"    "GSE3292"     "GSE41116"    "GSE65858"    "GSE85195"    "GSE95805"    "TCGA" "E-MTAB-1328" "GSE39366"    "Thurlow"    
#10 studies 

#Restrict to studies with at least 20 samples with survival data after removing HPV positive 
Survival_objects=Survival_objects.surv[keep]
explist=explist.surv[keep]
infolist=infolist.surv[keep]

lapply(Survival_objects, length) %>% unlist %>% sum
#1224

#Now run the meta-analysis in the way that you ran the original one
accessions.survival=names(Survival_objects)

#Analysis adjusted for HPV status
allcoxph=list()
for(i in 1:length(accessions.survival)){
  
  acc=accessions.survival[i]
  Survival_object=Survival_objects[[acc]]
  exp=as.matrix(explist[[acc]][,names(Survival_object)])
  #exp=exp[rownames(exp) %in% fm$gene,]
  info=infolist[[acc]][names(Survival_object),]
  
  survs=apply(exp,1, coxph.apply.provide.survobject.adj.age.sex)
  survs=survs[which(lapply(survs, length)==11)]
  survs=as.data.frame(abind(survs, along=1))
  survs=survs[!is.na(survs$z_gene),]
  
  survs=survs %>% dplyr::mutate_at(c('z_gene', 'z_age', 'N', 'z_sex', 'p_gene', 'p_age', 'p_sex', 'se_gene', 'se_age', 'se_sex'), as.numeric)
  survs$q_gene=p.adjust(survs$p_gene, method="fdr")
  survs$q_age=p.adjust(survs$p_age, method="fdr")
  survs$q_sex=p.adjust(survs$p_sex, method="fdr")
  
  survs$gene=rownames(survs)
  
  allcoxph[[i]]=survs
}
names(allcoxph)=accessions.survival

saveRDS(allcoxph, paste(Resultsdir, "Precog.HNSCC.coxph.survival.adjusted.age.sex.rds", sep=""))
allcoxph=readRDS(paste0(Resultsdir, "Precog.HNSCC.coxph.survival.adjusted.age.sex.rds"))

#Perform meta-analysis of prognostic scores using the Liptak method 
library(abind)
allcoxph.abind=as.data.frame(abind(allcoxph, along=1))
rownames(allcoxph.abind)=NULL

allcoxph.abind=allcoxph.abind%>% dplyr::mutate_at(c('z_gene', 'z_age', 'N', 'z_sex', 'p_gene', 'p_age', 'p_sex', 'se_gene', 'se_age', 'se_sex'), as.numeric)

allglms.abind %>% dplyr::select("DatasetID", "N") %>% 
  filter(DatasetID %in% names(allcoxph)) %>% unique %>% pull(N) %>% as.numeric() %>% sum
#1096 patients used in survival models

#restricting to genes that found within at least two studies, otherwise its not a meta-analysis
tab=table(allcoxph.abind$gene)
meta.genes=names(tab[which(tab>=2)])
#23521 genes for which at least two studies have the gene. Running meta analysis on all of these 

meta.genes=meta.genes %>% 
  as.data.frame() %>% 
  purrr::set_names("gene") 

meta.genes$metaz.gene=lapply(meta.genes$gene, function(gene){
  mdf=allcoxph.abind[allcoxph.abind$gene==gene,]
  metaz=sum(mdf$N*mdf$z_gene)/sqrt(sum(mdf$N^2))
  return(metaz)
}) %>% unlist

meta.genes$metaz.age=lapply(meta.genes$gene, function(gene){
  mdf=allcoxph.abind[allcoxph.abind$gene==gene,]
  metaz=sum(mdf$N*mdf$z_age)/sqrt(sum(mdf$N^2))
  return(metaz)
}) %>% unlist

meta.genes$metaz.sex=lapply(meta.genes$gene, function(gene){
  mdf=allcoxph.abind[allcoxph.abind$gene==gene,]
  metaz=sum(mdf$N*mdf$z_sex)/sqrt(sum(mdf$N^2))
  return(metaz)
}) %>% unlist

colnames(meta.genes)=c("gene","metaz_survival_gene_adj_age_sex", "metaz_survival_age_adj_age_sex", "metaz_survival_sex_adj_age_sex")

saveRDS(meta.genes, paste0(Resultsdir, "Precog.HNSCC.coxph.survival.adjusted.age.sex.meta.stats.rds"))
meta.genes=readRDS(paste0(Resultsdir, "Precog.HNSCC.coxph.survival.adjusted.age.sex.meta.stats.rds"))

fm=readRDS(paste0(Datadir, "Findmarkers_Puram_U54_allcombined.update.061323.rds"))

fm.surv.age.sex=fm %>% dplyr::left_join(meta.genes, by="gene")
fm=fm %>% dplyr::left_join(meta.genes, by="gene")
#saveRDS(fm, paste0(Datadir, "Findmarkers_Puram_U54_allcombined.update.100521.rds"))
#saveRDS(fm, paste0(Datadir, "Findmarkers_Puram_U54_allcombined.update.061323.rds"))

#Investigate further with respect to clusters. 
#Try same with LNM-does the prob persist
#Plot Z scores for all samples against z-scores for HPV negative analysis, i.e. samples that don't include HPV+ve OPH
PRECOG.HNSCC.All=readRDS(paste0(Resultsdir, "Precog.HNSCC.coxph.Liptak.metaz.allHNSCC.updated.06_09_2023.rds"))

PRECOG.HNSCC.All$metaz_survival_gene_adj_age_sex=meta.genes[match(PRECOG.HNSCC.All$gene, meta.genes$gene),"metaz_survival_gene_adj_age_sex"]

plot(PRECOG.HNSCC.All$metaz_survival_gene_adj_age_sex, PRECOG.HNSCC.All$metaz.coxph)
cor.test(PRECOG.HNSCC.All$metaz_survival_gene_adj_age_sex, PRECOG.HNSCC.All$metaz.coxph)
#cor=0.9306237 , p-value < 2.2e-16

file=paste0(Figuresdir, "Scatter_metazs_survival_vs_survival_adj_age_sex")
pdf(file=paste0(file,'.pdf',sep=''), height = 4.5,  width = 4.5, family = "Helvetica")
smoothScatter(PRECOG.HNSCC.All$metaz_survival_gene_adj_age_sex, PRECOG.HNSCC.All$metaz.coxph, xlab="Meta-z score (Adjusted for age & sex)", ylab="Meta-z score")
abline(lm(PRECOG.HNSCC.All$metaz.coxph~PRECOG.HNSCC.All$metaz_survival_gene_adj_age_sex), col="red", lty=2)
dev.off()

#Smooth scatter and correlations just for survival-associated genes 
file=paste0(Figuresdir, "Scatter_metazs_survival_vs_survival_adj_age_sex_significant_genes")
pdf(file=paste0(file,'.pdf',sep=''), height = 4.5,  width = 4.5, family = "Helvetica")
fm %>% dplyr::filter(!is.na(survival.direction)) %>% 
  with(smoothScatter(metaz_survival_gene_adj_age_sex, metaz.coxph, xlab="Meta-z score (Adjusted for age & sex)", ylab="Meta-z score"))
fm %>% dplyr::filter(!is.na(survival.direction)) %>% 
  with(abline(lm(metaz.coxph~metaz_survival_gene_adj_age_sex), col="red", lty=2))
dev.off()

fm %>% dplyr::filter(!is.na(survival.direction)) %>% 
  with(cor.test(metaz.coxph, metaz_survival_gene_adj_age_sex))
#cor=0.988171 , p-value < 2.2e-16

p=PRECOG.HNSCC.All %>% 
  ggplot2::ggplot(aes(metaz_survival_gene_adj_age_sex, metaz.coxph)) + geom_point()  + labs(x="Meta-z score (Adjusted for age & sex)", y="Meta-z score")

file=paste0(Figuresdir, "Scatter_metazs_survival_vs_survival_adj_age_sex")
pdf(file=paste0(file,'.pdf',sep=''), height = 4.5,  width = 4.5, family = "Helvetica")
p
dev.off()

#Boxplot of adjusted p-values by cluster 
dat=fm %>% 
  dplyr::filter(!is.na(survival.direction)) %>% 
  dplyr::mutate(Gene.cluster2=factor(ifelse(is.na(survival.gene.cluster.phenograph),"S7",as.character(survival.gene.cluster.phenograph)))) %>% 
  mutate(dir=factor(ifelse(metaz.coxph>(3.09),"Pro-survival","Anti-survival"), levels=c("Pro-survival","Anti-survival")))

#Boxplot of meta-zs after adjusting for age and sex
lnm.cols=gnuplot_colors[c(8:11,13:14)]
cols2=c(lnm.cols, "lightgrey")

p=ggplot(transform(dat, xjit=jitter(as.numeric(as.factor(Gene.cluster2)))), aes(x=as.factor(Gene.cluster2),y=metaz_survival_gene_adj_age_sex)) +
  geom_boxplot(outlier.shape=NA, aes(fill=Gene.cluster2)) + geom_point(size=0.8, aes(x=xjit)) + 
  labs( x="Gene cluster", y="Meta-z score (Adjusted for age + sex)") + 
  scale_fill_manual(values=cols2, name="LNM gene cluster", labels=c("S7"="No cluster")) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

file=paste0(Figuresdir, "jitter_clusters_survival_adjusted_age_sex")
pdf(file=paste0(file,'.pdf',sep=''), height = 6,  width = 8, family = "Helvetica")
p
dev.off()


##############################################################################################################################################
#Meta-analysis of genes associated with LNM adjusted for age and sex 
##############################################################################################################################################

#Get studies with enough LNM cases, women and men, age 
alldata2=readRDS("~/Documents/Projects/HNSCC_PreCog/data/Precog.HNSCC.all.exp.clin.data.curated.for.PRECOG.rds")

explist=alldata2$explist
infolist=alldata2$infolist

node_objects=readRDS(paste(Datadir, "Precog.HNSCC.lymph.node.metastasis.objects.all.HNSCC.rds", sep=""))

infolist=infolist[names(node_objects)]

accessions=names(infolist)

#restrict to studies with continuous age data 
infolist=infolist[lapply(infolist, function(x) is.numeric(x$COVAR_age)) %>% unlist %>% which() %>% names()]

tab.vars=lapply(infolist, function(x) table(x$COVAR_sex)) %>% abind::abind(along=2) %>% t %>% as.data.frame()
tab.vars$n.age=lapply(infolist, function(x) length(which(!is.na(x)))) %>% unlist
tab.vars$n.node=lapply(infolist, function(x) table(x$COVAR_N_status)) %>% abind::abind(along=2) %>% t %>% as.data.frame()

accessions.node=tab.vars %>% filter(female>5 & male>5) %>% rownames()
#"E-MTAB-1328" "GSE10121"    "GSE23558"    "GSE3292"     "GSE39366"    "GSE41116"    "GSE65858"    "GSE686"      "GSE78060"    "GSE85195"    "GSE95805"    "GSE9844"    "Thurlow"     "TCGA" 

accessions.node=setdiff(accessions.node, "GSE686") 
length(accessions.node) #13

lapply(infolist[accessions.node], nrow) %>% unlist %>% sum
#1316 patients

source(paste(RscriptsPath, "coxph.apply.R", sep=""))

####
allglms=list()
for(i in 1:length(accessions.node)){
  acc=accessions.node[i]
  exp=as.matrix(explist[[acc]])
  info=infolist[[acc]]
  
  glms=apply(exp,1, glm_node_adj_age_sex)
  
  glms=glms[which(lapply(glms, length)==12)]
  glms=as.data.frame(abind(glms, along=1))
  glms=glms[!is.na(glms$z_gene),]
  
  glms=glms %>% dplyr::mutate_at(c('z_gene', 'z_age', 'N', 'z_sex', 'p_gene', 'p_age', 'p_sex', 'se_gene', 'se_age', 'se_sex'), as.numeric)
  glms$q_gene=p.adjust(glms$p_gene, method="fdr")
  glms$q_age=p.adjust(glms$p_age, method="fdr")
  glms$q_sex=p.adjust(glms$p_sex, method="fdr")
  glms$gene=rownames(glms)
  
  allglms[[i]]=glms
}
names(allglms)=accessions.node

saveRDS(allglms, paste(Resultsdir, "Precog.HNSCC.glm.lnm.adjusted.age.sex.rds", sep=""))
allglms=readRDS(paste(Resultsdir, "Precog.HNSCC.glm.lnm.adjusted.age.sex.rds", sep=""))

####
#Perform meta-analysis of prognostic scores using the Liptak method 
library(abind)
allglms.abind=as.data.frame(abind(allglms, along=1))
rownames(allglms.abind)=NULL

allglms.abind=allglms.abind %>% dplyr::mutate_at(c('z_gene', 'z_age', 'N', 'z_sex', 'p_gene', 'p_age', 'p_sex', 'se_gene', 'se_age', 'se_sex'), as.numeric)

allglms.abind %>% dplyr::select("DatasetID", "N") %>% 
  filter(DatasetID %in% names(allglms)) %>% unique %>% pull(N) %>% as.numeric() %>% sum
#1181 patient samples used

#restricting to genes that found within at least two studies, otherwise its not a meta-analysis
tab=table(allglms.abind$gene)
meta.genes=names(tab[which(tab>=2)])
#23542 genes for which at least two studies have the gene. Running meta analysis on all of these 

meta.genes=meta.genes %>% 
  as.data.frame() %>% 
  purrr::set_names("gene") 

meta.genes$metaz.gene=lapply(meta.genes$gene, function(gene){
  mdf=allglms.abind[allglms.abind$gene==gene,]
  metaz=sum(mdf$N*mdf$z_gene)/sqrt(sum(mdf$N^2))
  return(metaz)
}) %>% unlist

meta.genes$metaz.age=lapply(meta.genes$gene, function(gene){
  mdf=allglms.abind[allglms.abind$gene==gene,]
  metaz=sum(mdf$N*mdf$z_age)/sqrt(sum(mdf$N^2))
  return(metaz)
}) %>% unlist

meta.genes$metaz.sex=lapply(meta.genes$gene, function(gene){
  mdf=allglms.abind[allglms.abind$gene==gene,]
  metaz=sum(mdf$N*mdf$z_sex)/sqrt(sum(mdf$N^2))
  return(metaz)
}) %>% unlist

colnames(meta.genes)=c("gene","metaz_lnm_gene_adj_age_sex", "metaz_lnm_age_adj_age_sex", "metaz_lnm_sex_adj_age_sex")

saveRDS(meta.genes, paste0(Resultsdir, "Precog.HNSCC.lnm.adjusted.age.sex.meta.stats.rds"))
meta.genes=readRDS(paste0(Resultsdir, "Precog.HNSCC.lnm.adjusted.age.sex.meta.stats.rds"))



fm=fm %>% dplyr::left_join(meta.genes, by="gene")

saveRDS(fm, paste0(Datadir, "Findmarkers_Puram_U54_allcombined.update.061323.rds"))

lnm.zscores.all=readRDS(paste0(Datadir, "zscores_lnm_all_genes.rds"))
lnm.zscores.all$metaz_lnm_gene_adj_age_sex=meta.genes[match(lnm.zscores.all$gene, meta.genes$gene),"metaz_lnm_gene_adj_age_sex"]

plot(lnm.zscores.all$metaz_lnm_gene_adj_age_sex, lnm.zscores.all$zval)
cor.test(lnm.zscores.all$metaz_lnm_gene_adj_age_sex, lnm.zscores.all$zval)
#cor=0.8147688 , p-value < 0.00000000000000022

file=paste0(Figuresdir, "Scatter_metazs_lnm_vs_lnm_adj_age_sex")
pdf(file=paste0(file,'.pdf',sep=''), height = 4.5,  width = 4.5, family = "Helvetica")
smoothScatter(lnm.zscores.all$metaz_lnm_gene_adj_age_sex, lnm.zscores.all$zval, xlab="Meta-z score (Adjusted for age & sex)", ylab="Meta-z score")
abline(lm(lnm.zscores.all$zval~lnm.zscores.all$metaz_lnm_gene_adj_age_sex), lty=2, col="red")
dev.off()

fm %>% dplyr::filter(!is.na(lnm.direction)) %>% 
  with(cor.test(metaz.lnm, metaz_lnm_gene_adj_age_sex))
#cor=0.9845168  , p-value < 2.2e-16

#Smooth scatter and correlations just for survival-associated genes 
file=paste0(Figuresdir, "Scatter_metazs_lnm_vs_lnm_adj_age_sex_significant_genes")
pdf(file=paste0(file,'.pdf',sep=''), height = 4.5,  width = 4.5, family = "Helvetica")
fm %>% dplyr::filter(!is.na(lnm.direction)) %>% 
  with(smoothScatter(metaz_lnm_gene_adj_age_sex, metaz.lnm, xlab="Meta-z score (Adjusted for age & sex)", ylab="Meta-z score"))
fm %>% dplyr::filter(!is.na(lnm.direction)) %>% 
  with(abline(lm(metaz.lnm~metaz_lnm_gene_adj_age_sex), col="red", lty=2))
dev.off()

p=lnm.zscores.all %>% 
  ggplot2::ggplot(aes(metaz_lnm_gene_adj_age_sex, zval)) + geom_point()  + labs(x="Meta-z score (Adjusted for age & sex)", y="Meta-z score")

file=paste0(Figuresdir, "Scatter_metazs_lnm_vs_lnm_adj_age_sex")
pdf(file=paste0(file,'.pdf',sep=''), height = 4.5,  width = 4.5, family = "Helvetica")
p
dev.off()

#Make boxplot of age and sex-adjusted meta-zscores for LNM gene clusters 
dat=fm %>% 
  dplyr::filter(!is.na(lnm.direction)) %>% 
  dplyr::mutate(Gene.cluster2=factor(ifelse(is.na(lnm.gene.cluster.phenograph),"S7",as.character(lnm.gene.cluster.phenograph)))) %>% 
  mutate(dir=factor(ifelse(metaz.lnm>(3.09),"Pro-LNM","Anti-LNM"), levels=c("Pro-LNM","Anti-LNM")))

#Add top ten genes in terms of meta-z score ranks
#Boxplot of meta-zs after adjusting for age and sex
lnm.cols=gnuplot_colors[c(8:11,13:14)]
cols2=c(lnm.cols, "lightgrey")

p=ggplot(transform(dat, xjit=jitter(as.numeric(as.factor(Gene.cluster2)))), aes(x=as.factor(Gene.cluster2),y=metaz_lnm_gene_adj_age_sex)) +
  geom_boxplot(outlier.shape=NA, aes(fill=Gene.cluster2)) + geom_point(size=0.8, aes(x=xjit)) + 
  labs( x="Gene cluster", y="Meta-z score (Adjusted for age + sex)") + 
  scale_fill_manual(values=cols2, name="LNM gene cluster", labels=c("S7"="No cluster")) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

file=paste0(Figuresdir, "jitter_clusters_lnm_adjusted_age_sex")
pdf(file=paste0(file,'.pdf',sep=''), height = 6,  width = 8, family = "Helvetica")
p
dev.off()


####################################################################################################
#Meta-analysis of survival-associated genes in HPV negative HNC only (Excluding HPV+ve OPH or possible HPV+ve OPH)
####################################################################################################

#Need to restrict to HPV negative for studies in which HPV status us known, then inlcude all cases in studies that can't contain HPV OPH
alldata2=readRDS("~/Documents/Projects/HNSCC_PreCog/data/Precog.HNSCC.all.exp.clin.data.curated.for.PRECOG.rds")

explist=alldata2$explist
infolist=alldata2$infolist

Survival_objects=readRDS(paste(Resultsdir, "Precog.HNSCC.survival.objects.all.HNSCC.rds", sep=""))
infolist=infolist[names(Survival_objects)]

#Restrict to studies that either don't include HPV positive OPH or those that do with annotation 
HPV.stat.studies=readRDS(paste0(Datadir, "HPV.status.studies.meta.analysis.rds"))

HPV.studies=c("GSE39366", "GSE65858", "Thurlow",  "TCGA") #These were identified from an earlier analysis as being studies with HPV and survival data 
No.HPV.studies=names(HPV.stat.studies[HPV.stat.studies %in% c("All.negative")])

#Restrict HPV.studies to HPV negative
infolist_hpv_studies=lapply(infolist[HPV.studies], function(x) x[x$COV_HNSCC_HPV_status=="negative",])
infolist=c(infolist[No.HPV.studies], infolist_hpv_studies)
infolist=infolist[intersect(names(infolist), names(Survival_objects))]
#Find out the min number of patients for survival analysis and restrict all analyses to those studies...

Survival_objects.surv=list()
infolist.surv=list()
explist.surv=list()

#Go through and restict to studies with at least 20 samples with survival data 
accessions=names(infolist)
for(acc in accessions){
  Survival_object=Survival_objects[[acc]]
  info=infolist[[acc]]
  exp=explist[[acc]]
  OverlapSamples=intersect(rownames(info), names(Survival_object))
  
  Survival_objects.surv[[acc]]=Survival_object[OverlapSamples]
  infolist.surv[[acc]]=info[OverlapSamples,]
  explist.surv[[acc]]=exp[,OverlapSamples]
}

keep=names(which(unlist(lapply(Survival_objects.surv, function(x) length(x)>=20))))
"GSE31056"        "GSE41613"        "GSE85195"        "GSE95805"        "GSE23558"        "GSE2379_GPL8300" "GSE41116"        "E-MTAB-1328"     "GSE27020"       "GSE39366"        "TCGA"            "GSE65858"        "Thurlow"     
#13 studies 

#Restrict to studies with at least 20 samples with survival data after removing HPV positive 
Survival_objects=Survival_objects.surv[keep]
explist=explist.surv[keep]
infolist=infolist.surv[keep]

#Now run the meta-analysis in the way that you ran the original one
accessions.survival=names(Survival_objects)
#OK 13 studies with over 20 or more samples with both gene expression data and survival

lapply(Survival_objects, length) %>% unlist %>% sum

#
RscriptsPath="~/Documents/scripts/rscripts2/"
source(paste(RscriptsPath, "coxph.apply.R", sep=""))

allcoxph=list()
for(i in 1:length(accessions.survival)){
  
  acc=accessions.survival[i]
  Survival_object=Survival_objects[[acc]]
  exp=explist[[acc]]
  
  OverlapSamples=intersect(names(Survival_object), colnames(exp))
  Survival_object=Survival_object[OverlapSamples]
  exp=exp[,OverlapSamples]
  
  survs=apply(exp,1, coxph.apply.provide.survobject.with.pvalue)
  survs=survs[which(lapply(survs, length)==5)]
  survs=as.data.frame(abind(survs, along=1))
  survs=survs[!is.na(survs$Z),]
  survs$Z=as.numeric(as.character(survs$Z))
  survs$P=as.numeric(as.character(survs$P))
  survs$Q=p.adjust(survs$P, method="fdr")
  survs$gene=rownames(survs)
  
  allcoxph[[i]]=survs
}
names(allcoxph)=accessions.survival

saveRDS(allcoxph, paste0(Resultsdir, "Precog.HNSCC.coxph.survival.HPV.negative.rds"))
allcoxph=readRDS(paste0(Resultsdir, "Precog.HNSCC.coxph.survival.HPV.negative.rds"))

#Now to meta-analysis to combine summary statistics for each study
library(abind)
allcoxph.abind=as.data.frame(abind(allcoxph, along=1))
rownames(allcoxph.abind)=NULL
allcoxph.abind=allcoxph.abind %>% dplyr::mutate_at(c('Z', 'N', 'SE', 'P', 'Q'), as.numeric)

#restricting to genes that found within at least two studies, otherwise its not a meta-analysis
tab=table(allcoxph.abind$gene)
meta.genes=names(tab[which(tab>=2)])
length(meta.genes) #23525 genes for which at least two studies have the gene. Running meta analysis on all of these 

metaz.coxph=lapply(meta.genes, function(x)  Liptak_combine_z(gene=x, df=allcoxph.abind)) %>% 
  purrr::set_names(meta.genes) %>% unlist

saveRDS(metaz.coxph, paste0(Resultsdir,"meta.z.scores.survival.HPV.negative.rds"))

metaz.coxph.annot=readRDS(paste0(Resultsdir, "Precog.HNSCC.coxph.Liptak.metaz.allHNSCC.updated.06_09_2023.rds"))
metaz.coxph.annot$metaz.coxph.negative=metaz.coxph.neg[match(metaz.coxph.annot$gene, names(metaz.coxph.neg))]

saveRDS(metaz.coxph.annot, paste0(Resultsdir, "Precog.HNSCC.coxph.Liptak.metaz.allHNSCC.updated.06_09_2023.rds"))

#add survival adjusted for HPV to the fm object 
fm=readRDS(paste0(Datadir, "Findmarkers_Puram_U54_allcombined.update.061323.rds"))

fm$metaz.coxph.HPV.negative=metaz.coxph[match(fm$gene, names(metaz.coxph))]

saveRDS(fm, paste0(Datadir, "Findmarkers_Puram_U54_allcombined.update.061323.rds"))

PRECOG.HNSCC.All=readRDS(paste0(Resultsdir, "Precog.HNSCC.coxph.Liptak.metaz.allHNSCC.updated.06_09_2023.rds"))

#PRECOG.HNSCC.All$metaz_survival_gene_adj_age_sex=meta.genes[match(PRECOG.HNSCC.All$gene, meta.genes$gene),"metaz_survival_gene_adj_age_sex"]

#Plot Z scores for all samples against z-scores for HPV negative analysis, i.e. samples that don't include HPV+ve OPH

PRECOG.HNSCC.All$metaz.coxph.HPV.negative=metaz.coxph[match(PRECOG.HNSCC.All$gene, names(metaz.coxph))]

plot(PRECOG.HNSCC.All$metaz.coxph.HPV.negative, PRECOG.HNSCC.All$metaz.coxph)
cor.test(PRECOG.HNSCC.All$metaz.coxph.HPV.negative, PRECOG.HNSCC.All$metaz.coxph)
#cor=0.8423222, p-value < 0.00000000000000022

file=paste0(Figuresdir, "Scatter_metazs_survival_all_vs_HPV_neg")
pdf(file=paste0(file,'.pdf',sep=''), height = 4.5,  width = 4.5, family = "Helvetica")
smoothScatter(PRECOG.HNSCC.All$metaz.coxph.HPV.negative, PRECOG.HNSCC.All$metaz.coxph, xlab="Meta-z score (HNC excluding HPV+ve OPC)", ylab="Meta-z score (All HNC)")
dev.off()

#Smooth scatter and correlations just for survival-associated genes 
file=paste0(Figuresdir, "Scatter_metazs_survival_all_vs_HPV_neg_significant_genes")
pdf(file=paste0(file,'.pdf',sep=''), height = 4.5,  width = 4.5, family = "Helvetica")
fm %>% dplyr::filter(!is.na(survival.direction)) %>% 
  with(smoothScatter(metaz.coxph.HPV.negative, metaz.coxph,  xlab="Meta-z score (HNC excluding HPV+ve OPC)", ylab="Meta-z score (All HNC)"))
fm %>% dplyr::filter(!is.na(survival.direction)) %>% 
  with(abline(lm(metaz.coxph~metaz.coxph.HPV.negative), col="red", lty=2))
dev.off()

fm %>% dplyr::filter(!is.na(survival.direction)) %>% 
  with(cor.test(metaz.coxph, metaz.coxph.HPV.negative))
#cor=0.9618233  , p-value < 2.2e-16

p=PRECOG.HNSCC.All %>% 
  ggplot2::ggplot(aes(metaz.coxph.HPV.negative, metaz.coxph)) + geom_point()  + labs(x="Meta-z score (HNC excluding HPV+ve OPC)", y="Meta-z score (All HNC)")

file=paste0(Figuresdir, "Scatter_metazs_survival_all_vs_HPV_neg_ggplot")
pdf(file=paste0(file,'.pdf',sep=''), height = 4.5,  width = 4.5, family = "Helvetica")
p
dev.off()

#Boxplots of meta-z scores in HPV negative cases only
dat=fm %>% 
  dplyr::filter(!is.na(survival.direction)) %>% 
  dplyr::mutate(Gene.cluster2=factor(ifelse(is.na(survival.gene.cluster.phenograph),"S7",as.character(survival.gene.cluster.phenograph)))) %>% 
  mutate(dir=factor(ifelse(metaz.coxph>(3.09),"Pro-survival","Anti-survival"), levels=c("Pro-survival","Anti-survival")))

#Boxplot of meta-zs after adjusting for age and sex
lnm.cols=gnuplot_colors[c(8:11,13:14)]
cols2=c(lnm.cols, "lightgrey")

p=ggplot(transform(dat, xjit=jitter(as.numeric(as.factor(Gene.cluster2)))), aes(x=as.factor(Gene.cluster2),y=metaz.coxph.HPV.negative)) +
  geom_boxplot(outlier.shape=NA, aes(fill=Gene.cluster2)) + geom_point(size=0.8, aes(x=xjit)) + 
  labs( x="Gene cluster", y="Meta-z score") + 
  scale_fill_manual(values=cols2, name="Survival gene cluster", labels=c("S7"="No cluster")) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

file=paste0(Figuresdir, "jitter_clusters_survival_hpv_negative")
pdf(file=paste0(file,'.pdf',sep=''), height = 6,  width = 8, family = "Helvetica")
p
dev.off()

####################################################################################################
#Meta-analysis of LNM-associated genes in HPV negative HNC only (Excluding HPV+ve OPH or possible HPV+ve OPH)
####################################################################################################

#Need to restrict to HPV negative for studies in which HPV status us known, then inlcude all cases in studies that can't contain HPV OPH
alldata2=readRDS("~/Documents/Projects/HNSCC_PreCog/data/Precog.HNSCC.all.exp.clin.data.curated.for.PRECOG.rds")

explist=alldata2$explist
infolist=alldata2$infolist

#get lymph node objects
node_objects=readRDS(paste(Datadir, "Precog.HNSCC.lymph.node.metastasis.objects.all.HNSCC.rds", sep=""))
infolist=infolist[names(node_objects)]

#Restrict to studies that either don't include HPV positive OPH or those that do with annotation 
HPV.stat.studies=readRDS(paste0(Datadir, "HPV.status.studies.meta.analysis.rds"))

HPV.studies=c("GSE33205", "GSE39366", "GSE65858", "Thurlow",  "TCGA")  #These were identified from an earlier analysis as being studies with HPV and lnm data 
No.HPV.studies=names(HPV.stat.studies[HPV.stat.studies %in% c("All.negative")])

#Restrict HPV.studies to HPV negative
infolist_hpv_studies=lapply(infolist[HPV.studies], function(x) x[x$COV_HNSCC_HPV_status=="negative",])
infolist=c(infolist[No.HPV.studies], infolist_hpv_studies)
infolist=infolist[intersect(names(infolist), names(Survival_objects))]
#Find out the min number of patients for survival analysis and restrict all analyses to those studies...

#Survival_objects.surv=list()
infolist.lnm=list()
explist.lnm=list()

#Go through and restict to studies with at least 20 samples with survival data 
accessions=names(infolist)
for(acc in accessions){
  info=infolist[[acc]]
  exp=explist[[acc]]
  
  OverlapSamples=intersect(rownames(info), colnames(exp))
  infolist.lnm[[acc]]=info[OverlapSamples,]
  explist.lnm[[acc]]=exp[,OverlapSamples]
}

keep=lapply(infolist.lnm, function(x) table(x$COVAR_N_status)) %>% 
  abind::abind(along=2) %>% t %>% as.data.frame %>% 
  purrr::set_names("N0","N1") %>% 
  dplyr::filter(N0>=5 & N1>=5) %>% 
  rownames()
#"GSE85195"    "GSE95805"    "GSE23558"    "GSE41116"    "E-MTAB-1328" "GSE39366"    "GSE65858"    "Thurlow"     "TCGA" 
#9 studies 

#Restrict to studies with at least 20 samples with survival data after removing HPV positive 
explist.node=explist.lnm[keep]
infolist.node=infolist.lnm[keep]

lapply(explist.node, ncol) %>% unlist %>% sum
#948 patient samples 

#Getting genes within fm 
#fm=readRDS(paste0(Datadir, "Findmarkers_Puram_U54_allcombined.update.100521.rds"))
#allgenes=fm$gene
allgenes=unique(unlist(lapply(explist.node, rownames)))
allgeneinstances=unlist(lapply(explist.node, rownames))

#Find the number of studies for each gene is available 
nstudies=lapply(allgenes, function(x) length(grep(paste("\\b",x, "\\b", sep=""), allgeneinstances)))
names(nstudies)=allgenes
#check a few randomly and they're all fine
#saveRDS(nstudies, paste0(Datadir, "nstudies.ge.rds"))
saveRDS(nstudies, paste0(Datadir, "nstudies.ge.HPVneg.rds"))

num=round(length(accessions.node)*.5)
genes.available.all=names(unlist(nstudies)[unlist(nstudies)>=num])
saveRDS(genes.available.all, paste(Resultsdir,"Precog.meta.analysis.genes.HPV.neg.rds", sep=""))

#16899 genes in at least half (n=10) studies 
#genes.available.all=fm$gene

source(paste(RscriptsPath, "meta.analysis.mean.difference.functions.R", sep=""))
node.list=make.meta.node_multistudies(clinlist=infolist.node, expressionlist=explist.node, studyaccessions=names(infolist.node), genes=genes.available.all, varfactor="COVAR_N_status")

saveRDS(node.list, paste(Resultsdir,"Precog.meta.analysis.node.status.HPV.negative.rds", sep=""))
node.list=readRDS(paste(Resultsdir,"Precog.meta.analysis.node.status.HPV.negative.rds", sep=""))

###
meta2=node.list$meta

metatab=cbind(
  unlist(lapply(meta2, function(x) x$pval)),
  unlist(lapply(meta2, function(x) x$b)),
  unlist(lapply(meta2, function(x) x$zval)))

metatab=as.data.frame(metatab)
dim(metatab[metatab$V3>3.09,])
dim(metatab[metatab$V3<(-3.09),])
colnames(metatab)=c("pval","beta","zval")
metatab$gene=rownames(metatab)

saveRDS(metatab, paste0(Datadir, "zscores_lnm_all_genes.HPV.negative.rds"))
metatab=readRDS(paste0(Datadir, "zscores_lnm_all_genes.HPV.negative.rds"))

fm=readRDS(paste0(Datadir, "Findmarkers_Puram_U54_allcombined.update.061323.rds"))
fm$metaz.lnm.HPV.negative=metatab[match(fm$gene, metatab$gene),"zval"]

saveRDS(fm, paste0(Datadir, "Findmarkers_Puram_U54_allcombined.update.061323.rds"))

fm %>% dplyr::filter(!is.na(lnm.direction)) %>% 
  with(cor.test(metaz.lnm, metaz.lnm.HPV.negative))
#cor=0.9684869  , p-value < 2.2e-16

#Smooth scatter and correlations just for survival-associated genes 
file=paste0(Figuresdir, "Scatter_metazs_lnm_all_vs_HPV_neg")
pdf(file=paste0(file,'.pdf',sep=''), height = 4.5,  width = 4.5, family = "Helvetica")
fm %>% dplyr::filter(!is.na(lnm.direction)) %>% 
  with(smoothScatter(metaz.lnm.HPV.negative, metaz.lnm, xlab="Meta-z score (HNC excluding HPV+ve OPC)", ylab="Meta-z score (All HNC)"))
fm %>% dplyr::filter(!is.na(lnm.direction)) %>% 
  with(abline(lm(metaz.lnm~metaz.lnm.HPV.negative), col="red", lty=2))
dev.off()


#Boxplot by gene cluster
dat=fm %>% 
  dplyr::filter(!is.na(lnm.direction)) %>% 
  dplyr::mutate(Gene.cluster2=factor(ifelse(is.na(lnm.gene.cluster.phenograph),"S7",as.character(lnm.gene.cluster.phenograph)))) %>% 
  mutate(dir=factor(ifelse(metaz.lnm>(3.09),"Pro-LNM","Anti-LNM"), levels=c("Pro-LNM","Anti-LNM")))

#Add top ten genes in terms of meta-z score ranks
#Boxplot of meta-zs after adjusting for age and sex
lnm.cols=gnuplot_colors[c(8:11,13:14)]
cols2=c(lnm.cols, "lightgrey")

p=ggplot(transform(dat, xjit=jitter(as.numeric(as.factor(Gene.cluster2)))), aes(x=as.factor(Gene.cluster2),y=metaz.lnm.HPV.negative)) +
  geom_boxplot(outlier.shape=NA, aes(fill=Gene.cluster2)) + geom_point(size=0.8, aes(x=xjit)) + 
  labs( x="Gene cluster", y="Meta-z score") + 
  scale_fill_manual(values=cols2, name="LNM gene cluster", labels=c("S7"="No cluster"))+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

file=paste0(Figuresdir, "jitter_clusters_lnm_hpv_negative")
pdf(file=paste0(file,'.pdf',sep=''), height = 6,  width = 8, family = "Helvetica")
p
dev.off()












































