
###########################################################################################
#Processing and analyzing the GSE26549 (Saintigny et al.) dataset
############################################################################################

library(affy)
library(preprocessCore)
library(WGCNA)
library(GEOquery)

setwd("/oak/stanford/groups/andrewg/users/kbren/projects/Puram/data/")
expdir="/oak/stanford/groups/andrewg/users/kbren/projects/Puram/data/"

install.packages(file.path(expdir, "hugene10sthsentrezg.db_25.0.0.tar.gz"), repos = NULL)
install.packages(file.path(expdir, "hugene10sthsentrezgcdf_25.0.0.tar.gz"), repos = NULL)
install.packages(file.path(expdir, "hugene10sthsentrezgprobe_25.0.0.tar.gz"), repos = NULL)

##
celdir="/data/GSE26549/"

gse="GSE26549"
filePaths = getGEOSuppFiles(gse, baseDir = expdir)

#processing all hgu133plus2 datasets
library(oligo)
library(hugene10sthsentrezg.db)
library(hugene10sthsentrezgcdf)
library(hugene10sthsentrezgprobe)

celdir="/data/GSE26549/GSE26549_RAW.tar"
untar(celdir)
celdir="/data/GSE26549/"

data=affy::ReadAffy(celfile.path=celdir, cdfname="hugene10sthsentrezgcdf") 
eset=affy::mas5(data)

###
mas5.ALL=exprs(eset)
colnames(mas5.ALL)=gsub(" ","", colnames(mas5.ALL))
colnames(mas5.ALL)=gsub("[.].*","", colnames(mas5.ALL))
colnames(mas5.ALL)=gsub("_.*","", colnames(mas5.ALL))

#Remove Affy control probes, which are at the end of the document and mess with the adding gene symbols to the rows
#mas5.ALL=mas5.ALL[-grep("AFFX-", rownames(mas5.ALL)),]
#Format values to 5 decimal places
mas5.ALL=format(mas5.ALL, digits=5)
mas5.ALL=as.matrix(mas5.ALL)

mas5.ALL2=apply(mas5.ALL, 1, function(x) as.numeric(as.character(x)))
mas5.ALL2=as.data.frame(t(mas5.ALL2))
colnames(mas5.ALL2)=colnames(mas5.ALL)

#get clinical info 
library(GEOquery)
gse="GSE26549"
Sys.setenv(VROOM_CONNECTION_SIZE=500072)
getSet=getGEO(gse)
info=pData(getSet[[1]])

setdiff(rownames(info), colnames(mas5.ALL2))
setdiff(colnames(mas5.ALL2), rownames(info))

info=info[match(colnames(mas5.ALL2), rownames(info)),]
#All are oral luekoplakia so I guess it is safe to quantile normalize
#quantile normalize
mas5.ALL.qn=preprocessCore::normalize.quantiles(as.matrix(mas5.ALL2))
dimnames(mas5.ALL.qn)=dimnames(mas5.ALL2)

#collapse probes to genes 
probes.ALL=row.names(mas5.ALL)
entrezIDmap=hugene10sthsentrezgENTREZID
entrezID.ALL = unlist(mget(probes.ALL, entrezIDmap))
ID.ALL = unlist(mget(probes.ALL, entrezIDmap))
collapserows.annot=as.data.frame(cbind(ID.ALL, entrezID.ALL))
collapserows.annot$probe=rownames(collapserows.annot)
collapserows.annot=na.omit(collapserows.annot[,c("entrezID.ALL", "probe")])

mas5.collapse=WGCNA::collapseRows(datET=mas5.ALL.qn, rowGroup=collapserows.annot$entrezID.ALL, rowID=collapserows.annot$probe)
mas5.c=mas5.collapse$datETcollapsed
mas5.clog2=log2(mas5.c)

#Now z-score the genes
mas5.z=apply(mas5.clog2, 1, function(x)  (x-mean(x))/sd(x))
mas5.z=as.data.frame(t(mas5.z))
all(colnames(mas5.z)==rownames(info))

allinfo=list(info=info, mas5.z=mas5.z)
saveRDS(allinfo, paste0(celdir, "processed_exp_z_scores_", gse, ".rds"))

###############################################################################
#Analyzing LNM genes in OPLs
###############################################################################

library(RColorBrewer)

allinfo=readRDS("~/Documents/Projects/HNSCC_PreCog/data/processed_exp_z_scores_GSE26549.rds")
exp=allinfo$mas5.z
info=allinfo$info

#fm=readRDS(paste0(Datadir, "Findmarkers_Puram_U54_allcombined.update.100521.rds"))
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

info$histology=plyr::revalue(factor(info$`histology at baseline (breakdown):ch1`), c("NA"=NA))
info$OCFS_Time=as.numeric(as.character(info$`oral cancer-free survival time (years):ch1`))*12
info$OS_Status=as.numeric(factor(info$`outcome:ch1`))

info$`oral cancer-free survival time (years):ch1`=as.numeric(info$`oral cancer-free survival time (years):ch1`)

df=cbind(info[which(info$`oral cancer-free survival time (years):ch1`<=5),c("OCFS_Time", "OS_Status")])

library(survival)
Survival_Object=Surv(df$OCFS_Time,df$OS_Status)

sigs=c(
  makesig(genetype="gene",group="lnm.direction"),
  makesig(genetype="gene",group="lnm.gene.cluster.phenograph"))

sigs=lapply(sigs, function(x) x[x %in% rownames(exp)])

sigtab=as.data.frame(abind::abind(lapply(sigs, function(sig) colMeans(exp[sig,])), along=2))
all(rownames(sigtab)==rownames(info))

sigtab=cbind(sigtab, info)

cols=RColorBrewer::brewer.pal(4,"YlOrRd")

t1=sigtab[,c("Anti-LNM", "histology")]
colnames(t1)=c("Exp","Histology")
t1$Signature=rep("Anti-LNM", nrow(t1))

t2=sigtab[,c("Pro-LNM", "histology")]
colnames(t2)=c("Exp","Histology")
t2$Signature=rep("Pro-LNM", nrow(t2))

tab=as.data.frame(rbind(t1, t2))
tab=na.omit(tab)

p=ggplot(tab, aes(x=Histology, y=Exp)) + geom_boxplot(outlier.shape=NA, aes(fill=Histology)) + geom_jitter(shape=16, size=0.6, position=position_jitter(0.2)) + facet_grid(~Signature)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + labs(x = "Histology (Premalignant disease stage)", y="Scaled mean gene expression")+ scale_fill_manual(name = "Histology", values = c("#FFFFB2", "#FECC5C", "#FD8D3C", "#E31A1C"), labels = c("hyperplasia" = "Hyperplasia", "mild dysplasia" = "Mild dysplasia", "moderate dysplasia" = "Moderate dysplasia", "severe dysplasia"="Severe dysplasia"))+ theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
p=p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white", colour = "black",size = 0.5, linetype = "solid"))

file=paste0(Figuresdir, "boxplot_GSE26549_oral_dysplasia_grade_ggplot_Anti_pro_update.160623")
pdf(file=paste0(file,'.pdf',sep=''), height = 3,  width = 5, family = "Helvetica")
p
dev.off()

#
t1=sigtab[,c("Anti-LNM", "histology")]
colnames(t1)=c("Exp","Histology")
t1$Signature=rep("Anti-LNM", nrow(t1))

t2=sigtab[,c("L4", "histology")]
colnames(t2)=c("Exp","Histology")
t2$Signature=rep("L4", nrow(t2))

tab=as.data.frame(rbind(t1, t2))
tab=na.omit(tab)

tab$Signature2=plyr::revalue(as.factor(tab$Signature), c("L1"="Anti-LNM cluster L1", "L4"="Pro-LNM cluster L4"))

p=ggplot(tab, aes(x=Histology, y=Exp)) + geom_boxplot(outlier.shape=NA, aes(fill=Histology), lwd=0.2) + geom_jitter(shape=16, size=0.4, position=position_jitter(0.2)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + labs(x = "Histology (Premalignant disease stage)", y="Scaled mean gene expression")+ scale_fill_manual(name = "Histology", values = c("#FFFFB2", "#FECC5C", "#FD8D3C", "#E31A1C"), labels = c("hyperplasia" = "Hyperplasia", "mild dysplasia" = "Mild dysplasia", "moderate dysplasia" = "Moderate dysplasia", "severe dysplasia"="Severe dysplasia"))+ theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+ facet_grid(~Signature2)
p=p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white", colour = "black",size = 0.5, linetype = "solid"))

file=paste0(Figuresdir, "boxplot_GSE26549_oral_dysplasia_grade_ggplot_long_labels_update.160623")
pdf(file=paste0(file,'.pdf',sep=''), height = 3,  width = 5, family = "Helvetica")
p
dev.off()

file=paste0(Figuresdir, "boxplot_GSE26549_oral_dysplasia_grade_ggplot_Anti_L4_update.160623")
pdf(file=paste0(file,'.pdf',sep=''), height = 3,  width = 5, family = "Helvetica")
p
dev.off()


t1=sigtab[,c("Anti-LNM", "histology")]
colnames(t1)=c("Exp","Histology")
t1$Signature=rep("Anti-LNM", nrow(t1))

t2=sigtab[,c("L4", "histology")]
colnames(t2)=c("Exp","Histology")
t2$Signature=rep("Pro-LNM", nrow(t2))

tab=as.data.frame(rbind(t1, t2))
tab=na.omit(tab)

p=ggplot(tab, aes(x=Histology, y=Exp)) + geom_boxplot(outlier.shape=NA, aes(fill=Histology)) + geom_jitter(shape=16, size=0.6, position=position_jitter(0.2)) + facet_grid(~Signature)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + labs(x = "Histology (Premalignant disease stage)", y="Scaled mean gene expression")+ scale_fill_manual(name = "Histology", values = c("#FFFFB2", "#FECC5C", "#FD8D3C", "#E31A1C"), labels = c("hyperplasia" = "Hyperplasia", "mild dysplasia" = "Mild dysplasia", "moderate dysplasia" = "Moderate dysplasia", "severe dysplasia"="Severe dysplasia"))+ theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
p=p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white", colour = "black",size = 0.5, linetype = "solid"))

clusts=levels(fm$lnm.gene.cluster.phenograph)

tab=lapply(clusts, function(clust) {sigtab[,c(clust, "histology")] %>% purrr::set_names("Exp","Histology") %>% 
    mutate(Signature=rep(clust, nrow(sigtab)))}) %>% 
  abind::abind(along=1) %>% as.data.frame() %>% na.omit() %>% 
  mutate(Exp=as.numeric(Exp), Histology=factor(Histology, levels=levels(info$histology)), Signature=factor(Signature))


p=ggplot(tab, aes(x=Histology, y=Exp)) + geom_boxplot(outlier.shape=NA, aes(fill=Histology)) + geom_jitter(shape=16, size=0.6, position=position_jitter(0.2)) + facet_wrap(~Signature)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + labs(x = "Histology (Premalignant disease stage)", y="Scaled mean gene expression")+ scale_fill_manual(name = "Histology", values = c("#FFFFB2", "#FECC5C", "#FD8D3C", "#E31A1C"), labels = c("hyperplasia" = "Hyperplasia", "mild dysplasia" = "Mild dysplasia", "moderate dysplasia" = "Moderate dysplasia", "severe dysplasia"="Severe dysplasia"))+ theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
p=p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white", colour = "black",size = 0.5, linetype = "solid"))

p

file=paste0(Figuresdir, "boxplot_GSE26549_oral_dysplasia_grade_ggplot_LNM_clusters")
pdf(file=paste0(file,'.pdf',sep=''), height = 6,  width = 10, family = "Helvetica")
p
dev.off()

sigtab=lapply(names(sigs)[3:8],function(x){summary(glm(as.numeric(sigtab$histology)~sigtab[,x]))$coefficients[c(2,8)]}) %>% 
  abind::abind(along=2) %>% t %>% as.data.frame() %>% purrr::set_names("Estimate","P") %>% 
  dplyr::mutate(Signature=names(sigs)[3:8]) %>% 
  mutate(Q=p.adjust(P, method="fdr"))

saveRDS(sigtab, paste0(Resultsdir,"glm_association_LNM_sigs_with_dysplasia.rds"))
