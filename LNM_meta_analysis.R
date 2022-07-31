
##############################
#Auxiliary scripts for meta-analyses
##############################

library(metafor)
library(esc)

#script to get statistics that are needed to perform random effects model from each study 
Get_meta_vars=function(dataframe=dataframe, var=as.factor(var)){
        dat=array(NA,c(nrow(dataframe),6))
        rownames(dat)=rownames(dataframe)
        colnames(dat)=c("grp1m", "grp1sd", "grp1n", "grp2m", "grp2sd", "grp2n")
        var=as.factor(var)
          
        for( i in 1:nrow(dataframe)){
        dat[i,1]=mean(as.numeric(as.character(dataframe[i,which(var==levels(var)[2])])), na.rm=T)
        dat[i,2]=sd(as.numeric(as.character(dataframe[i,which(var==levels(var)[2])])), na.rm=T)
        dat[i,3]=length(which(!is.na(as.numeric(as.character(dataframe[i,which(var==levels(var)[2])])))))
        dat[i,4]=mean(as.numeric(as.character(dataframe[i,which(var==levels(var)[1])])), na.rm=T)
        dat[i,5]=sd(as.numeric(as.character(dataframe[i,which(var==levels(var)[1])])), na.rm=T)
        dat[i,6]=length(which(!is.na(as.numeric(as.character(dataframe[i,which(var==levels(var)[1])])))))
        } 
          
    dat=as.data.frame(dat, drop=FALSE)
    return(dat)
}


#Wrapper for combine_esc_sd needed to combine statistics acorss studies for a variable number of studies (Between 5 and 21 studies)
combine_esc_KB=function(g){
          if(length(g)==5){
     out=combine_esc(g[[1]], g[[2]], g[[3]], g[[4]], g[[5]])
          } else {
          if(length(g)==6){
     out=combine_esc(g[[1]], g[[2]], g[[3]], g[[4]], g[[5]], g[[6]])
          } else {
          if(length(g)==7){
     out=combine_esc(g[[1]], g[[2]], g[[3]], g[[4]], g[[5]], g[[6]], g[[7]])
          } else {
          if(length(g)==8){
     out=combine_esc(g[[1]], g[[2]], g[[3]], g[[4]], g[[5]], g[[6]], g[[7]], g[[8]])
          } else {
         if(length(g)==9){
     out=combine_esc(g[[1]], g[[2]], g[[3]], g[[4]], g[[5]], g[[6]], g[[7]], g[[8]], g[[9]])
          } else {
         if(length(g)==10){
     out=combine_esc(g[[1]], g[[2]], g[[3]], g[[4]], g[[5]], g[[6]], g[[7]], g[[8]], g[[9]], g[[10]])
          } else {
          if(length(g)==11){
      out=combine_esc(g[[1]], g[[2]], g[[3]], g[[4]], g[[5]], g[[6]], g[[7]], g[[8]], g[[9]], g[[10]], 
          g[[11]])
           } else {
           if(length(g)==12){
      out=combine_esc(g[[1]], g[[2]], g[[3]], g[[4]], g[[5]], g[[6]], g[[7]], g[[8]], g[[9]], g[[10]], 
          g[[11]], g[[12]])
            } else {
            if(length(g)==13){
      out=combine_esc(g[[1]], g[[2]], g[[3]], g[[4]], g[[5]], g[[6]], g[[7]], g[[8]], g[[9]], g[[10]], 
          g[[11]], g[[12]], g[[13]])
             } else {
             if(length(g)==14){
      out=combine_esc(g[[1]], g[[2]], g[[3]], g[[4]], g[[5]], g[[6]], g[[7]], g[[8]], g[[9]], g[[10]], 
          g[[11]], g[[12]], g[[13]], g[[14]])
             } else {
             if(length(g)==15){
      out=combine_esc(g[[1]], g[[2]], g[[3]], g[[4]], g[[5]], g[[6]], g[[7]], g[[8]], g[[9]], g[[10]], 
          g[[11]], g[[12]], g[[13]], g[[14]], g[[15]])
              } else {
              if(length(g)==16){
      out=combine_esc(g[[1]], g[[2]], g[[3]], g[[4]], g[[5]], g[[6]], g[[7]], g[[8]], g[[9]], g[[10]], 
          g[[11]], g[[12]], g[[13]], g[[14]], g[[15]], g[[16]])
               } else {
               if(length(g)==17){
      out=combine_esc(g[[1]], g[[2]], g[[3]], g[[4]], g[[5]], g[[6]], g[[7]], g[[8]], g[[9]], g[[10]], 
          g[[11]], g[[12]], g[[13]], g[[14]], g[[15]], g[[16]], g[[17]])
                } else {
                if(length(g)==18){
      out=combine_esc(g[[1]], g[[2]], g[[3]], g[[4]], g[[5]], g[[6]], g[[7]], g[[8]], g[[9]], g[[10]], 
          g[[11]], g[[12]], g[[13]], g[[14]], g[[15]], g[[16]], g[[17]], g[[18]])
                } else {
                if(length(g)==19){
      out=combine_esc(g[[1]], g[[2]], g[[3]], g[[4]], g[[5]], g[[6]], g[[7]], g[[8]], g[[9]], g[[10]], 
          g[[11]], g[[12]], g[[13]], g[[14]], g[[15]], g[[16]], g[[17]], g[[18]], g[[19]])
                } else {
                if(length(g)==20){
      out=combine_esc(g[[1]], g[[2]], g[[3]], g[[4]], g[[5]], g[[6]], g[[7]], g[[8]], g[[9]], g[[10]], 
          g[[11]], g[[12]], g[[13]], g[[14]], g[[15]], g[[16]], g[[17]], g[[18]], g[[19]], g[[20]])
                } else {
                if(length(g)==21){
      out=combine_esc(g[[1]], g[[2]], g[[3]], g[[4]], g[[5]], g[[6]], g[[7]], g[[8]], g[[9]], g[[10]], 
          g[[11]], g[[12]], g[[13]], g[[14]], g[[15]], g[[16]], g[[17]], g[[18]], g[[19]], g[[20]], g[[21]])
                   }
                  }
                 }
                }
               }
              }
             }
            }
           }
          }
         }
        }
       }
      }
     }
    }
   }
return(out)  
}


#Run random effects meta-analysis acorss studies
make.meta.node_multistudies=function(clinlist, expressionlist, studyaccessions, genes, varfactor){
  
EffectSizes=list()

for(i in 1:length(studyaccessions)){
  
acc=studyaccessions[i] 

exp=expressionlist[[which(names(expressionlist) %in% acc)]]
info=clinlist[[which(names(clinlist) %in% acc)]]

gt=Get_meta_vars(dataframe=exp, var=info[,varfactor])
gt2=lapply(c(1:nrow(gt)), function(k) esc::esc_mean_sd(grp1m = gt[k,1], grp1sd = gt[k,2], grp1n = gt[k,3], grp2m = gt[k,4], grp2sd = gt[k,5], grp2n = gt[k,6], es.type = "g", study = paste("Study ",acc, sep="")))
names(gt2)=rownames(gt)

EffectSizes[[i]]=gt2
}
names(EffectSizes)=studyaccessions

#now perform meta analysis for each gene in genes 
meta.genes=list()

for(f in 1:length(genes)){
gene=genes[f]
#get the indices of studies in which the gene exists
indices=c(which(unlist(lapply(EffectSizes, function(x) gene %in% names(x)))))
#get the effect sizes for studies that have this gene 
genestudies=lapply(EffectSizes[c(indices)], function(x) x[[which(names(x) %in% gene)]])

#Just need to hack this so it will run for any number of studies 
mydat2=combine_esc_KB(genestudies)

#mydat2=combine_esc(lapply(1:length(genestudies), function(x) c(genestudies[[x]])))
#here is a problem with the code. can't figure out how to pass multiple objects combine_esc as list

rma.res=myTryCatch(metafor::rma(yi = es, sei = se, method = "REML", data = mydat2))

if(is.null(rma.res$error)){
  meta.genes[[f]]=rma.res$value
 } else {
  meta.genes[[f]]=rma.res$error
}

}
names(meta.genes)=genes

EffectSizesList=list(meta.genes, EffectSizes)
names(EffectSizesList)=c("meta","effect.sizes")

return(EffectSizesList)
}


#Perfrom coxph to test association of individual gene with survival 
coxph.apply.provide.survobject.with.pvalue=function(x){
  library(survival)
  library(abind)
  
  out <- tryCatch(
  {
  cox=coxph(Survival_object~x)
  z = as.list(coef(cox)/sqrt(diag(vcov(cox))))[[1]]
  se = summary(cox)$coefficients[3]
  n = cox$n
  p=summary(cox)$coefficients[5]
  df = data.frame(Z = z, SE = se, N = n, DatasetID = acc, P=p)
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

#Liptak's meta-Z score to combine per-study z-scores into meta-Z score 
Liptak_combine_z=function(gene, df){
mdf=df[df$gene==gene,]
metaz=sum(mdf$N*mdf$Z, na.rm=T)/sqrt(sum(mdf$N^2, na.rm=T))
return(metaz)
}

#########################################################
##########################################################
#Meta-analysis of survival associated genes 
##########################################################
##########################################################

#Get cleaned lists of gene expression and clinical data for all studies 
alldata=readRDS("~/Documents/Projects/HNSCC_PreCog/data/Precog.HNSCC.all.exp.clin.data.rds")
explist=alldata$explist
zscores.list=explist #expression
infolist=alldata$infolist #clinical data

#Finding the survival measures that are available for analysis within each study 
survlist=list()
for(i in 1:length(infolist)){
  info=infolist[[i]]
  timevars=colnames(info)[grep("_Time", colnames(info))]
  status=colnames(info)[grep("_Status", colnames(info))]
  survars=list(status, timevars)
  names(survars)=c("status", "time")
  print(names(infolist)[i])
  print("OS_Time" %in%  timevars)
  survlist[[i]]=survars
}
names(survlist)=names(infolist)

#make infolist.surv
survlist=survlist[as.numeric(which(unlist(lapply(survlist, function(x) length(x$status)!=0 & length(x$time)!=0))))]

#Find out which datasets have overall survival. For those that don't, find out which alternative survival variables they have 
survvartype=unlist(lapply(survlist, function(x) ifelse("OS_Time" %in% x$time & "OS_Status" %in% x$status, "OS", gsub("_Status","", x$status))))
survvartype[which(names(survvartype)=="GSE39366")]="RFS"

accessions=names(survvartype)
infolist.surv=infolist[accessions]

for(i in 1:length(accessions)){
acc=accessions[i]
print(acc)
survvar=paste(survvartype, "_Status", sep="")[i]
survtimevar=paste(survvartype, "_Time", sep="")[i]
infolist.surv[[which(names(infolist.surv) %in% acc)]][,survvar]=as.numeric(as.character(infolist.surv[[which(names(infolist.surv) %in% acc)]][,survvar]))
infolist.surv[[which(names(infolist.surv) %in% acc)]][,survtimevar]=as.numeric(as.character(infolist.surv[[which(names(infolist.surv) %in% acc)]][,survtimevar]))
}
#

#Restrict the meta-analysis to studies that had at least five patients that were censored and five patients that had an event.
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

#Make survival objects for each study using the survival package
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
saveRDS(Survival_objects, paste(Resultsdir, "Precog.HNSCC.survival.objects.all.HNSCC.rds", sep=""))

infolist.surv=infolist[which(names(infolist) %in% accessions.survival)]
accessions.survival=names(infolist.surv)

#Get survival objects and use these to do the meta analysis 
Survival_objects=readRDS(paste(Datadir, "Precog.HNSCC.survival.objects.all.HNSCC.rds", sep=""))

#Restrict to studies with at least 20 patients with survival data 
keep=names(which(unlist(lapply(Survival_objects, function(x) length(x)>=20))))
#17 studies
Survival_objects=Survival_objects[keep]
accessions.survival=names(Survival_objects)

accessions.survival=intersect(accessions.survival, names(explist))
#Down to 17 studies with both gene expression results and survival data 

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

#Restricting to studies with at least 20 patients with survival and gene expression data
#This is in additon to earlier restrction that there are at least 5 events and 5 censorships within each study 
accessions.survival=names(which(unlist(lapply(coxphlengths, function(x) length(x)>=20))))
#OK 17 studies with over 20 or more samples with both gene expression data and survival

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
saveRDS(allcoxph, paste(Resultsdir, "Precog.HNSCC.coxph.survival.allstudies.allHNSCC.updated.03.30.2010.rds", sep=""))

#Apply Liptak's weighted meta-Z test to per-study z-scores in order to calculate meta-Z score for association of genes with survival across studies
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

#restricting to genes that found within at least two studies
tab=table(allcoxph.abind$gene)
meta.genes=names(tab[which(tab>=2)])

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

#add the number of studies for which each gene was available
metaz.coxph.annot$n.studies=tab[match(metaz.coxph.annot$gene, names(tab))]

#Identify genes that were signficantly associated with survival (Absolute meta-z >=3.09)
threshold="3.09"
metaz.coxph.annot$signficant=ifelse(abs(metaz.coxph.annot$metaz.coxph)>=threshold,"yes","no")

saveRDS(metaz.coxph.annot, paste(Resultsdir, "Precog.HNSCC.coxph.Liptak.metaz.allHNSCC.updated.03_30_2020.rds", sep=""))


#################################################################################
#Unsupervised clustering all signficantly prognostic genes across all datasets with a sufficient number of genes
####################################################################################

pvals.sig=metaz.coxph.annot$gene[metaz.coxph.annot$signficant=="yes"]

#Make a dataframe indicating the studies each survival-associated gene is represented in 
alldata=readRDS("~/Documents/Projects/HNSCC_PreCog/data/Precog.HNSCC.all.exp.clin.data.rds")
explist=alldata$explist
zscores.list=explist
infolist=alldata$infolist

allgenes=lapply(explist, rownames)
#for each dataset get the intersection with the signficant gene list 
allgenes2=lapply(allgenes, function(x) intersect(pvals.sig, x))

genes=unique(unlist(allgenes2))

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

#Now get genes that are present in all of the studies that have 80% of genes (The genes that will be clustered into phenograph clusters)
gene.present.df.allgenes=gene.present.df[,studiesallgenes]
gene.present.df.allgenes=gene.present.df.allgenes[which(apply(gene.present.df.allgenes, 1, function(x) all(x=="yes"))),]

explist.sig=explist[which(names(explist) %in% colnames(gene.present.df.allgenes))]
ex.commongenes=lapply(explist.sig, function(x) x[match(genes, rownames(x)),])
ex.commongenes.df=abind(ex.commongenes, along=2)
ex.commongenes.df=na.omit(ex.commongenes.df)
ex.commongenes.df=ex.commongenes.df

library(Rphenograph)
pheno <- Rphenograph(ex.commongenes.df, k = 45)

library(umap)
clust.umap = umap(ex.commongenes.df)

#Make umap of survival gene phenograph cluster 
saveRDS(clust.umap, paste0(Resultsdir, "UMAP_PRECOG_HNSCC_survival_gene_clusters.rds"))

saveRDS(pheno, paste0(Resultsdir, "Precog.prognostic.phenograph.clusters.Liptak.allHNSCC.rds"))

membership=as.character(membership(Rphenograph_out[[2]]))
df=data.frame(gene=as.character(rownames(ex.commongenes.df)), cluster=membership)

#add phenograph cluster to survival genes results object 
metaz.coxph.annot$gene.clust.phenograph=df[match(metaz.coxph.annot$gene, "cluster"]

saveRDS(metaz.coxph.annot, paste(Resultsdir, "Precog.HNSCC.coxph.Liptak.metaz.allHNSCC.updated.03_30_2020.rds", sep=""))

#########################################################
##########################################################
#Meta-analysis of lymph node metastasis-associated genes 
##########################################################
##########################################################

#Making 'node_objects', a list of vectors indicating the LNM status of each patient within each study
alldata=readRDS(paste(Datadirclin, "Precog.HNSCC.all.exp.clin.data.rds", sep=""))
infolist=alldata$infolist

#First get all gene expression and clinical data
explist=alldata$explist
infolist3=alldata$infolist

accessions=names(infolist3)

#excluding studies with insufficient numbers of genes 
excluded_studies=c("GSE30788_discovery", "GSE30788_validation")
accessions=setdiff(accessions, excluded_studies)

varfactor="COVAR_N_status"

#format node status and node number for all studies
has.node.num=names(which(unlist(lapply(infolist, function(x) !all(is.na(x[,varfactor])) & !all(is.na(x[,varfactor2]))))))

#restrict to studies with at least five LNM+ve and five LNM-ve rpimary HNCs
node.enough=list()
accessions=names(infolist)

for(i in 1:length(accessions)){
acc=accessions[i]
info=infolist[[which(names(infolist) %in% acc)]]
print(acc)
if(varfactor %in% colnames(info)){
  print(summary(as.factor(info[, varfactor])))
  node.enough[[acc]]=length(which(info[, varfactor]==levels(info[, varfactor])[1]))>=5 & length(which(info[, varfactor]==levels(info[, varfactor])[2]))>=5
} else {
  print("no node status")
  node.enough[[acc]]=FALSE
  }
}
accessions.node=names(which(unlist(node.enough)))
infolist.node=infolist[accessions.node]
#21 studies with at least five cases with and without lymph nodes 

#Make lymph node metastasis objects with just primary tumor
#This is a just list of vectors indicating the LNM status of each primary HNC in the study, which can be used for the meta-analysis
accessions.node=names(infolist.node)

node_objects=list()
for(i in 1:length(infolist.node)){
  acc=accessions.node[i]
  info=infolist.node[[acc]]
  info=info[info$COV_HNSCC_tissue_type=="tumor-primary"|info$COV_HNSCC_tissue_type=="tumor-unspecified",]
  node_object=info[,c("SAMPLID","COVAR_N","COVAR_N_status")]
  node_object=node_object[!is.na(node_object$COVAR_N_status) & !is.na(node_object$SAMPLID),]
  rownames(node_object)=node_object$SAMPLID
  node_objects[[acc]]=node_object
}
saveRDS(node_objects, paste(Datadirclin, "Precog.HNSCC.lymph.node.metastasis.objects.all.HNSCC.rds", sep=""))

##################################################################
#Meta analysis of genes associated with lymph node metastasis (LNM) status
#####################################################################

#get gene expression data and clinical data  
alldata=readRDS(paste(Datadir, "Precog.HNSCC.all.exp.clin.data.rds", sep=""))
infolist=alldata$infolist

#list of gene expression datsets
explist=alldata$explist

#
accessions=names(explist)
excluded_studies=c("GSE30788_discovery", "GSE30788_validation")
accessions=setdiff(accessions, excluded_studies)

#Get study IDs
accessions.expnode=intersect(accessions, names(node_objects))
#21 datasets

#restricting further to datasets with at least five LNM+ve and five LNM-ve patients with gene exprssion data 
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
#21 datasets left

################################################
#make matching lists of expression and LNM status data objects
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

#restrict to studies with gene in half or more studies. Otherwise top genes found by random effects models are those that are only in one ot two studies 
allgenes=unique(unlist(lapply(explist.node, rownames)))
#For each gene, find the number of studies the gene is in 
allgeneinstances=unlist(lapply(explist.node, rownames))

#Find the number of studies for each gene is available 
nstudies=lapply(allgenes, function(x) length(grep(paste("\\b",x, "\\b", sep=""), allgeneinstances)))
names(nstudies)=allgenes
saveRDS(nstudies, paste0(Datadir, "nstudies.ge.rds"))

num=round(length(accessions.node)*.5)
genes.available.all=names(unlist(nstudies)[unlist(nstudies)>=num])
#16899 genes in at least half (n=10) studies 

lnm.list=make.meta.node_multistudies(expressionlist=explist.node, studyaccessions=accessions.var, clinlist=varlist_objects, genes=genes.available.all, varfactor="COVAR_N_status")

saveRDS(lnm.list, paste(Resultsdir,"Precog.meta.analysis.node.status.all.HNSCC.updated.03.30.2010.rds", sep=""))

library(devtools)

meta2=lnm.list$meta

metatab=cbind(
unlist(lapply(meta2, function(x) x$pval)),
unlist(lapply(meta2, function(x) x$b)),
unlist(lapply(meta2, function(x) x$zval)))

metatab=as.data.frame(metatab)
dim(metatab[metatab$V3>3.09,])
dim(metatab[metatab$V3<(-3.09),])
colnames(metatab)=c("pval","beta","zval")

#restrict to genes that are signficantly associated with LNM based on z-score threshold
metatab=metatab[abs(metatab$zval)>=3.09,]

#annotate genes
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

saveRDS(metatab, paste0(Resultsdir,"lnm_associated_genes_zscore.rds"))


#####################################################
#Clustering LNM-associated genes based on coexpression using Phenograph
#####################################################

alldata=readRDS("~/Documents/Projects/HNSCC_PreCog/data/Precog.HNSCC.all.exp.clin.data.rds")
explist=alldata$explist
zscores.list=explist
infolist=alldata$infolist

#
allgenes=lapply(explist, rownames)
#for each dataset get the intersection with the signficant gene list 
allgenes2=lapply(allgenes, function(x) intersect(metatab$gene, x))

#Make dataframe indicating the presence or absence of each gene within each study 
genes=unique(unlist(allgenes2))
#872 genes altogether 

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

#Now gene genes that are present in all of the studies that have 80% of genes (i.e., the LNM-associated genes that will be clustered)
gene.present.df.allgenes=gene.present.df[,studiesallgenes]
gene.present.df.allgenes=gene.present.df.allgenes[which(apply(gene.present.df.allgenes, 1, function(x) all(x=="yes"))),]
#743 genes in 20 studies 

explist.sig=explist[which(names(explist) %in% colnames(gene.present.df.allgenes))]
ex.commongenes=lapply(explist.sig, function(x) x[match(genes, rownames(x)),])
ex.commongenes.df=abind(ex.commongenes, along=2)
ex.commongenes.df=na.omit(ex.commongenes.df)
ex.commongenes.df=ex.commongenes.df

#Find LNM gene clusters using Phenograph 
library(Rphenograph)
Rphenograph_out <- Rphenograph(ex.commongenes.df, k = 45)

library(umap)
#make UMAP of LNM-associated genes 
clust.umap = umap(ex.commongenes.df)
saveRDS(clust.umap, paste0(Resultsdir, "UMAP_PRECOG_HNSCC_LNM_gene_clusters.rds"))

membership=as.character(membership(Rphenograph_out[[2]]))
df=data.frame(gene=as.character(rownames(ex.commongenes.df)), cluster=membership)

val2=metatab[match(df$gene, metatab$gene),]$lnm_direction
val2=as.numeric(plyr::revalue(metatab[match(df$gene, metatab$gene),]$lnm_direction, c("Anti-LNM"="1","Pro-LNM"="2")))

metatab$cluster=factor(df[match(metatab$gene, df$gene),"cluster"])
saveRDS(metatab, paste0(Resultsdir,"lnm_associated_genes_zscore.rds"))


#############################################
#########################################
#Meta-analysis to idnetify genes associated with tumor grade
#########################################
############################################

#Run meta-analysis of genes associated with level of differentiation in PRECOG
#have formatted differentiation level
alldata=readRDS("~/Documents/Projects/HNSCC_PreCog/data/Precog_HNSCC_clinical_data/Precog.HNSCC.all.exp.clin.data.rds")
explist=alldata$explist
infolist=alldata$infolist

#get study grade/level of differentiation
varfactor="COV_HNSCC_study_grade"
accessions.grade=names(which(unlist(lapply(infolist, function(x) !all(is.na(x[,varfactor]))))))

infolist.grade=infolist[accessions.grade]
explist.grade=explist[accessions.grade]

genes=rownames(explist.grade[[1]])[1:10]

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
saveRDS(glmlists, paste0(Resultsdir, "glms.all.genes.all.studies.histological.grade.rds"))

#Combine p-values using Liptak's test
allglm.abind=as.data.frame(abind(glmlists, along=1))
rownames(allglm.abind)=NULL
allglm.abind$Z=as.numeric(as.character(allglm.abind$Z))
allglm.abind$N=as.numeric(as.character(allglm.abind$N))

#restricting to genes that found within at least two studies, otherwise its not a meta-analysis
tab=table(allglm.abind$gene)
meta.genes=names(tab[which(tab>=2)])
#25058 genes for which at least two studies have the gene. Running meta analysis on all of these 

#Calculate meta-Z scores
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

threshold="3.09"
metaz.glm.annot$signficant=ifelse(abs(metaz.glm.annot$metaz.glm)>=threshold,"yes","no")
summary(as.factor(metaz.glm.annot$signficant))
#5047/25058 genes

metaz.glm.annot$direction=rep(NA, nrow(metaz.glm.annot))
metaz.glm.annot$direction[metaz.glm.annot$metaz.glm<=-3.09]="Anti-grade"
metaz.glm.annot$direction[metaz.glm.annot$metaz.glm>=3.09]="Pro-grade"

saveRDS(metaz.glm.annot, paste(Resultsdir, "Precog.HNSCC.glm.grade.Liptak.metaz.allHNSCC.rds", sep=""))
#











