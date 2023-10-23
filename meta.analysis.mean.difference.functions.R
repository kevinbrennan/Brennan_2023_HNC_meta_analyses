
#Make function to extract variables needed to meta-analysis, for every gene
library(metafor)
library(esc)

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

#Function that will perform a meta-analyis to get difference in mean expression of a set of genes between two levels of a categorical variable such as HPV pos neg, lymph node status, yes, no, etc.

#expressionlist is a lits of gene expression dataframes 
#clinlist is a list of clinical dataframes that have a column var, which is the variable I want to perform meta-analysis
#function assumes that ncol explist and nrow clinilist are matching for each matched dataframe
#function also assumes that the variable var is a factor with two levels and that the levels are in the same order
#studyaccessions is a character vector of accession or names for each study that links explist to clinlist
#genes is a character vector of genes that are with all studies, so rownames in each study in explist

make.meta=function(clinlist, expressionlist,studyaccessions, genes, varfactor){
  
EffectSizes=list()

for(i in 1:length(studyaccessions)){
  
acc=studyaccessions[i] 

exp=expressionlist[[which(names(expressionlist) %in% acc)]]
info=clinlist[[which(names(clinlist) %in% acc)]]

gt=Get_meta_vars(dataframe=exp[genes,], var=info[,varfactor])
gt2=lapply(c(1:nrow(gt)), function(k) esc_mean_sd(grp1m = gt[k,1], grp1sd = gt[k,2], grp1n = gt[k,3], grp2m = gt[k,4], grp2sd = gt[k,5], grp2n = gt[k,6], es.type = "g", study = paste("Study ",acc, sep="")))
names(gt2)=rownames(gt)

EffectSizes[[i]]=gt2
}
names(EffectSizes)=studyaccessions

#now perform meta analysis for each gene in genes 
meta.genes=list()

for(f in 1:length(genes)){
gene=genes[f]
genestudies=lapply(EffectSizes, function(x) x[[which(names(x) %in% gene)]])
#mydat2=combine_esc(lapply(1:length(genestudies), function(x) c(genestudies[[x]])))
#here is a problem with the code. can't figure out how to pass multiple objects combine_esc as list
mydat2=combine_esc(genestudies[[1]], genestudies[[2]], genestudies[[3]], genestudies[[4]], genestudies[[5]])

#perform meta analysis. random effects model is hard coded in here 

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



#version of function to get meta variables when the input is a vector of genes that are part of a signature 
Get_meta_vars_signature=function(dataframe, var, siggenes){

dataframe=as.data.frame(dataframe)  
var=as.factor(var)
mean1=mean(colMeans(dataframe[siggenes,which(var==levels(var)[1])], na.rm=T), na.rm=T)
sd1=sd(colMeans(dataframe[siggenes,which(var==levels(var)[1])], na.rm=T), na.rm=T)
length1=length(colMeans(dataframe[siggenes,which(var==levels(var)[1])], na.rm=T))
mean2=mean(colMeans(dataframe[siggenes,which(var==levels(var)[2])], na.rm=T), na.rm=T)
sd2=sd(colMeans(dataframe[siggenes,which(var==levels(var)[2])], na.rm=T), na.rm=T)
length2=length(colMeans(dataframe[siggenes,which(var==levels(var)[2])], na.rm=T))

metavars=c(mean1, sd1, length1, mean2, sd2, length2)
names(metavars)=c("mean1","sd1","length1","mean2","sd2","length2")

return(metavars)
}

#version of function to get meta analysis when the input is a vector of genes that are part of a signature 


make.meta.signature=function(clinlist, expressionlist, studyaccessions, siggenes, varfactor){
  
EffectSizes=list()
clinlist=clinlist[studyaccessions]
expressionlist=expressionlist[studyaccessions]

for(i in 1:length(studyaccessions)){
  
acc=studyaccessions[i] 
exp=as.data.frame(expressionlist[[which(names(expressionlist) %in% acc)]])
info=as.data.frame(clinlist[[which(names(clinlist) %in% acc)]])

gt=Get_meta_vars_signature(exp, var=info[,varfactor], siggenes)
gt2=esc_mean_sd(grp1m = gt[1], grp1sd = gt[2], grp1n = gt[3], grp2m = gt[4], grp2sd = gt[5], grp2n = gt[6], es.type = "g", study = paste("Study ",acc, sep=""))

EffectSizes[[i]]=gt2
}
names(EffectSizes)=studyaccessions

mydat2=combine_esc(EffectSizes[[1]], EffectSizes[[2]], EffectSizes[[3]], EffectSizes[[4]], EffectSizes[[5]],
                   EffectSizes[[6]], EffectSizes[[7]], EffectSizes[[8]], EffectSizes[[9]], EffectSizes[[10]],
                   EffectSizes[[11]], EffectSizes[[12]], EffectSizes[[13]], EffectSizes[[14]], EffectSizes[[15]], 
                   EffectSizes[[16]], EffectSizes[[17]], EffectSizes[[18]], EffectSizes[[19]], EffectSizes[[20]])
#Need to add trycatch to metafor to handle convergene error

meta=metafor::rma(yi = es, sei = se, method = "REML", data = mydat2)

return(meta)
}


#Version of TryCatch that logs values errors and warnings
#Means that I can retrieve the warning or error message and go back and deal with the values that throw up and error 
myTryCatch <- function(expr) {
  warn <- err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr, error=function(e) {
      err <<- e
      NULL
      return(NA)
    }), warning=function(w) {
      warn <<- w
      invokeRestart("muffleWarning")
    })
  list(value=value, warning=warn, error=err)
}

#myTryCatch(func(a))

#get=lapply(str, function(y)  myTryCatch(func(as.numeric(as.character(y)))))
#get2=unlist(lapply(get, function(x) x$value))
#errorindex=which(is.na(get2))
#grep("converge", as.character(get[[errorindex]]$warning))

#could use grepl to find indices of genes for which rma could not converge 


###############################
#make meta node 
#################################

make.meta.node=function(clinlist, expressionlist,studyaccessions, genes, varfactor){
  
EffectSizes=list()

for(i in 1:length(studyaccessions)){
  
acc=studyaccessions[i] 

exp=expressionlist[[which(names(expressionlist) %in% acc)]]
info=clinlist[[which(names(clinlist) %in% acc)]]

gt=Get_meta_vars(dataframe=exp[genes,], var=info[,varfactor])
gt2=lapply(c(1:nrow(gt)), function(k) esc_mean_sd(grp1m = gt[k,1], grp1sd = gt[k,2], grp1n = gt[k,3], grp2m = gt[k,4], grp2sd = gt[k,5], grp2n = gt[k,6], es.type = "g", study = paste("Study ",acc, sep="")))
names(gt2)=rownames(gt)

EffectSizes[[i]]=gt2
}
names(EffectSizes)=studyaccessions

#now perform meta analysis for each gene in genes 
meta.genes=list()

for(f in 1:length(genes)){
gene=genes[f]
genestudies=lapply(EffectSizes, function(x) x[[which(names(x) %in% gene)]])
#mydat2=combine_esc(lapply(1:length(genestudies), function(x) c(genestudies[[x]])))
#here is a problem with the code. can't figure out how to pass multiple objects combine_esc as list
mydat2=combine_esc(genestudies[[1]], genestudies[[2]], genestudies[[3]], genestudies[[4]], genestudies[[5]], genestudies[[6]], genestudies[[7]], genestudies[[8]], genestudies[[9]],
                   genestudies[[10]], genestudies[[11]], genestudies[[12]], genestudies[[13]], genestudies[[14]], genestudies[[15]], genestudies[[16]], genestudies[[17]], genestudies[[18]])

#perform meta analysis. random effects model is hard coded in here 

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







make.meta.node_multistudies=function(clinlist, expressionlist, studyaccessions, genes, varfactor){
  
EffectSizes=list()

for(i in 1:length(studyaccessions)){
  
acc=studyaccessions[i] 

exp=expressionlist[[acc]]
info=clinlist[[acc]]

gt=Get_meta_vars(dataframe=exp, var=info[,varfactor])
gt2=lapply(c(1:nrow(gt)), function(k) esc_mean_sd(grp1m = gt[k,1], grp1sd = gt[k,2], grp1n = gt[k,3], grp2m = gt[k,4], grp2sd = gt[k,5], grp2n = gt[k,6], es.type = "g", study = paste("Study ",acc, sep="")))
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

#Ridiculous wrapper needed for combine_esc for between 5 and 21 studies. Need to extend to include numbers of studies outside this range
combine_esc_KB=function(g){
  if(length(g)==2){
    out=combine_esc(g[[1]], g[[2]])
  } else {
    if(length(g)==3){
      out=combine_esc(g[[1]], g[[2]], g[[3]])
    } else {
      if(length(g)==4){
        out=combine_esc(g[[1]], g[[2]], g[[3]], g[[4]])
      } else {
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
      }
    }
  }
  return(out)  
}


#Version that that runs meta-analysis but just increases the number of iterations in the rma algorithm to make sure it converges
make.meta.node_multistudies_maxiter1000=function(clinlist, expressionlist, studyaccessions, genes, varfactor){
  
EffectSizes=list()

for(i in 1:length(studyaccessions)){
  
acc=studyaccessions[i] 

exp=expressionlist[[which(names(expressionlist) %in% acc)]]
info=clinlist[[which(names(clinlist) %in% acc)]]

gt=Get_meta_vars(dataframe=exp, var=info[,varfactor])
gt2=lapply(c(1:nrow(gt)), function(k) esc_mean_sd(grp1m = gt[k,1], grp1sd = gt[k,2], grp1n = gt[k,3], grp2m = gt[k,4], grp2sd = gt[k,5], grp2n = gt[k,6], es.type = "g", study = paste("Study ",acc, sep="")))
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
rma.res=myTryCatch(metafor::rma(yi = es, sei = se, method = "REML", data = mydat2, control=list(maxiter=1000)))

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







######################################################
#functions for CIBERSORT meta-analysis with binary variables such as HPV status, sex, and lymph node status 
####################################################

make.meta.binary.var.cibersort=function(clinlist, ciblist, studyaccessions, cells, varfactor){
  
EffectSizes=list()

for(i in 1:length(studyaccessions)){
  
acc=studyaccessions[i] 

exp=ciblist[[which(names(ciblist) %in% acc)]]
info=clinlist[[which(names(clinlist) %in% acc)]]

gt=Get_meta_vars_cib(dataframe=exp, var=info[,varfactor])

gt2=lapply(c(1:nrow(gt)), function(k) esc_mean_sd(grp1m = gt[k,1], grp1sd = gt[k,2], grp1n = gt[k,3], grp2m = gt[k,4], grp2sd = gt[k,5], grp2n = gt[k,6], es.type = "g", study = paste("Study ",acc, sep="")))
names(gt2)=rownames(gt)

EffectSizes[[i]]=gt2
}
names(EffectSizes)=studyaccessions

#now perform meta analysis for each gene in genes 
meta.cells=list()

for(f in 1:length(cellTypes)){
cellType=cellTypes[f]
#get the indices of studies in which the gene exists
indices=c(which(unlist(lapply(EffectSizes, function(x) cellType %in% names(x)))))
#get the effect sizes for studies that have this gene 

genestudies=lapply(EffectSizes[c(indices)], function(x) x[[which(names(x) %in% cellType)]])

#Just need to hack this so it will run for any number of studies 
mydat2=combine_esc_KB(genestudies)

#mydat2=combine_esc(lapply(1:length(genestudies), function(x) c(genestudies[[x]])))
#here is a problem with the code. can't figure out how to pass multiple objects combine_esc as list

rma.res=myTryCatch(metafor::rma(yi = es, sei = se, method = "REML", data = mydat2))

if(is.null(rma.res$error)){
  meta.cells[[f]]=rma.res$value
 } else {
  meta.cells[[f]]=rma.res$error
}

}
names(meta.cells)=cellTypes

EffectSizesList=list(meta.cells, EffectSizes)
names(EffectSizesList)=c("meta","effect.sizes")

return(EffectSizesList)
}


################
library(metafor)
library(esc)

Get_meta_vars_cib=function(dataframe=dataframe, var=as.factor(var)){
        dat=array(NA,c(length(cellTypes),6))
        rownames(dat)=cellTypes
        colnames(dat)=c("grp1m", "grp1sd", "grp1n", "grp2m", "grp2sd", "grp2n")
        var=as.factor(var)
          
        for( i in 1:length(cellTypes)){
        cellType=cellTypes[i]
        dat[i,1]=mean(as.numeric(as.character(dataframe[which(var==levels(var)[1]),cellType])), na.rm=T)
        dat[i,2]=sd(as.numeric(as.character(dataframe[which(var==levels(var)[1]),cellType])), na.rm=T)
        dat[i,3]=length(which(!is.na(as.numeric(as.character(dataframe[which(var==levels(var)[1]),cellType])))))
        dat[i,4]=mean(as.numeric(as.character(dataframe[which(var==levels(var)[2]),cellType])), na.rm=T)
        dat[i,5]=sd(as.numeric(as.character(dataframe[which(var==levels(var)[2]),cellType])), na.rm=T)
        dat[i,6]=length(which(!is.na(as.numeric(as.character(dataframe[which(var==levels(var)[2]),cellType])))))
        } 
            
    dat=as.data.frame(dat, drop=FALSE)
    return(dat)
}



        















