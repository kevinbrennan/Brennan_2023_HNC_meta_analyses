
##############################################################
##############################################################
#Processing datasets for meta-analysis 
##############################################################
##############################################################

RscriptsPath="" #Directory for auxiliart R scripts
source(paste(RscriptsPath, "ProcessRawAffyForPrecog.R", sep=""))
expdir="" #Directory that stores raw expression data foles for all studies
celdir.output=paste(expdir, "zscores_for_genes_expression_analysis_affy/", sep="")
dir.create(paste(celdir.output,"log2/", sep=""))
#dir.create(celdir.output)

#Loding the list of clinical datasets, which were manually curated 
infolist=readRDS(paste(Datadirclin, "Precog.HNSCC.clinical.data.list.all.samples.rds"))

##############################################################
#Auxiliary scripts for processing datasets
##############################################################

library(preprocessCore)
library(WGCNA)
library(affy)
library(GEOquery)

#load packages needed to process various affy arrays from brainarray
library(hgu133plus2hsentrezg.db)
library(hgu133plus2hsentrezgcdf)
library(hgu133plus2hsentrezgprobe)

library(hgu133ahsentrezg.db)
library(hgu133ahsentrezgcdf)
library(hgu133ahsentrezgprobe)

library(hgu95ahsentrezg.db)
library(hgu95ahsentrezgcdf)
library(hgu95ahsentrezgprobe)

library(hgu95av2hsentrezg.db)
library(hgu95av2hsentrezgcdf)
library(hgu95av2hsentrezgprobe)

library(hgu133ahsentrezg.db)
library(hgu133ahsentrezgcdf)
library(hgu133ahsentrezgprobe)

library(hgu133a2hsentrezg.db)
library(hgu133a2hsentrezgcdf)
library(hgu133a2hsentrezgcdf)

library(hgu133plus2hsentrezg.db)
library(hgu133plus2hsentrezgcdf)
library(hgu133plus2hsentrezgprobe)

library(hta20hsentrezg.db)
library(hta20hsentrezgcdf)
library(hta20hsentrezgprobe)

#Function to process raw affy CEL files for affy array using supplied brainarray cdfs. Output is gene expression data in liear space. Not Log2 trasnfromed or z-scored
ProcessRawAffyForPrecog=function(gse, cdf, entrezIDmap, celdir.output, clininfo, celdirectory){

data=ReadAffy(celfile.path=celdirectory, cdfname=cdf) 
eset=affy::mas5(data)
#apply mas5
mas5.ALL=exprs(eset)
colnames(mas5.ALL)=gsub(" ","", colnames(mas5.ALL))
colnames(mas5.ALL)=gsub("[.].*","", colnames(mas5.ALL))
colnames(mas5.ALL)=gsub("_.*","", colnames(mas5.ALL))

#Remove Affy control probes, which are at the end of the document and mess with the adding gene symbols to the rows
mas5.ALL=mas5.ALL[-grep("AFFX-", rownames(mas5.ALL)),]
#Format values to 5 decimal places
mas5.ALL=format(mas5.ALL, digits=5)
mas5.ALL=as.matrix(mas5.ALL)

mas5.ALL2=apply(mas5.ALL, 1, function(x) as.numeric(as.character(x)))
mas5.ALL2=as.data.frame(t(mas5.ALL2))
colnames(mas5.ALL2)=colnames(mas5.ALL)

#need to restrict to tumor samples before quantile normalizing 
#get processed clinical data associated datset to restrict samples to primary tumors before quantile normalization 
#files.gse=names(exp.data.list)
info=clininfo[[which(names(clininfo)==gse)]]
IDs=info[which(info$COV_HNSCC_tissue_type=="tumor-unspecified"|info$COV_HNSCC_tissue_type=="tumor-primary"),"SAMPLID"]
mas5.ALL2=mas5.ALL2[,colnames(mas5.ALL2) %in% IDs]

#quantile normalize
mas5.ALL.qn=normalize.quantiles(as.matrix(mas5.ALL2))
dimnames(mas5.ALL.qn)=dimnames(mas5.ALL2)

#collapse probes to genes 
probes.ALL=row.names(mas5.ALL)
entrezID.ALL = unlist(mget(probes.ALL, entrezIDmap))
ID.ALL = unlist(mget(probes.ALL, entrezIDmap))
collapserows.annot=as.data.frame(cbind(ID.ALL, entrezID.ALL))
collapserows.annot$probe=rownames(collapserows.annot)
collapserows.annot=na.omit(collapserows.annot[,c("entrezID.ALL", "probe")])

mas5.collapse=collapseRows(datET=mas5.ALL.qn, rowGroup=collapserows.annot$entrezID.ALL, rowID=collapserows.annot$probe)
mas5.c=mas5.collapse$datETcollapsed

mas5.clog2=log2(mas5.c)

#Now z-score the genes
mas5.z=apply(mas5.clog2, 1, function(x)  (x-mean(x))/sd(x))
mas5.z=as.data.frame(t(mas5.z))

saveRDS(mas5.z, paste(celdir.output, "processed_exp_z_scores_", gse, ".rds", sep=""))
saveRDS(mas5.clog2, paste(celdir.output,"log2/", "processed_exp_log2_", gse, ".rds", sep=""))

return(mas5.z)
}

#Version of ProcessRawAffyForPrecog from the quantile normalization step onward
#For datasets that require manual processing of the first steps
ProcessRawAffyForPrecogFromQuantialNormalization=function(mas5.processed.df, gse, entrezIDmap, celdir.output){

#quantile normalize
mas5.ALL.qn=normalize.quantiles(as.matrix(mas5.processed.df))
dimnames(mas5.ALL.qn)=dimnames(mas5.processed.df)

#collapse probes to genes 
probes.ALL=row.names(mas5.ALL.qn)
entrezID.ALL = unlist(mget(probes.ALL, entrezIDmap))
ID.ALL = unlist(mget(probes.ALL, entrezIDmap))
collapserows.annot=as.data.frame(cbind(ID.ALL, entrezID.ALL))
collapserows.annot$probe=rownames(collapserows.annot)
collapserows.annot=na.omit(collapserows.annot[,c("entrezID.ALL", "probe")])

mas5.collapse=collapseRows(datET=mas5.ALL.qn, rowGroup=collapserows.annot$entrezID.ALL, rowID=collapserows.annot$probe)
mas5.c=mas5.collapse$datETcollapsed

mas5.clog2=log2(mas5.c)

#Now z-score the genes
mas5.z=apply(mas5.clog2, 1, function(x)  (x-mean(x))/sd(x))
mas5.z=as.data.frame(t(mas5.z))

saveRDS(mas5.z, paste(celdir.output, "processed_exp_z_scores_", gse, ".rds", sep=""))
saveRDS(mas5.clog2, paste(celdir.output,"log2/", "processed_exp_log2_", gse, ".rds", sep=""))

return(mas5.z)
}

########################################################
#processing publicly accessible gene expression datasets that were analyzed within the meta-analyses
########################################################

#process all Affymetrix hgu133a datasets from GEO
cdf="hgu133ahsentrezgcdf"
entrezIDmap=hgu133ahsentrezgENTREZID

GPL96_accessions=c("GSE3524", "GSE2280", "GSE27020")

for(i in 1:length(GPL96_accessions)){
  gse=GPL96_accessions[i]
  celdir=paste(expdir, gse, "/", gse, "_RAW", "/", sep="")
  assign(paste(gse, "_processed", sep=""), ProcessRawAffyForPrecog(gse, cdf, entrezIDmap, celdir.output, clininfo=infolist, celdirectory = celdir))
}

#process Affymetrix hgu133plus2 datasets from GEO

cdf="hgu133plus2hsentrezgcdf"
entrezIDmap=hgu133plus2hsentrezgENTREZID

GPL570_accessions=c("GSE9844", "GSE42743", "GSE30784", "GSE41613", "GSE78060", "GSE31056", "GSE3292", "GSE6791")

for(i in 1:length(GPL570_accessions)){
  gse=GPL570_accessions[i]
  celdir=paste(expdir, gse, "/", gse, "_RAW", "/", sep="")
  assign(paste(gse, "_processed", sep=""), ProcessRawAffyForPrecog(gse, cdf, entrezIDmap, celdir.output, clininfo=infolist, celdirectory = celdir))
}

#Process datasets that were generated using other Affymetrix arrays

#GSE6631 Kuriakose
gse="GSE6631"
celdir=paste(expdir, gse, "/", gse, "_RAW", "/", sep="")
cdf="hgu95av2hsentrezgcdf"
entrezIDmap=hgu95av2hsentrezgENTREZID
celdirectory=celdir

assign(paste(gse, "_processed", sep=""), ProcessRawAffyForPrecog(gse, cdf, entrezIDmap, celdir.output, clininfo=infolist, celdirectory = celdir))


#Jung E-MTAB-1328
gse="E-MTAB-1328"
celdir=paste(expdir, gse, "/", sep="")
cdf="hgu133plus2hsentrezgcdf"
entrezIDmap=hgu133plus2hsentrezgENTREZID

data=ReadAffy(celfile.path=celdir, cdfname=cdf) 
eset=affy::mas5(data)
#apply mas5
mas5.ALL=exprs(eset)
colnames(mas5.ALL)=gsub("_U133.*","", colnames(mas5.ALL))

mas5.ALL=mas5.ALL[-grep("AFFX-", rownames(mas5.ALL)),]
mas5.ALL=format(mas5.ALL, digits=5)
mas5.ALL=as.matrix(mas5.ALL)

mas5.ALL2=apply(mas5.ALL, 1, function(x) as.numeric(as.character(x)))
mas5.ALL2=as.data.frame(t(mas5.ALL2))
colnames(mas5.ALL2)=colnames(mas5.ALL)

#problem with trailing zeros
IDs %in% colnames(mas5.ALL)
gsub("_0","_", IDs) %in% colnames(mas5.ALL)

#need to restrict to tumor samples before quantile normalizing 
#get processed clinical data associated datset to restrict samples to primary tumors before quantile normalization 
#files.gse=names(exp.data.list)
clininfo=infolist
info=clininfo[[which(names(clininfo)==gse)]]
IDs=info[which(info$COV_HNSCC_tissue_type=="tumor-unspecified"|info$COV_HNSCC_tissue_type=="tumor-primary"),"SAMPLID"]
mas5.ALL2=mas5.ALL2[,colnames(mas5.ALL2) %in% IDs]

assign(paste(gse, "_processed", sep=""), ProcessRawAffyForPrecogFromQuantialNormalization(mas5.ALL2, gse, entrezIDmap, celdir.output))


#GSE23036 Pavon
gse="GSE23036"
celdir=paste(expdir, gse, "/", gse, "_RAW", "/", sep="")

cdf="hgu133a2hsentrezgcdf"
entrezIDmap=hgu133a2hsentrezgENTREZID

assign(paste(gse, "_processed", sep=""), ProcessRawAffyForPrecog(gse, cdf, entrezIDmap, celdir.output, clininfo=infolist, celdirectory = celdir))


############################################
####################################
#Processing datasets that needed to be curated manually 
############################################
############################################

gse="GSE2379"

#Find out which samples belong to which array
getSet=getGEO(gse)
dat2=exprs(getSet$`GSE2379-GPL8300_series_matrix.txt.gz`)
files_GPL8300=colnames(dat2)

cdf="hgu95av2hsentrezgcdf"
entrezIDmap=hgu95av2hsentrezgENTREZID

celdir=paste(expdir, gse, "/", gse, "_RAW", "/", sep="")
files=list.files(celdir, pattern=".CEL")
filenames=gsub("[.].*","", files)

CELfiles_GPL8300=files[filenames %in% files_GPL8300]

celdirectory=celdir
#try instead of changing directory
setwd(celdir)
data=ReadAffy(filenames=CELfiles_GPL8300, cdfname=cdf) 
eset=affy::mas5(data)
#apply mas5
mas5.ALL=exprs(eset)
colnames(mas5.ALL)=gsub(" ","", colnames(mas5.ALL))
colnames(mas5.ALL)=gsub("[.].*","", colnames(mas5.ALL))
colnames(mas5.ALL)=gsub("_.*","", colnames(mas5.ALL))

#Remove Affy control probes, which are at the end of the document and mess with the adding gene symbols to the rows
mas5.ALL=mas5.ALL[-grep("AFFX-", rownames(mas5.ALL)),]
#Format values to 5 decimal places
mas5.ALL=format(mas5.ALL, digits=5)
mas5.ALL=as.matrix(mas5.ALL)

mas5.ALL2=apply(mas5.ALL, 1, function(x) as.numeric(as.character(x)))
mas5.ALL2=as.data.frame(t(mas5.ALL2))
colnames(mas5.ALL2)=colnames(mas5.ALL)

#need to restrict to tumor samples before quantile normalizing 
#get processed clinical data associated datset to restrict samples to primary tumors before quantile normalization 
#files.gse=names(exp.data.list)
clininfo=infolist


gse="GSE2379_GPL8300"
info=clininfo[[which(names(clininfo)==gse)]]
IDs=info[which(info$COV_HNSCC_tissue_type=="tumor-unspecified"|info$COV_HNSCC_tissue_type=="tumor-primary"),"SAMPLID"]
mas5.ALL2=mas5.ALL2[,colnames(mas5.ALL2) %in% IDs]

assign(paste(gse, "_processed", sep=""), ProcessRawAffyForPrecogFromQuantialNormalization(mas5.ALL2, gse, entrezIDmap, celdir.output))

cdf="hgu95ahsentrezgcdf"
symbolmap=hgu95ahsentrezgSYMBOL
entrezIDmap=hgu95ahsentrezgENTREZID

dat1=exprs(getSet$`GSE2379-GPL91_series_matrix.txt.gz`)
files_GPL91=colnames(dat1)

CELfiles_GPL91=files[filenames %in% files_GPL91]
CELfiles_GPL91.full=paste(celdir, CELfiles_GPL91, sep="")

data=ReadAffy(filenames=CELfiles_GPL91, cdfname=cdf) 

eset=affy::mas5(data)
#apply mas5
mas5.ALL=exprs(eset)
colnames(mas5.ALL)=gsub("[.].*","", colnames(mas5.ALL))
#Remove Affy control probes, which are at the end of the document and mess with the adding gene symbols to the rows
mas5.ALL=mas5.ALL[-grep("AFFX-", rownames(mas5.ALL)),]
#Format values to 5 decimal places
mas5.ALL=format(mas5.ALL, digits=5)
mas5.ALL=as.matrix(mas5.ALL)

mas5.ALL2=apply(mas5.ALL, 1, function(x) as.numeric(as.character(x)))
mas5.ALL2=as.data.frame(t(mas5.ALL2))
colnames(mas5.ALL2)=colnames(mas5.ALL)

clininfo=infolist

gse="GSE2379_GPL91"
info=clininfo[[which(names(clininfo)==gse)]]
IDs=info[which(info$COV_HNSCC_tissue_type=="tumor-unspecified"|info$COV_HNSCC_tissue_type=="tumor-primary"),"SAMPLID"]
mas5.ALL2=mas5.ALL2[,colnames(mas5.ALL2) %in% IDs]

assign(paste(gse, "_processed", sep=""), ProcessRawAffyForPrecogFromQuantialNormalization(mas5.ALL2, gse, entrezIDmap, celdir.output))

######################################################
Datasets that require oligo package# 
#########################################################

gse="GSE95805"
celdir=paste(expdir, gse, "/", gse, "_RAW", "/", sep="")

library(pd.hta.2.0)
library(oligo)
celFiles <- list.celfiles(celdir, full.names=TRUE)
dat=oligo::read.celfiles(celFiles)

eset=oligo::rma(dat, background=TRUE, normalize=TRUE, subset=NULL, target="core")

#does biomaRt have 
exp=exprs(eset)
pdata <- pData(eset) # data.frame of phenotypic information.

colnames(exp)=gsub(" ","", colnames(exp))
colnames(exp)=gsub("[.].*","", colnames(exp))
colnames(exp)=gsub("_.*","", colnames(exp))

#quantile normalization step
info=clininfo[[which(names(clininfo)==gse)]]
IDs=info[which(info$COV_HNSCC_tissue_type=="tumor-unspecified"|info$COV_HNSCC_tissue_type=="tumor-primary"),"SAMPLID"]
exp=exp[,colnames(exp) %in% IDs]

#quantile normalize
exp.qn=normalize.quantiles(as.matrix(exp))
dimnames(exp.qn)=dimnames(exp)

library(affycoretools)
#eset <- annotateEset(eset, pd.hta.2.0)
eset2=annotateEset(eset, pd.hta.2.0, columns = c("PROBEID", "ENTREZID", "SYMBOL", "GENENAME"), multivals = "first")

#annotate
fdata=fData(eset2)
fdata=fdata[!is.na(fdata$ID),]

library(biomaRt)
ensembl=useMart("ensembl", dataset="hsapiens_gene_ensembl")

#trying to cnver to gene names uing pd.hta.2.0 then biomaRT
refseqgenes=getBM(attributes = c("refseq_mrna","entrezgene_id","hgnc_symbol"), filters = "refseq_mrna", values = fdata$ID, mart = ensembl)
refseqgenes=refseqgenes[match(fdata$ID, refseqgenes$refseq_mrna),]
fdata=cbind(fdata, refseqgenes)

collapserows.annot=as.data.frame(fdata[,c("PROBEID","entrezgene_id")])
colnames(collapserows.annot)=c("probe","entrezgene_id")

collapserows.annot=na.omit(collapserows.annot)
collapserows.annot=collapserows.annot[nchar(collapserows.annot$entrezgene_id)>0 & nchar(collapserows.annot$probe)>0,]

exp2=exp.qn[rownames(exp.qn) %in% collapserows.annot$probe,]

library(WGCNA)
collapse=collapseRows(datET=exp2, rowGroup=collapserows.annot$entrezgene_id, rowID=collapserows.annot$probe)
dat3=collapse$datETcollapsed
colnames(dat3)=gsub("_.*","", colnames(dat3))

dat3.clog2=log2(dat3)
saveRDS(dat3.clog2, paste(celdir.output,"log2/", "processed_exp_log2_", gse, ".rds", sep=""))

#Now z-score the genes
dat3.z=apply(dat3.clog2, 1, function(x)  (x-mean(x))/sd(x))
dat3.z=as.data.frame(t(dat3.z))

saveRDS(dat3.z, paste(celdir.output, "processed_exp_z_scores_", gse, ".rds", sep=""))

#Thurlow
gse="Thurlow"

RscriptsPath=dirs[2]
source(paste(RscriptsPath, "Process_data_filter_missing_KB.R", sep=""))

RscriptsPath=dirs[2]
source(paste(RscriptsPath, "ProcessRawAffyLinearSpaceDeconvolution.R", sep=""))

library(hgu133plus2hsentrezg.db)
library(hgu133plus2hsentrezgcdf)
library(hgu133plus2hsentrezgprobe)
cdfprobes=hgu133plus2hsentrezgprobe

cdf="hgu133plus2hsentrezgcdf"
symbolmap=hgu133plus2hsentrezgSYMBOL
entrezIDmap=hgu133plus2hsentrezgENTREZID

#Accessed normalized counts dataset from authors of paper 
expdata=readRDS(paste(Datadirexp, "Head.and.neck.cancer.Primary.Thurlow_HNSCC.HGU133Plus2_EntrezCDF.data.RData", sep=""))

rownames(expdata)=expdata$Name
expdata=expdata[,3:ncol(expdata)]
expdata=expdata[-grep("AFFX-", rownames(expdata)),]
expdata=Process_data_filter_missing_KB(expdata)

#quantile normalization step
info=clininfo[[which(names(clininfo)==gse)]]
IDs=info[which(info$COV_HNSCC_tissue_type=="tumor-unspecified"|info$COV_HNSCC_tissue_type=="tumor-primary"),"SAMPLID"]
exp=expdata[,colnames(expdata) %in% IDs]

#quantile normalize
exp.qn=normalize.quantiles(as.matrix(exp))
dimnames(exp.qn)=dimnames(exp)
expdata=exp.qn

#collapse probes to genes 
probes.ALL=row.names(expdata)
probes.ALL=probes.ALL[probes.ALL %in% cdfprobes$Probe.Set.Name]
entrez.ALL=unlist(mget(probes.ALL, entrezIDmap))

collapserows.annot=as.data.frame(as.matrix(entrez.ALL), drop=FALSE)
collapserows.annot$probe=rownames(collapserows.annot)
colnames(collapserows.annot)=c("entrezgene_id","probe")
collapserows.annot=na.omit(collapserows.annot[,c("entrezgene_id", "probe")])

expdata.qn.collapse=collapseRows(datET=expdata, rowGroup=collapserows.annot$entrezgene_id, rowID=collapserows.annot$probe)
expdata.c=expdata.qn.collapse$datETcollapsed

dat3.clog2=log2(expdata.c)
saveRDS(dat3.clog2, paste(celdir.output,"log2/", "processed_exp_log2_", gse, ".rds", sep=""))

#Now z-score the genes
dat3.z=apply(dat3.clog2, 1, function(x)  (x-mean(x))/sd(x))
dat3.z=as.data.frame(t(dat3.z))

saveRDS(dat3.z, paste(celdir.output, "processed_exp_z_scores_", gse, ".rds", sep=""))

#GSE41116 Pickering
library(R.utils)
library(pd.hta.2.0)
library(oligo)
library(affycoretools)
library(huex10sttranscriptcluster.db)
library(affycore)

#Affymetrix Human Exon 1.0 ST Array
gse="GSE41116"
celdir=paste(expdir, gse, "/", gse, "_RAW", "/", sep="")

cels=list.files(celdir, pattern=".CEL.gz")
setwd(celdir)
lapply(cels, function(x) gunzip(x, remove=FALSE))
lapply(list.files(celdir, pattern=".gz"), file.remove)

celFiles <- list.celfiles(celdir, full.names=TRUE)
dat=oligo::read.celfiles(celFiles)

eset=oligo::rma(dat, background=TRUE, normalize=TRUE, subset=NULL, target="core")

eset2=annotateEset(eset, huex10sttranscriptcluster.db, columns = c("PROBEID", "ENTREZID", "SYMBOL", "GENENAME"), multivals = "first")
fdata=fData(eset2)

exp=exprs(eset2)
colnames(exp)=gsub("_.*","",colnames(exp))

exp=format(exp, digits=5)
exp=as.matrix(exp)

exp2=apply(exp, 1, function(x) as.numeric(as.character(x)))
exp2=as.data.frame(t(exp2))
colnames(exp2)=colnames(exp)
exp=exp2

#restrict to tumors
info=clininfo[[which(names(clininfo)==gse)]]
IDs=info[which(info$COV_HNSCC_tissue_type=="tumor-unspecified"|info$COV_HNSCC_tissue_type=="tumor-primary"),"SAMPLID"]
exp=exp[,colnames(exp) %in% IDs]

#quantile normalize
exp.qn=normalize.quantiles(as.matrix(exp))
dimnames(exp.qn)=dimnames(exp)

collapserows.annot=as.data.frame(fdata[,c("ENTREZID","PROBEID")])
colnames(collapserows.annot)=c("entrezgene_id","probe")
collapserows.annot=na.omit(collapserows.annot[,c("entrezgene_id", "probe")])

library(WGCNA)
exp.collapse=collapseRows(datET=exp.qn, rowGroup=collapserows.annot$entrezgene_id, rowID=collapserows.annot$probe)
expdata.c=exp.collapse$datETcollapsed
#17284 genes

dat3.clog2=log2(expdata.c)
saveRDS(dat3.clog2, paste(celdir.output,"log2/", "processed_exp_log2_", gse, ".rds", sep=""))

#Now z-score the genes
dat3.z=apply(dat3.clog2, 1, function(x)  (x-mean(x))/sd(x))
dat3.z=as.data.frame(t(dat3.z))

saveRDS(dat3.z, paste(celdir.output, "processed_exp_z_scores_", gse, ".rds", sep=""))

#GSE33205 Fertig
gse="GSE33205"
celdir=paste(expdir, gse, "/", gse, "_RAW", "/", sep="")

#When gunzip works, remove zipped files
lapply(list.files(celdir, pattern=".gz"), file.remove)

celFiles <- list.celfiles(celdir, full.names=TRUE)
dat=oligo::read.celfiles(celFiles)

eset=oligo::rma(dat, background=TRUE, normalize=TRUE, subset=NULL, target="core")

#eset <- annotateEset(eset, pd.hta.2.0)
eset2=annotateEset(eset, huex10sttranscriptcluster.db, columns = c("PROBEID", "ENTREZID", "SYMBOL", "GENENAME"), multivals = "first")
fdata=fData(eset2)

exp=exprs(eset2)
colnames(exp)=gsub("[.].*","",colnames(exp))

exp=format(exp, digits=5)
exp=as.matrix(exp)

exp2=apply(exp, 1, function(x) as.numeric(as.character(x)))
exp2=as.data.frame(t(exp2))
colnames(exp2)=colnames(exp)
exp=exp2

#restrict to tumors
info=clininfo[[which(names(clininfo)==gse)]]
IDs=info[which(info$COV_HNSCC_tissue_type=="tumor-unspecified"|info$COV_HNSCC_tissue_type=="tumor-primary"),"SAMPLID"]
exp=exp[,colnames(exp) %in% IDs]

#quantile normalize
exp.qn=normalize.quantiles(as.matrix(exp))
dimnames(exp.qn)=dimnames(exp)

collapserows.annot=as.data.frame(fdata[,c("ENTREZID","PROBEID")])
colnames(collapserows.annot)=c("entrezgene_id","probe")
collapserows.annot=na.omit(collapserows.annot[,c("entrezgene_id", "probe")])

library(WGCNA)
exp.collapse=collapseRows(datET=exp.qn, rowGroup=collapserows.annot$entrezgene_id, rowID=collapserows.annot$probe)
expdata.c=exp.collapse$datETcollapsed
#17284 genes

dat3.clog2=log2(expdata.c)
saveRDS(dat3.clog2, paste(celdir.output,"log2/", "processed_exp_log2_", gse, ".rds", sep=""))

#Now z-score the genes
dat3.z=apply(dat3.clog2, 1, function(x)  (x-mean(x))/sd(x))
dat3.z=as.data.frame(t(dat3.z))

saveRDS(dat3.z, paste(celdir.output, "processed_exp_z_scores_", gse, ".rds", sep=""))

######################################################
###################################################
#make single large list of gene expression datasets with matched clinical datasets, to which meta-analyses can easily be applied 
###################################################
#######################################################

celdir.output=paste(expdir, "zscores_for_genes_expression_analysis_affy/", sep="")

affy.processed.files=list.files(celdir.output)
affy.processed.files=setdiff(affy.processed.files, "log2")
accessions.affy.processed=gsub("processed_exp_z_scores_","", gsub(".rds","", affy.processed.files))
#accessions.affy.processed=setdiff(accessions.affy.processed, "log2")

setwd(celdir.output)
exp.list.affy=lapply(affy.processed.files, readRDS)
names(exp.list.affy)=accessions.affy.processed

#Get non-affy datasets 
expdir.nonaffy=paste(expdir, "zscores_for_genes_expression_analysis_non_affy/", sep="")

nonaffy.processed.files=list.files(expdir.nonaffy)
nonaffy.processed.files=setdiff(nonaffy.processed.files, "log2")
accessions.nonaffy.processed=gsub("processed_exp_z_scores_","", gsub(".rds","", nonaffy.processed.files))

setwd(expdir.nonaffy)
exp.list.nonaffy=lapply(nonaffy.processed.files, readRDS)
names(exp.list.nonaffy)=accessions.nonaffy.processed

exp.list.combined=c(exp.list.affy, exp.list.nonaffy)
#30 datasets 

#add Previously processed TCGA data 
expression.all.datasets.zscores.updated=readRDS(paste(Datadirexp, "expression.all.datasets.zscores.updated.rds", sep=""))
exp.list.combined[["TCGA"]]=expression.all.datasets.zscores.updated$TCGA
#Also missing GSE23558, which is the last datset I proessed
exp.list.combined[["GSE23558"]]=expression.all.datasets.zscores.updated$GSE23558

zscores.list=exp.list.combined
clininfolist=readRDS(paste(Datadirclin, "Precog.HNSCC.clinical.data.list.all.samples.rds", sep=""))

accessions=names(clininfolist)
infolist=list()
explist=list()

for(i in 1:length(accessions)){
  acc=accessions[i]
  info=clininfolist[[which(names(clininfolist) %in% acc)]]
  
  info.tumor.IDs=info[info$COV_HNSCC_tissue_type=="tumor-primary"|info$COV_HNSCC_tissue_type=="tumor-unspecified", "SAMPLID"]
  exp=zscores.list[[which(names(zscores.list) %in% acc)]]
  
  OverlapSamples=intersect(info.tumor.IDs, colnames(exp))
  exp=exp[,OverlapSamples]
  info=info[info$SAMPLID %in% OverlapSamples,]
  info=info[match(OverlapSamples, info$SAMPLID),]
  rownames(info)=info$SAMPLID
  
  infolist[[i]]=info
  explist[[i]]=exp
}

names(infolist)=accessions
names(explist)=accessions

alldata=list(explist, infolist)
names(alldata)=c("explist","infolist")

saveRDS(alldata, paste0(Datadirclin, "Precog.HNSCC.all.exp.clin.data.rds")


###########################################################################
###########################################################################
#Processing the Stanford bulk RNA-Seq dataset
###########################################################################
###########################################################################

###########################################
#Run trimgalore! to trim and filter reads
###########################################

cd dir/python-virtual-environments
source dir/env/bin/activate
trim_galore=dir/TrimGalore-0.6.0/trim_galore
export PATH=dir/TrimGalore-0.6.0/:$PATH
pigz=dir/NGStools/cell_ranger/cellranger-3.1.0/miniconda-cr-cs/4.3.21-miniconda-cr-cs-c10/bin/pigz
parallel=/bin/parallel

#Setup needed to run trim_galore

samplesDir=dir/samples/FASTQ
cd $samplesDir

ls -1 *_R1.fastq.gz | cut -d_ -f1 | sort | uniq | parallel -j 4 'trim_galore -q 20 --stringency 3 --gzip --length 20  --paired {}_R1.fastq.gz {}_R2.fastq.gz --output_dir dir/samples/trim_galore_output'

################################
#Running kallisto to pseudoalign reads
##################################

toolsDir=dir/NGStools
kallisto=dir/kallisto_linux-v0.46.0/kallisto

refDir=dir/RNASeq/
KALLISTO_INDEX=/kallisto.gencodev34.idx #Point to transcriptome index file 

IDsuffix1="_R1_val_1"
IDsuffix2="_R2_val_2"

samplesDir=dir/samples/trim_galore_output

outDir=dir/samples/kallisto

for ID in SCC7174-E SCC7174-F SCC7174-L SCC7184-E SCC7184-F SCC7184-T SCC7190-F SCC7190-T SCC7195-E SCC7195-F SCC7195-L SCC7195-T SCC7174-T SCC7184-L SCC7190-E SCC7197-E SCC7197-F SCC7197-L SCC7197-T SCC7200-E SCC7200-L SCC7202-E SCC7202-F SCC7202-L SCC7202-T SCC7204-L SCC7204-T SCC7205-E SCC7205-F SCC7205-L SCC7205-T SCC7207-E SCC7207-L SCC7207-T SCC7190-L SCC7209-L SCC7209-T SCC7211-E SCC7211-F SCC7211-L SCC7211-T SCC7224-F SCC7224-L SCC7224-T SCC7228-E SCC7228-F SCC7228-L SCC7228-T SCC7151-T SCC7151-L
do
echo "Running kallisto on $ID"
$kallisto quant -i $refDir/kallisto.gencodev34.idx -o $outDir/${ID} ${samplesDir}/${ID}${IDsuffix1}.fq.gz ${samplesDir}/${ID}${IDsuffix2}.fq.gz &> $outDir/${ID}_kallisto.log
done

##############################
#run multiQC on samples before and after pseudoalignment 
###############################

cd dir/python-virtual-environments
source env/bin/activate

TrimGaloreoutDir=dir/samples/trim_galore_output

multiqc $TrimGaloreoutDir1 $TrimGaloreoutDir2  $outDir 

######################
#Read in aligned reads and convert to gene-level raw expression data using tximport
########################

R
library(tximport)

#get tx2gene file
tx2g=readRDS(paste0(dir, "gencode.v34.tx2gene.rds"))

dir=dir/samples/kallisto
allfiles=list.files(dir)

Kallisto.txim = tximport(allfiles, type="kallisto",tx2gene = tx2g)

saveRDS(Kallisto.txim, paste0(dir, "Kallisto.txim.gencode34.U54.all.rds"))
write.csv(Kallisto.txim,file=paste0(dir, "Kallisto.txim.gencode34.U54.all.txt"),row.names = T)

################################################
#Apply DeSeq2 to normlize the counts 
#################################################

library(DESeq2)

#Get clinical data for patient samples 
patients=readRDS(paste0("U54_data_clinical_formatted_all.rds"))
sampleTable=patients[match(colnames(Kallisto.txim$counts), rownames(patients)),]

#make normalized counts. Normalizing for batch
ddsb=DESeqDataSetFromTximport(Kallisto.txim,sampleTable, design=~batch)
ddsb=ddsb[rowSums(counts(ddsb))>1, ]

#Make normalized counts
dds=estimateSizeFactors(ddsb)
normalized_counts=counts(dds, normalized=TRUE)

saveRDS(normalized_counts, paste0(dir, "U54.normalized.counts.allsamples.gene.names.rds"))


##################################################################
##################################################################
#Processing and basic analysis of the Stanford scRNA-Seq dataset
##################################################################
##################################################################

#Starting from cellranger-aligned read files

#srun -p gpu --gres=gpu:2 --time=3:0:0 -c 1 --pty bash 
ml load R/3.6.1
R

library(Seurat)
library(scran)

#Get Cell Ranger-aligned files and read into Seurat
dir="/dir/cellranger_aligned/"
files=list.files(dir)

ddirs=paste0(dir, files, "/outs/filtered_feature_bc_matrix/")

#Read into Seurat object and generate a sparse matrix
readall=lapply(ddirs, function(ddir) Read10X(data.dir = ddir, gene.column = 2, unique.features = TRUE))
names(readall)=files

d10x.data <- sapply(1:length(ddirs), function(x){
  d10x <- Read10X(data.dir = ddirs[x], gene.column = 2, unique.features = TRUE)
  colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),ids[x],sep="-")
  d10x
})
experiment.data <- do.call("cbind", d10x.data)

experiment.aggregate <- CreateSeuratObject(
  experiment.data,
  project = "scRNA CCSB HNSC", 
  min.cells = 10,
  min.genes = 200,
  names.field = 2,
  names.delim = "\\-")

saveRDS(experiment.aggregate, paste0(dir, "CCSB.HNSC.scRNASeq.raw.Seurat.object.rds"))

#Calculate the proportion of mitochondrial RNA
mito.genes <- grep("^MT-", rownames(experiment.aggregate), value = T)
percent.mito <- Matrix::colSums(experiment.aggregate[mito.genes, ]) / Matrix::colSums(experiment.aggregate)

# AddMetaData adds columns to object@data.info, and is a great place to stash QC stats
experiment.aggregate <- AddMetaData(
  object = experiment.aggregate,
  metadata = percent.mito,
  col.name= "percent.mito")

ids=paste0("scc", gsub(".*scc", "", colnames(experiment.aggregate)))
experiment.aggregate=AddMetaData(experiment.aggregate, ids, col.name = "sample")

#add metadata indicating the dissociation protocol and primary/lpymph node metastasis status of the tumor based on the sample names
experiment.aggregate=AddMetaData(experiment.aggregate, ifelse(gsub(".*-", "", md2$sample)=="E", "Enzymatic","Mechanical"), col.name = "dissociation")
location=factor(ifelse(substr(gsub("-.*", "", md2$sample), nchar(gsub("-.*", "", md2$sample)), nchar(gsub("-.*", "", md2$sample)))=="B", "Met", "Primary"), levels=c("Primary", "Met"))
experiment.aggregate=AddMetaData(experiment.aggregate, location, col.name = "location")

#Removing cells with % mitochondrial RNA >0.035
keep=colnames(experiment.aggregate[,experiment.aggregate$percent.mito<0.35])
experiment.aggregate.keep=subset(experiment.aggregate, cells=keep)

#removing cells with <200 or greater tha 4000 unique features
low.det.v3 <- WhichCells(experiment.aggregate.keep, expression = nFeature_RNA < 200)
high.det.v3 <- WhichCells(experiment.aggregate.keep, expression = nFeature_RNA > 4100)

experiment.aggregate=subset(experiment.aggregate.keep, cells=setdiff(WhichCells(experiment.aggregate.keep),c(low.det.v3,high.det.v3)))

#Restrict to primary HNCs that were processed using enzymatic dissociation 
primE=experiment.aggregate[,experiment.aggregate$location=="Primary" & experiment.aggregate$dissociation=="Enzymatic"]
primE$sample=factor(primE$sample)

#Normalize data for each sample and generate an integrated object using the "standard" seurat analysis workflow
samples=levels(as.factor(primE$sample))

samples.list=SplitObject(primE, split.by = "sample")
names(samples.list)=samples

reference.list=samples.list[samples]

for (i in 1:length(reference.list)) {
    reference.list[[i]] <- NormalizeData(reference.list[[i]], verbose = FALSE)
    reference.list[[i]] <- FindVariableFeatures(reference.list[[i]], selection.method = "vst", 
        nfeatures = 2000, verbose = FALSE)
}

#Ppecify the number of neighbors to use when filtering anchors
k.filter <- min(200, min(sapply(reference.list, ncol)))
patient.anchors <- Seurat::FindIntegrationAnchors(reference.list, k.filter = k.filter)
#

#Integrate samples into a single Seurat object
int=IntegrateData(anchorset = patient.anchors, k.weight=k.filter)

#Run PCA
all.genes = rownames(int)
int=ScaleData(int, features = all.genes)
int=RunPCA(int, npcs=50)

Figuresdir=""
file=paste0(Figuresdir, "Plot.elbowplot.primary_enzymatic_50PCs")
png(file=paste(file,'.png',sep=''), units="in", width=6, height=6, res=600)
ElbowPlot(int)
dev.off()
#Visualize elbow plot to determine the best number of principa components

##
npcs=50 #In practice, Gentles lab members have observed more parsimonious clustering when using a larger number of pcs. Using 50
int = FindNeighbors(int, dims = 1:npcs)

#Finding unsupervised cell clusters at various different levels of resolution. 
#For downstream anaysis, will analyze the number that best fits the data based on visualization of umap and tsne plots
int.res.0.3 = FindClusters(int, resolution = 0.3)
int.res.0.4 = FindClusters(int, resolution = 0.4)
int.res.0.5 = FindClusters(int, resolution = 0.5)
int.res.0.6 = FindClusters(int, resolution = 0.6)
int.res.0.7 = FindClusters(int, resolution = 0.7)
int.res.0.8 = FindClusters(int, resolution = 0.8)
int.res.0.9 = FindClusters(int, resolution = 0.9)
int.res.1 = FindClusters(int, resolution = 1)

#add highres clusters as metadat to Seurat object
int@meta.data[,"FindClusters.res0.3"]=Idents(int.res.0.3)
int@meta.data[,"FindClusters.res0.4"]=Idents(int.res.0.4)
int@meta.data[,"FindClusters.res0.5"]=Idents(int.res.0.5)
int@meta.data[,"FindClusters.res0.6"]=Idents(int.res.0.6)
int@meta.data[,"FindClusters.res0.7"]=Idents(int.res.0.7)
int@meta.data[,"FindClusters.res0.8"]=Idents(int.res.0.8)
int@meta.data[,"FindClusters.res0.9"]=Idents(int.res.0.9)
int@meta.data[,"FindClusters.res1"]=Idents(int.res.1)

#Run umap abd tSNE
library(umap)
int=RunUMAP(int, dims = 1:npcs)
int=RunTSNE(int, dims = 1:npcs)


#Annotate cell cycle phase
DefaultAssay(int)="RNA"

s.genes=cc.genes$s.genes
g2m.genes=cc.genes$g2m.genes

int=CellCycleScoring(int, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)


#Add EMT score
ga=as.data.frame(as.matrix(GetAssayData(object = int)))

genes_EMT_epi=c("CDH1","DSP","TJP1")
genes_EMT_mes=c("VIM","CDH2","FOXC2","SNAI1","SNAI2","TWIST1","GSC","FN1","ITGB6","MMP2","MMP3","MMP9","SOX10")

emt=colSums(ga[genes_EMT_mes,], na.rm=T)-colSums(ga[genes_EMT_epi,], na.rm=T)
int$EMT.score=emt


#Add CyTOTRACE score (Requires identification of malignant cells)
cells=colnames(intfull[,intfull$cell.type=="Malignant"])
#1494 cells
epi=subset(intfull, cells=cells)

#On linux server, navigate into python virtual environment
cd /home/kbren/python3-virtual-environments
source env/bin/activate

R 

Sys.setenv(RETICULATE_PYTHON="/home/kbren/python3-virtual-environments")
library(CytoTRACE)
library(reticulate)
library(Seurat)

counts=epi@assays$RNA@counts
counts=as.matrix(counts)
results=CytoTRACE(counts)
CytoTRACE=results$CytoTRACE
  
int@meta.data$CytoTRACE.epithelial=CytoTRACE[match(colnames(int), names(CytoTRACE))]

saveRDS(int, paste0(dir, "CCSB_scRNASeq_HNSCC_primary_enzymatic_intgrated_50PCs.rds"))
#This is the integrated object to which all downstream analysis of scRNA-Seq data will be applied


############################################################
############################################################
#Downloading and processing the Puram scRNA-Seq dataset
############################################################
############################################################

library(GEOquery)

gse="GSE103322"
#Sys.setenv(VROOM_CONNECTION_SIZE=500072)
getSet=getGEO(gse)
info=pData(getSet[[1]])
dat=exprs(getSet[[1]])

#gene expression matrix is in supp files gse="GSE103322"
cellDir=""
dir.create(cellDir)

gse="GSE103322"
filePaths = getGEOSuppFiles(gse, baseDir = cellDir)
dir=paste0(cellDir, "/","GSE103322","/")
setwd(dir)
gunzip("GSE103322_HNSCC_all_data.txt.gz")

library(R.utils)
gunzip("GSE103322_HNSCC_all_data.txt.gz", remove=FALSE)

supp=read.table(paste0(dir, "GSE103322_HNSCC_all_data.txt"), sep="\t", header=T,  fill=T, quote="", comment.char="")
rownames(supp)=supp$X
supp=supp[,2:ncol(supp)]

supp.pheno=supp[1:5,]
supp.pheno=as.data.frame(t(supp.pheno))

exp=supp[6:nrow(supp),]
rownames(exp)=gsub("'","",rownames(exp))

##
supp.pheno=data.frame(ID=rownames(supp.pheno), supp.pheno)
rownames(supp.pheno)=NULL
colnames(supp.pheno)=c("ID", "enzyme", "lymph_node", "malignant", "non.malignant", "CellType")
pheno$sample=gsub("_.*","",supp.pheno$ID)

pheno$sample2=gsub("HN|HNSCC","MEEI", pheno$sample)
length(levels(as.factor(pheno$sample2)))

dir=""
info.Puram1=readRDS(paste0(dir, Puram.clinical.info.rds"))
  
ap=info.Puram1[,c("COVAR_N_status","COVAR_N","COV_HNSCC_study_grade","Primary_Site_side_omitted","COVAR_sex", "COVAR_age")]

pheno=cbind(pheno, ap[match(pheno$sample2, rownames(ap)),])
rownames(pheno)=pheno$ID

#remove sample "MEEI" as I can't find clinical data to match this patient. Also very small numbers of cells
pheno=pheno[pheno$sample2!="MEEI",]
exp=exp[,match(pheno$ID, colnames(exp))]

patients=levels(as.factor(pheno$sample))

sum=summary(as.factor(pheno$sample))
#13 patients with over 200 cells 

#How about just primary cells 
sum2=summary(as.factor(pheno[pheno$lymph_node=="Non-Lymph Node","sample2"]))
sum2[sum2>200] #9 patients with at least 200 cells in primary tumor 

pat.200=names(sum2[sum2>200])

#################################################
#Making an integrate object with all samples that have at least 200 cells
################################################

library(Seurat)
library(scran)

exobj = CreateSeuratObject(counts = exp, project = "Puram_all_with_ln", min.cells = 10, min.features = 200, meta.data = pheno)
summary(as.factor(exobj$sample2))

#Restricting to samples with at least 50 cells overall. This excludes MEEI7, which has 7 cells
patients=names(which(summary(as.factor(exobj$sample2))>=50))
#
patientcells=exobj@meta.data[exobj@meta.data$sample2 %in% patients, "ID"]
exobj=subset(exobj, cells=patientcells)

#Adding supplementary data (Cinical data) to integrated objecy
supp.match=supp[,match(colnames(exobj), colnames(supp))]
rownames(supp.match)[1:5]=gsub(" ",".",gsub("  ","",rownames(supp.match)[1:5]))

add.meta=rownames(supp.match)[1:5]

for(i in add.meta){
  x=as.character(supp.match[i,])
  exobj@meta.data[,paste0("Puram_supp_", i)]=x
}

#restrict to primary tumor cells
exobj.prim=exobj[,exobj$lymph_node=="Non-Lymph Node"]
#4275 cells

#make integrated object with all patients that have at least 200 cells
exobj.prim.200=exobj.prim[,exobj.prim$sample2 %in% pat.200]
#3524 cells

#make integrated object
patients=levels(as.factor(exobj$sample2))
exp.list=SplitObject(exobj.prim.200, split.by = "sample2")

for (i in 1:length(exp.list)) {
    exp.list[[i]] <- NormalizeData(exp.list[[i]], verbose = FALSE)
    exp.list[[i]] <- FindVariableFeatures(exp.list[[i]], selection.method = "vst", 
        nfeatures = 2000, verbose = FALSE)
}

k.filter <- min(200, min(sapply(exp.list, ncol)))

patient.anchors <- FindIntegrationAnchors(object.list = exp.list, k.filter = k.filter)

exp.integrated <- IntegrateData(anchorset = patient.anchors, k.weight=k.filter)

int=exp.integrated

#Making cell type annotation that indicated unclassifed cells (That were not classified either as malignant or non-malignant by Puram et al.)
int$Puram_supp_unclassified=factor(ifelse(which(pheno$non.malignant=="0" & pheno$malignant=="Non Malignant","1","0"))
int$cell_type_with_unclassified=factor(ifelse(int$Puram_supp_unclassified==1, "Unclassified", as.character(int$CellType)))
#

#Running PCA
all.genes = rownames(int)
int=ScaleData(int, features = all.genes)
int=RunPCA(int, npcs=50)

file=paste0(dir, "Elbowplot.integated_Puram_prim_patients200cells")
png(file=paste(file,'.png',sep=''), units="in", width=6, height=6, res=600)
ElbowPlot(int)
dev.off()

#Decided on 20 PCs based on the Elbow plot
npcs=20

int = FindNeighbors(int, dims = 1:npcs)

####
int.res.0.3 = FindClusters(int, resolution = 0.3)
int.res.0.4 = FindClusters(int, resolution = 0.4)
int.res.0.5 = FindClusters(int, resolution = 0.5)
int.res.0.6 = FindClusters(int, resolution = 0.6)
int.res.0.7 = FindClusters(int, resolution = 0.7)
int.res.0.8 = FindClusters(int, resolution = 0.8)
int.res.0.9 = FindClusters(int, resolution = 0.9)
int.res.1 = FindClusters(int, resolution = 1)

#add highres clusters as metadat to Seurat object
int@meta.data[,"FindClusters.res0.3"]=Idents(int.res.0.3)
int@meta.data[,"FindClusters.res0.4"]=Idents(int.res.0.4)
int@meta.data[,"FindClusters.res0.5"]=Idents(int.res.0.5)
int@meta.data[,"FindClusters.res0.6"]=Idents(int.res.0.6)
int@meta.data[,"FindClusters.res0.7"]=Idents(int.res.0.7)
int@meta.data[,"FindClusters.res0.8"]=Idents(int.res.0.8)
#Names match colnames(patients.integrated$RNA)

#
library(umap)
int=RunUMAP(int, dims = 1:npcs)

#save the integrated object
saveRDS(int, paste0(dir, "Integated_Puram_prim_patients200cells.rds"))

#################################################################
#################################################################
#Processing GSE26549 (Saintigny et al. Oral premalignant lesion (Leukoplakia) dataset)
#################################################################
#################################################################

gse="GSE26549"
filePaths = GEOquery::getGEOSuppFiles(gse, baseDir = expdir)

#processing all hgu133plus2 datasets
library(oligo)
library(hugene10sthsentrezg.db)
library(hugene10sthsentrezgcdf)
library(hugene10sthsentrezgprobe)

celdir=paste0(celdir, "/GSE26549/")
setwd(celdir)
untar("GSE26549_RAW.tar")

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


