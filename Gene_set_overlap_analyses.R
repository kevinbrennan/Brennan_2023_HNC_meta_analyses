library(msigdbr)
source("~/Documents/Projects/HNSCC_PreCog/data/Pecog_HNSCC_directories.R")
source(paste0(RscriptsPath, "basic_scripts.R"))

require(XML)
#data <- xmlParse(paste0(Datadir, "msigdb_v7.5.1.xml"))
data <- xmlParse(paste0(Datadir, "msigdb_v2023.1.Hs.xml"))
xml_data=xmlToList(data)

xml_data_list=lapply(xml_data, function(x) gsub("  ", "", as.character(x), fixed = TRUE)[16])
xml_data_list_names=lapply(xml_data, function(x) gsub("  ", "", as.character(x), fixed = TRUE)[1])
names(xml_data_list)=unlist(xml_data_list_names)

code_list=lapply(xml_data, function(x) gsub("  ", "", as.character(x), fixed = TRUE)[12])
code_list_names=lapply(xml_data, function(x) gsub("  ", "", as.character(x), fixed = TRUE)[1])
names(code_list)=unlist(code_list_names)
subcategory_list=lapply(xml_data, function(x) gsub("  ", "", as.character(x), fixed = TRUE)[13])
names(subcategory_list)=unlist(code_list_names)
#gene ontology categories are C5_BP (biological process, n=527), C5_CC (cellular component, n=159), C5_MF (molecular function, n=178)

MEMBERS_EZID_list=lapply(xml_data, function(x) gsub("  ", "", as.character(x), fixed = TRUE)[21])
names(MEMBERS_EZID_list)=unlist(code_list_names)
MEMBERS_EZID_list=lapply(MEMBERS_EZID_list, function(x) unlist(strsplit(x,",")))

category_list=lapply(xml_data, function(x) gsub("  ", "", as.character(x), fixed = TRUE)[12])
names(category_list)=unlist(code_list_names)

names1=names(unlist(category_list)[which(unlist(category_list) %in% c("C1", "C2", "C5", "C6", "H"))])
#names2=names(unlist(subcategory_list)[which(unlist(subcategory_list)!="CGP")])
names2=names(unlist(subcategory_list)[which(!unlist(subcategory_list) %in% c("CGP", "CP:WIKIPATHWAYS", "VAX"))])
names=intersect(names1, names2)
#18833

MEMBERS_EZID=MEMBERS_EZID_list[names]

###################################################
#Applying fGSEA to LNM-associated genes
####################################################

#lnm.list=readRDS(paste0(Resultsdir, "Precog.meta.analysis.node.status.all.HNSCC.updated.03.30.2010.rds"))
lnm.list=readRDS(paste(Resultsdir,"Precog.meta.analysis.node.status.all.HNSCC.updated.06.09.23.rds", sep=""))

meta2=lnm.list$meta

metatab=cbind(
  unlist(lapply(meta2, function(x) x$pval)),
  unlist(lapply(meta2, function(x) x$beta)),
  unlist(lapply(meta2, function(x) x$zval)))
colnames(metatab)=c("P","Estimate","Z")
metatab=as.data.frame(metatab)

#@Running gesa based on estimate
ranks.Z=metatab$Z
names(ranks.Z)=rownames(metatab)

library(fgsea)

fgseaRes.lnm=fgsea(pathways = MEMBERS_EZID, 
                   stats    = ranks.Z,
                   minSize  = 15,
                   maxSize  = 500, 
                   eps=0)

saveRDS(fgseaRes.lnm, paste0(Resultsdir,"GSEA.LNM.genes.updated.061423.rds"))
fgseaRes.lnm=readRDS(paste0(Resultsdir,"GSEA.LNM.genes.updated.061423.rds"))

fgseaRes.lnm=fgseaRes.lnm[order(fgseaRes.lnm$pval),]
fgseaRes.lnm.sig=as.data.frame(fgseaRes.lnm[which(fgseaRes.lnm$padj<0.05),])

leading.edge.symbols=na.omit.char(mapIds(org.Hs.eg.db, keys=unlist(fgseaRes.lnm.sig$leadingEdge[1]), column="SYMBOL", keytype="ENTREZID", multiVals="first"))
#Both AURKA and AURKB in there. Look at Quah paper for mention of AURK

fgseaRes=rbind(fgseaRes.lnm[fgseaRes.lnm$NES>0,][1:10,],fgseaRes.lnm[fgseaRes.lnm$NES<0,][1:10,])

file=paste(Figuresdir, "GSEA_lnm_genes_top10_pos_neg_updated_061423", sep="")
pdf(file=paste0(file,'.pdf',sep=''), height = 4,  width = 7, family = "Helvetica")
ggplot(fgseaRes, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=-log10(padj))) +
  coord_flip() +
  labs(x="MSigDB gene set", y="Normalized Enrichment Score",
       title="Survival-associated gene enrichments") + 
  theme_minimal() + scale_fill_continuous(name = "-log10 padj")
dev.off()

p1=ggplot(fgseaRes, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=-log10(padj))) +
  coord_flip() +
  labs(x="MSigDB gene set", y="Normalized Enrichment Score",
       title="LNM-associated gene enrichments") + 
  theme_minimal() + scale_fill_continuous(name = "-log10 padj")


#fGSEA adding p53-DREAM target genes
library(AnnotationDbi)
library("org.Hs.eg.db")

#Add signatures from figure 3A
p53=as.data.frame(readxl::read_xlsx(paste0(Datadir, "Fischer_TP53_targets_highconfidence.xlsx")))

#conv=mapIds(org.Hs.eg.db, keys=per$Symbol, column="SYMBOL", keytype="SYMBOL", multiVals="first")
conv=na.omit.char(mapIds(org.Hs.eg.db, keys=p53$`Gene Symbol`, column="ENTREZID", keytype="SYMBOL", multiVals="first"))
p53.ENTREZID=conv

#p53 dream complex genes from Fischer et al 2016 :PMID: 26384566
p53d=as.data.frame(readxl::read_xlsx(paste0(Datadir, "P53_DREAM_Fisher_targets.xlsx")))
p53d=p53d[p53d$`p53-p21-DREAM-CDE/CHR target gene`=="1",]

conv=na.omit.char(mapIds(org.Hs.eg.db, keys=p53d$`Ensembl ID`, column="ENTREZID", keytype="ENSEMBL", multiVals="first"))
p53d.ENTREZID=conv

#Get cell cycle genes 
per=as.data.frame(readxl::read_xlsx(paste0(Datadir, "Dominguez_periodic_genes.xlsx")))
conv=mapIds(org.Hs.eg.db, keys=per$Symbol, column="ENTREZID", keytype="SYMBOL", multiVals="first")
per$ENTREZID=conv

G2M.ENTREZID=na.omit.char(per[which(per$`Core 67`=="Yes" & per$Phase=="G2-M"),"ENTREZID"])
G1S.ENTREZID=na.omit.char(per[which(per$`Core 67`=="Yes" & per$Phase=="G1-S"),"ENTREZID"])

siglist=list(p53=p53.ENTREZID, p53.DREAM=p53d.ENTREZID, G1S.markers=G1S.ENTREZID, G2M.markers=G2M.ENTREZID)
siglist=lapply(siglist, as.character)

MEMBERS_EZID2=c(MEMBERS_EZID, siglist)

fgseaRes.lnm.with.p53.DREAM=fgsea(pathways = MEMBERS_EZID2, 
                   stats    = ranks.Z,
                   minSize  = 15,
                   maxSize  = 500, 
                   eps=0)

saveRDS(fgseaRes.lnm.with.p53.DREAM, paste0(Resultsdir,"GSEA.LNM.genes.with.p53.DREAM.updated.061423.rds"))
fgseaRes.lnm.with.p53.DREAM=readRDS(paste0(Resultsdir,"GSEA.LNM.genes.with.p53.DREAM.updated.061423.rds"))

fgseaRes.lnm.with.p53.DREAM %>% arrange(pval) %>% head %>% as.data.frame() 

fgseaRes.lnm.with.p53.DREAM=fgseaRes.lnm.with.p53.DREAM[order(fgseaRes.lnm.with.p53.DREAM$pval),]

fgseaRes=rbind(fgseaRes.lnm.with.p53.DREAM[fgseaRes.lnm.with.p53.DREAM$NES>0,][1:10,],fgseaRes.lnm.with.p53.DREAM[fgseaRes.lnm.with.p53.DREAM$NES<0,][1:10,])

file=paste(Figuresdir, "GSEA_lnm_genes_top10_pos_neg_with_p53_DREAM_updated_061423", sep="")
pdf(file=paste0(file,'.pdf',sep=''), height = 4,  width = 7, family = "Helvetica")
ggplot(fgseaRes, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=-log10(padj))) +
  coord_flip() +
  labs(x="MSigDB gene set", y="Normalized Enrichment Score",
       title="Survival-associated gene enrichments") + 
  theme_minimal() + scale_fill_continuous(name = "-log10 padj")
dev.off()

file=paste(Figuresdir, "Enrichment_plot_p53_DREAM_LNM_updated_061423", sep="")
pdf(file=paste0(file,'.pdf',sep=''), height = 4,  width = 6, family = "Helvetica")
plotEnrichment(MEMBERS_EZID2[[head(fgseaRes.lnm.with.p53.DREAM[order(pval), ], 1)$pathway]], ranks.Z) + labs(title=head(fgseaRes.lnm.with.p53.DREAM[order(pval), ], 1)$pathway)
dev.off()


###################################################
#Applying fGSEA to survival-associated genes
####################################################

#PRECOG.HNSCC.All=readRDS(paste0(Resultsdir, "Precog.HNSCC.coxph.survival.allstudies.allHNSCC.updated.06_09_2023.rds"))
metaz.coxph.annot=readRDS(paste(Resultsdir, "Precog.HNSCC.coxph.Liptak.metaz.allHNSCC.updated.06_09_2023.rds", sep=""))

ranks.Z=metaz.coxph.annot$metaz.coxph
names(ranks.Z)=rownames(metaz.coxph.annot)

fgseaRes.Z=fgsea(pathways = MEMBERS_EZID, 
                 stats    = ranks.Z,
                 minSize  = 15,
                 maxSize  = 500, 
                 eps=0)
fgseaRes.Z=fgseaRes.Z[order(fgseaRes.Z$pval),]

fgseaRes.Z.sig=as.data.frame(fgseaRes.Z[which(fgseaRes.Z$padj<0.05),])

saveRDS(fgseaRes.Z, paste0(Resultsdir,"GSEA.survival.genes.Zscore.updated.061423.rds"))
fgseaRes.Z=readRDS(paste0(Resultsdir,"GSEA.survival.genes.Zscore.updated.061423.rds"))

fgseaRes.Z.surv=fgseaRes.Z

fgseaRes.Z.surv$pathway.edited=gsub("GOBP_|GOCC_|HP_|GOMF_|REACTOME_|HALLMARK_|IMMUNE_RESPONSE_", "", fgseaRes.Z.surv$pathway)

fgseaRes.Z.surv=fgseaRes.Z.surv[fgseaRes.Z.surv$pathway!="GOBP_IMMUNE_RESPONSE_REGULATING_CELL_SURFACE_RECEPTOR_SIGNALING_PATHWAY",]

fgseaRes.Z.surv$pathway.edited=gsub("REACTOME_REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_IGF_TRANSPORT_AND_UPTAKE_BY_INSULIN_LIKE_GROWTH_FACTOR_BINDING_PROTEINS_IGFBPS","REACTOME_REGULATION_OF_IGF2_TRANSPORT_AND_UPTAKE_BY_IGFBPS*", fgseaRes.Z.surv$pathway)     

fgseaRes=rbind(fgseaRes.Z.surv[fgseaRes.Z.surv$NES>0,][1:10,],fgseaRes.Z.surv[fgseaRes.Z.surv$NES<0,][1:10,])

file=paste(Figuresdir, "GSEA_survival_genes_top10_pos_neg_updated_061423", sep="")
pdf(file=paste0(file,'.pdf',sep=''), height = 4,  width = 10, family = "Helvetica")
ggplot(fgseaRes, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=-log10(padj))) +
  coord_flip() +
  labs(x="MSigDB gene set", y="Normalized Enrichment Score",
       title="Survival-associated gene enrichments") + 
  theme_minimal() + scale_fill_continuous(name = "-log10 padj")
dev.off()

p2=ggplot(fgseaRes, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=-log10(padj))) +
  coord_flip() +
  labs(x="MSigDB gene set", y="Normalized Enrichment Score",
       title="Survival-associated gene enrichments") + 
  theme_minimal() + scale_fill_continuous(name = "-log10 padj")

library(cowplot)

both2 <- align_patches(p1, p2)
p1x <- ggdraw(both2[[1]])
p2x <- ggdraw(both2[[2]])

file=paste(Figuresdir, "GSEA_survival_genes_top10_pos_neg_aligned_short_labs_updated_061523", sep="")
pdf(file=paste0(file,'.pdf',sep=''), height = 4,  width = 8, family = "Helvetica")
p2x
dev.off()

file=paste(Figuresdir, "GSEA_lnm_genes_top10_pos_neg_aligned_short_labs_updated_061523", sep="")
pdf(file=paste0(file,'.pdf',sep=''), height = 4,  width = 8, family = "Helvetica")
p1x
dev.off()

###############################################################################
#Run GSOA for LNM and survival-associated genes 
################################################################################

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

siglist=c(makesig(genetype="symbol_AnnotationDbi",group="lnm.direction"),
          makesig(genetype="symbol_AnnotationDbi",group="lnm.gene.cluster.phenograph"))

siglist=lapply(siglist, function(x) x[nchar(x)>0])

library(msigdbr)

genesets=msigdbr(species = "Homo sapiens", category = NULL, subcategory = NULL)
genesets=genesets[genesets$gs_cat %in% c("C1", "C2", "C5", "C6", "H") & !genesets$gs_subcat %in% c("CGP"),]

#arraygenes=readRDS("~/Documents/Projects/HNSCC_PreCog/results/Precog.meta.analysis.genes.all.HNSCC.rds")
#Easier to report as symbols 
#arraygenes=readRDS(paste(Resultsdir,"Precog.meta.analysis.genes.all.HNSCC.updated.06.09.23.rds", sep=""))
lnm.list=readRDS(paste(Resultsdir,"Precog.meta.analysis.node.status.all.HNSCC.updated.06.09.23.rds", sep=""))
arraygenes=names(lnm.list$meta)

gs=as.data.frame(genesets)
arraysymbols=as.character(gs[match(arraygenes, gs$entrez_gene),"human_gene_symbol"])
arraysymbols=arraysymbols[!is.na(arraysymbols)]

gs=data.frame(term=genesets$gs_name, gene=genesets$gene_symbol)
length(levels(as.factor(gs$term)))#18993

gs.res=lapply(siglist, function(x) clusterProfiler::enricher(gene = x, TERM2GENE = gs, universe=arraysymbols))

saveRDS(gs.res, paste0(Resultsdir, "Prognostic.sigs.msigdbr_LNM_clusters.updated.06.09.23.rds"))
gs.res=readRDS(paste0(Resultsdir, "Prognostic.sigs.msigdbr_LNM_clusters.updated.06.09.23.rds"))

##########################################################
#Survival genes
##########################################################

siglist=c(makesig(genetype="symbol_AnnotationDbi",group="survival.direction"),
          makesig(genetype="symbol_AnnotationDbi",group="survival.gene.cluster.phenograph"))

siglist=lapply(siglist, function(x) x[nchar(x)>0])

genesets=msigdbr(species = "Homo sapiens", category = NULL, subcategory = NULL)
genesets=genesets[genesets$gs_cat %in% c("C1", "C2", "C5", "C6", "H") & !genesets$gs_subcat %in% c("CGP"),]

metaz.coxph.annot=readRDS(paste0(Resultsdir, "Precog.HNSCC.coxph.Liptak.metaz.allHNSCC.updated.06_09_2023.rds"))
arraygenes=metaz.coxph.annot$gene

gs=as.data.frame(genesets)
arraysymbols=as.character(gs[match(arraygenes, gs$entrez_gene),"human_gene_symbol"])
arraysymbols=arraysymbols[!is.na(arraysymbols)]

gs=data.frame(term=genesets$gs_name, gene=genesets$gene_symbol)
length(levels(as.factor(gs$term)))#18993

gs.res=lapply(siglist, function(x) clusterProfiler::enricher(gene = x, TERM2GENE = gs, universe=arraysymbols))

saveRDS(gs.res, paste0(Resultsdir, "Prognostic.sigs.msigdbr_survival_clusters.updated.06.09.23.rds"))

###############################
#Make GSOA table for supplementary table 
###############################

gs.res.surv=readRDS(paste0(Resultsdir, "Prognostic.sigs.msigdbr_survival_clusters.updated.06.09.23.rds"))
gs.res.lnm=readRDS(paste0(Resultsdir, "Prognostic.sigs.msigdbr_LNM_clusters.updated.06.09.23.rds"))

gs.res=c(gs.res.surv, gs.res.lnm)
gs.res=lapply(gs.res, as.data.frame)
n_gr=unlist(lapply(gs.res, nrow))

#Remove signatures with no enrichments 
gs.res=gs.res[-which(n_gr==0)]
n_gr=unlist(lapply(gs.res, nrow))

n_gr=ifelse(n_gr<10, n_gr, 10)

gs.res2=lapply(1:length(n_gr), function(i) {gs.res[[i]][1:n_gr[i],]})
names(gs.res2)=names(gs.res)

for(i in 1:length(gs.res2)){
  gs.res2[[i]]$Prognostic_signature=rep(names(gs.res2)[i], nrow(gs.res2[[i]]))
}

dat=as.data.frame(abind::abind(gs.res2, along=1))
dat=dat[!is.na(dat$qvalue),]
rownames(dat)=NULL
dat=dat[,setdiff(colnames(dat), c("Description", "Count"))]
dat=dat[,c(8,1:7)]

gs=as.data.frame(genesets)
dat$gs_description=gs[match(dat$ID, gs$gs_name),"gs_description"]

saveRDS(dat, paste0(Resultsdir, "table.gsea-prog.sigs.msigdbr.updated.062223.rds"))
dat=readRDS(paste0(Resultsdir, "table.gsea-prog.sigs.msigdbr.updated.062223.rds"))
dat$GeneRatio=paste0(" ", dat$GeneRatio)

write.table(dat, file=paste0(Resultsdir, "table.gsea.prog.sigs.msigdbr.updated.062223.txt"), sep="\t", quote=F, row.names = F)

