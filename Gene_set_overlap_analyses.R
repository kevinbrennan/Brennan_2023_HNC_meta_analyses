
#########################################################################
#########################################################################
#Gene set overlap analysis (Hypergeometric analysis applied to prognostic signatures) using the msigdbr package
#########################################################################
#########################################################################

library(msigdbr)

genesets=msigdbr(species = "Homo sapiens", category = NULL, subcategory = NULL)
#Restrict to useful (Adequately annotated) collections of gene sets 
genesets=genesets[genesets$gs_cat %in% c("C1", "C2", "C5", "C6", "H") & !genesets$gs_subcat %in% c("CGP"),]

#Get all genes that were represented on all gene expression datasets included on the meta-analysis (The gene 'universe')
arraygenes=readRDS(paste0(dir, "Precog.meta.analysis.genes.all.HNSCC.rds"))

gs2=as.data.frame(genesets)
arraysymbols=as.character(gs2[match(arraygenes, gs2$entrez_gene),"human_gene_symbol"])
arraysymbols=arraysymbols[!is.na(arraysymbols)]

#Get prognostic gene signatures 
siglist=
c(makesig(genetype="hgnc_symbol",group="lnm.direction"),
makesig(genetype="hgnc_symbol",group="lnm.gene.cluster.phenograph"),
makesig(genetype="hgnc_symbol",group="survival.direction"),
makesig(genetype="hgnc_symbol",group="survival.gene.cluster.phenograph"))

gs=data.frame(term=genesets$gs_name, gene=genesets$gene_symbol)
length(levels(as.factor(gs$term)))#18437

gs.res=lapply(siglist, function(x) clusterProfiler::enricher(gene = x, TERM2GENE = gs, universe=arraysymbols))

saveRDS(gs.res, paste0(Resultsdir, "Prognostic.sigs.msigdbr_removed_C7.rds"))

gs.res2=lapply(gs.res, function(x) x[1:10,])

for(i in 1:length(gs.res2)){
  gs.res2[[i]]$Prognostic_signature=rep(names(gs.res2)[i], nrow(gs.res2[[i]]))
}dat=as.data.frame(abind::abind(gs.res2, along=1))
dat=dat[!is.na(dat$qvalue),]
rownames(dat)=NULL
dat=dat[,setdiff(colnames(dat), c("Description", "Count"))]
dat=dat[,c(8,1:7)]

#Adding infromation for datasets that was downloaded as an xml file from the 'Downloads' page of the MSigDB website
require(XML)
data=xmlParse(paste0(Datadir, "msigdb_v7.4.xml"))
xml_data=xmlToList(data)

xml_data_list=lapply(xml_data, function(x) gsub("  ", "", as.character(x), fixed = TRUE)[16])
xml_data_list_names=lapply(xml_data, function(x) gsub("  ", "", as.character(x), fixed = TRUE)[1])
names(xml_data_list)=unlist(xml_data_list_names)

code_list=lapply(xml_data, function(x) gsub("  ", "", as.character(x), fixed = TRUE)[12])
code_list_names=lapply(xml_data, function(x) gsub("  ", "", as.character(x), fixed = TRUE)[1])
names(code_list)=unlist(code_list_names)

subcategory_list=lapply(xml_data, function(x) gsub("  ", "", as.character(x), fixed = TRUE)[13])
names(subcategory_list)=unlist(code_list_names)

MEMBERS_EZID_list=lapply(xml_data, function(x) gsub("  ", "", as.character(x), fixed = TRUE)[21])
names(MEMBERS_EZID_list)=unlist(code_list_names)
MEMBERS_EZID_list=lapply(MEMBERS_EZID_list, function(x) unlist(strsplit(x,",")))

names(xml_data)=names(MEMBERS_EZID_list)

category_list=lapply(xml_data, function(x) gsub("  ", "", as.character(x), fixed = TRUE)[12])
names(category_list)=unlist(code_list_names)

dat$`Short description`=xml_data_list[dat$ID]
dat$`Short description`=as.character(dat$`Short description`)

dat=dat[,c(1:7,9,8)]
saveRDS(dat, paste0(Resultsdir, "table.gsea-prog.sigs.msigdbr.rds"))

dat$GeneRatio=paste0(dat$GeneRatio,".")
dat$BgRatio=paste0(dat$BgRatio,".")

write.table(dat, file=paste0(Resultsdir, "table.gsea-prog.sigs.msigdbr.txt"), sep="\t", quote=F, row.names = F)

############################################################
#Hypergeometric analysis to test for overlap between prognostic gene signatures and previously reported gene sets related to P53 and cell cycle (Shown in Figure 3)
##########################################################

#Get prognostic gene signatures 
sigs=c(
makesig(genetype="gene",group="survival.direction"),
makesig(genetype="gene",group="survival.gene.cluster.phenograph"),
makesig(genetype="gene",group="lnm.direction"),
makesig(genetype="gene",group="lnm.gene.cluster.phenograph"))

#Getting vector of all genes that were included in gene expression datasets that were used in the meta-analyses 
arraygenes=readRDS("~/Documents/Projects/HNSCC_PreCog/results/Precog.meta.analysis.genes.all.HNSCC.rds")
##Get p53 target genes from 

#Get sets that will be compared with prognostic signatures using hypergeometric analysis 

#Look at overlap with P53 targets #Fischer 2017
p53=as.data.frame(readxl::read_xlsx(paste0(Datadir, "Fischer_TP53_targets_highconfidence.xlsx")))
#Convert gene symbiles to entrez IDs
conv=na.omit.char(mapIds(org.Hs.eg.db, keys=p53$`Gene Symbol`, column="ENTREZID", keytype="SYMBOL", multiVals="first"))
p53.ENTREZID=conv

#Get cell cycle (Periodic) genes 
per=as.data.frame(readxl::read_xlsx(paste0(Datadir, "Dominguez_periodic_genes.xlsx")))
conv=mapIds(org.Hs.eg.db, keys=per$Symbol, column="ENTREZID", keytype="SYMBOL", multiVals="first")
per$ENTREZID=conv

G2M.ENTREZID=na.omit.char(per[which(per$`Core 67`=="Yes" & per$Phase=="G2-M"),"ENTREZID"])
G1S.ENTREZID=na.omit.char(per[which(per$`Core 67`=="Yes" & per$Phase=="G1-S"),"ENTREZID"])

#Get P53-activated genes
p53=as.data.frame(readxl::read_xlsx(paste0(Datadir, "Fischer_TP53_targets_highconfidence.xlsx")))
conv=na.omit.char(mapIds(org.Hs.eg.db, keys=p53$`Gene Symbol`, column="ENTREZID", keytype="SYMBOL", multiVals="first"))
p53.ENTREZID=conv

#Make list of gene sets
siglist=list(p53.DREAM=p53d.ENTREZID, P53.activated=p53.ENTREZID, G1S.markers=G1S.ENTREZID, G2M.markers=G2M.ENTREZID)
siglist=lapply(siglist, as.character)

#Run hypergeometric analysis using the phyper function 
pvals=lapply(siglist, function(genes2) unlist(lapply(sigs, function(genes1) phyper(length(intersect(genes1,genes2))-1,length(intersect(genes1, arraygenes)),length(setdiff(arraygenes, genes1)),length(genes2),lower.tail=F, log.p = FALSE))))
pvals=abind::abind(pvals, along=1)

#Get numbers and percentages of overlapping genes
mat.percent.prognostic.sig=lapply(siglist, function(genes2) unlist(lapply(sigs, function(genes1) length(intersect(genes1, genes2))/length(genes1)*100)))
mat.percent.prognostic.sig=abind::abind(mat.percent.prognostic.sig, along=1)
mat.ns=abind::abind(lapply(siglist, function(genes2) unlist(lapply(sigs, function(genes1) length(intersect(genes1, genes2))))), along=1)

df=data.frame(pvals=pvals, percent.prognostic.sig=mat.percent.prognostic.sig, prognostic.sig=names(pvals), n=mat.ns)
df$gene.sig=unlist(lapply(names(siglist), function(x) rep(x, length(sigs))))
df$gene.sig=factor(plyr::revalue(as.factor(df$gene.sig), c("p53.DREAM"="P53-DREAM-targets", "P53.activated"="P53-activated", "G1S.markers"="G1/S-markers", "G2M.markers"="G2/M-markers")), levels=c("G1/S-markers", "G2/M-markers", "P53-DREAM-targets", "P53-activated"))
df$prognostic.sig=factor(df$prognostic.sig, levels=c("Anti-survival", "Pro-survival", "S1", "S2", "S3", "S4", "S5", "S6", "Anti-LNM", "Pro-LNM", "L1", "L2", "L3", "L4", "L5", "L6"))
df$pvals=-log10(df$pvals)

#Make dotplot showing signficance and extent of overlap between prognostic gene signatures and P53/cell cycle-related genesets (Included in Figure 3A)
p=ggplot(df, aes(prognostic.sig, pvals))+geom_point(aes(color = n, size = percent.prognostic.sig))+ facet_wrap(~gene.sig)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+labs(y="-log10 p-value (Hypergeomteric test)", x="Prognostic signature")
p=p + scale_color_gradient(name="N overlapping\ngenes")+ scale_size_continuous(name="% prognostic\nsignature\ngenes")+ ylim(0, 47)
p=p + geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red")
#remove grey background
p=p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white", colour = "black",size = 0.5, linetype = "solid"))

file=paste0(Figuresdir, "Enrichment_plot_pval_p53_sigs_with_P53_activated")
pdf(file=paste0(file,'.pdf',sep=''), height = 6,  width = 6, family = "Helvetica")
p
dev.off()














