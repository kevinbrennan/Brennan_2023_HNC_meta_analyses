
##########################
##########################
#Analysis GSE26549 (Saintigny et al. Oral premalignant lesion (Leukoplakia) dataset)
#For dataset processing, see script 'Processing gene expression datasets'
##########################
##########################

library(RColorBrewer)

allinfo=readRDS(paste0(dir, "processed_exp_z_scores_GSE26549.rds"))
exp=allinfo$mas5.z
info=allinfo$info

lMN.sigs=makesig(genetype="gene", group="lnm.gene.cluster.phenograph")
lMN.sigs=lapply(lMN.sigs, function(x) x[x %in% rownames(exp)])
info$histology=plyr::revalue(factor(info$`histology at baseline (breakdown):ch1`), c("NA"=NA))
info$OCFS_Time=as.numeric(as.character(info$`oral cancer-free survival time (years):ch1`))*12
info$OS_Status=as.numeric(factor(info$`outcome:ch1`))

sigL4=lMN.sigs$L4
sigL1=lMN.sigs$L1

info$`oral cancer-free survival time (years):ch1`=as.numeric(info$`oral cancer-free survival time (years):ch1`)

df=cbind(info[which(info$`oral cancer-free survival time (years):ch1`<=5),c("OCFS_Time", "OS_Status")])

library(survival)
Survival_Object=Surv(df$OCFS_Time,df$OS_Status)

sigs=c(
makesig(genetype="gene",group="survival.direction"),
makesig(genetype="gene",group="survival.gene.cluster.phenograph"),
makesig(genetype="gene",group="lnm.direction"),
makesig(genetype="gene",group="lnm.gene.cluster.phenograph"))

sigs=lapply(sigs, function(x) x[x %in% rownames(exp)])

sigtab=as.data.frame(abind::abind(lapply(sigs, function(sig) colMeans(exp[sig,])), along=2))
all(rownames(sigtab)==rownames(info))

sigtab=cbind(sigtab, info)

cols=RColorBrewer::brewer.pal(4,"YlOrRd")

L4=lMN.sigs$L4
L1=lMN.sigs$L1

#Linear regression to test the association of gene clusters L1 and L4 with histology
summary(glm(as.numeric(info$histology)~colMeans(exp[sigL4,])))
#4.02e-12 
summary(glm(as.numeric(info$histology)~colMeans(exp[sigL1,])))
#0.039

#Plot levels of L1 and L4 gene signatures in samples stratified by histology
t1=sigtab[,c("L1", "histology")]
colnames(t1)=c("Exp","Histology")
t1$Signature=rep("L1", nrow(t1))

t2=sigtab[,c("L4", "histology")]
colnames(t2)=c("Exp","Histology")
t2$Signature=rep("L4", nrow(t2))

tab=as.data.frame(rbind(t1, t2))
tab=na.omit(tab)

tab$Signature2=plyr::revalue(as.factor(tab$Signature), c("L1"="Anti-LNM cluster L1", "L4"="Pro-LNM cluster L4"))

p=ggplot(tab, aes(x=Histology, y=Exp)) + geom_boxplot(outlier.shape=NA, aes(fill=Histology)) + geom_jitter(shape=16, size=0.4, position=position_jitter(0.2)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + labs(x = "Histology (Premalignant disease stage)", y="Scaled mean gene expression")+ scale_fill_manual(name = "Histology", values = c("#FFFFB2", "#FECC5C", "#FD8D3C", "#E31A1C"), labels = c("hyperplasia" = "Hyperplasia", "mild dysplasia" = "Mild dysplasia", "moderate dysplasia" = "Moderate dysplasia", "severe dysplasia"="Severe dysplasia"))+ theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+ facet_grid(~Signature2)
p=p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white", colour = "black",size = 0.5, linetype = "solid"))

file=paste0(Figuresdir, "boxplot_GSE26549_oral_dysplasia_grade_ggplot_long_labels")
pdf(file=paste0(file,'.pdf',sep=''), height = 2.8,  width = 5, family = "Helvetica")
p
dev.off()

#######################################################################################
#survival analysis to test the association of lymph node metastasis gene score with oral cancer-free survival 
#######################################################################################

library(survival)

Survival_Object=Surv(sigtab$OCFS_Time, sigtab$OS_Status)

#LNM gene score is defined as the mean of genes within cluster L4 minus the mean of genes within cluster L1
var=sigtab$L4-sigtab$L1

myBreaks= quantile(var,c(0,.5,1))
Sigdicot = cut(var,
              breaks=myBreaks, 
              include.lowest=TRUE,
              labels=c("1","2"))
names(Sigdicot)=names(var)

#Coxph model
model=coxph(Survival_Object~Sigdicot)
summary(model) #p=0.038
#
exp(model$coefficients) #HR= 2.104392 

#Make kaplan meier plot
mfit= survfit(Survival_Object ~ Sigdicot)

file=paste0(Figuresdir, "Kaplan_OC_free_survival_median_cluster_L4minusL1_GSE26549")
pdf(file=paste0(file,'.pdf',sep=''), height = 5,  width = 5, family = "Helvetica")
plot(mfit, col=c(1:2), ylab="Survival probability", xlab="Months to event (OC diagnosis or last follow-up)")
legend("bottomleft", legend=c("Low", "High"),col=c("black", "red"), lty=1, cex=0.8, title="LNM gene score level\n(Dichotomized)", bty="n")
dev.off()

#Plotting correlation between LNM gene clusters L1 and L4 in OPLs
#Define color vector to represent histological categories
sigtab$histology_col=as.character(plyr::revalue(sigtab$histology, c("hyperplasia"="#FFFFB2", "mild dysplasia"="#FECC5C", "moderate dysplasia"="#FD8D3C", "severe dysplasia"="#E31A1C")))

colors=RColorBrewer::brewer.pal(4, "RdBu")

df=sigtab[,c("L1","L4","histology")]
df=na.omit(df)

cols=RColorBrewer::brewer.pal(4,"YlOrRd")

p=ggplot(df, aes(x=L1, y=L4, color=histology)) +geom_point(size=0.8, shape=19)+ scale_color_manual(values=cols, name = "Histology", labels = c("hyperplasia" = "Hyperplasia", "mild dysplasia" = "Mild dysplasia", "moderate dysplasia" = "Moderate dysplasia", "severe dysplasia"="Severe dysplasia"))+geom_smooth(method=lm, se=FALSE, linetype="dashed",color="black")+labs(x="Scaled mean L1 gene expression",y="Scaled mean L4 gene expression")

file=paste0(Figuresdir, "Scatter_L1_Vs_L4_OPLs")
png(file=paste(file,'.png',sep=''), units="in", width=5, height=4, res=600)
p + annotate(geom="text", x=-1, y=-0.75, label=paste0("r=",rval),color="black")
dev.off()
















