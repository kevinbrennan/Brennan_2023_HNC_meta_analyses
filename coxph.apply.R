library(metafor)

coxph.apply=function(x){
  library(survival)
  library(abind)
  
  out <- tryCatch(
  {
  cox=coxph(Surv(info[,survtimevar], info[,survstatvar])~x)
  z = as.list(coef(cox)/sqrt(diag(vcov(cox))))[[1]]
  se = summary(cox)$coefficients[3]
  n = cox$n
  df = data.frame(Z = z, SE = se, N = n, DatasetID = acc)
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

func=function(x){
  df=as.data.frame(x[match(intersectgenes, rownames(x)),])
  df$gene=rownames(df)
  return(df)
}

make.coxph.intersect=function(x){
  df=as.data.frame(x[match(intersectgenes, rownames(x)),])
  df$gene=intersectgenes
  df=df[!is.na(df$DatasetID),]
  return(df)
}




meta.coxph=function(gene) {
MetaGene=MetaFrame[MetaFrame$gene==gene,]
d=as.numeric(as.character(MetaGene$Z))
se=as.numeric(as.character(MetaGene$SE))
names1=as.character(MetaGene$DatasetID)
out <- tryCatch(
  {
  g=metafor::rma(d, sei=se, data=MetaGene)
  },
   error=function(cond) {
        message("Here's the original error message:")
            message(cond)
                 return(NA)
   }, 
   warning=function(cond) {
       message("Here's the original warning message:")
            message(cond)
            #Returning NA for error becuase it doesn't calculate a summary statistic if the algorithm doesn't converge, yet it is treated as a warning
                 return(NA)
   }
)
return(out)
}

#Version of coxph. apply if you already have the survival object 
coxph.apply.provide.survobject=function(x){
  library(survival)
  library(abind)
  
  out <- tryCatch(
  {
  cox=coxph(Survival_object~x)
  z = as.list(coef(cox)/sqrt(diag(vcov(cox))))[[1]]
  se = summary(cox)$coefficients[3]
  n = cox$n
  df = data.frame(Z = z, SE = se, N = n, DatasetID = acc)
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

#coxph.apply.provide.survobject

############################
#same as above but record p.value

#Version of coxph. apply if you already have the survival object 
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


#Version of coxph adjusted for HPV status with interaction term, to be used if provifing survival object
#Collecting coefficient to help interpret interaction term
coxph.apply.provide.survobject.with.pvalue_adj_HPV_interaction=function(x){
  library(survival)
  library(abind)
  
  out <- tryCatch(
    {
      cox=coxph(Survival_object~x*info$COV_HNSCC_HPV_status)
      coeffs=summary(cox)$coefficients
      
      coef_gene=coeffs[1]
      coef_hpv=coeffs[2]
      coef_geneXhpv=coeffs[3]
      
      z_gene=coeffs[10]
      z_hpv=coeffs[11]
      z_geneXhpv=coeffs[12]
      
      p_gene=coeffs[13]
      p_hpv=coeffs[14]
      p_geneXhpv=coeffs[15]
      
      se_gene=coeffs[7]
      se_hpv=coeffs[8]
      se_geneXhpv=coeffs[9]
      
      n = cox$n
      
      df = data.frame(z_gene, z_hpv, N = n, z_geneXhpv, p_gene, p_hpv, p_geneXhpv, se_gene, se_hpv, se_geneXhpv, coef_gene, coef_hpv, coef_geneXhpv, DatasetID = acc)
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
#

#Glm genes associated with LNM ajusted for HPV with interaction term 
glm_node_adj_HPV_with_interaction=function(x){
  library(survival)
  library(abind)
  
  out <- tryCatch(
    {
      mod=glm(info$COVAR_N_status~x*info$COV_HNSCC_HPV_status, family = "binomial")
      coeffs=summary(mod)$coefficients
      
      est_gene=coeffs[2]
      est_hpv=coeffs[3]
      est_geneXhpv=coeffs[4]
      
      z_gene=coeffs[10]
      z_hpv=coeffs[11]
      z_geneXhpv=coeffs[12]
      
      p_gene=coeffs[14]
      p_hpv=coeffs[15]
      p_geneXhpv=coeffs[16]
      
      se_gene=coeffs[6]
      se_hpv=coeffs[7]
      se_geneXhpv=coeffs[8]
      
      n=stats::nobs(mod)
      
      df = data.frame(z_gene, z_hpv, z_geneXhpv, p_gene, p_hpv, p_geneXhpv, se_gene, se_hpv, se_geneXhpv, est_gene, est_hpv, est_geneXhpv, DatasetID = acc, N = n)
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


#Coxph survival meta analysis adjusted for age and sex
coxph.apply.provide.survobject.adj.age.sex=function(x){
  library(survival)
  library(abind)
  
  out <- tryCatch(
    {
      cox=coxph(Survival_object~x + info$COVAR_age + info$COVAR_sex)
      coeffs=summary(cox)$coefficients
      
      z_gene=coeffs[10]
      z_age=coeffs[11]
      z_sex=coeffs[12]
      
      p_gene=coeffs[13]
      p_age=coeffs[14]
      p_sex=coeffs[15]
      
      se_gene=coeffs[7]
      se_age=coeffs[8]
      se_sex=coeffs[9]
      
      n = cox$n
      
      df = data.frame(z_gene, z_age, N = n, z_sex, p_gene, p_age, p_sex, se_gene, se_age, se_sex,DatasetID = acc, P=p)
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


###

#glm association of genes with LNM status adjusted for age and sex
coxph.apply.provide.survobject.adj.age.sex=function(x){
  library(survival)
  library(abind)
  
  out <- tryCatch(
    {
      cox=coxph(Survival_object~x + info$COVAR_age + info$COVAR_sex)
      coeffs=summary(cox)$coefficients
      
      z_gene=coeffs[10]
      z_age=coeffs[11]
      z_sex=coeffs[12]
      
      p_gene=coeffs[13]
      p_age=coeffs[14]
      p_sex=coeffs[15]
      
      se_gene=coeffs[7]
      se_age=coeffs[8]
      se_sex=coeffs[9]
      
      n = cox$n
      
      df = data.frame(z_gene, z_age, N = n, z_sex, p_gene, p_age, p_sex, se_gene, se_age, se_sex,DatasetID = acc)
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

#
Liptak_combine_z=function(gene, df){
  mdf=df[df$gene==gene,]
  metaz=sum(mdf$N*mdf$Z)/sqrt(sum(mdf$N^2))
  return(metaz)
}

