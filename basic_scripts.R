
#Basic operations

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

get.matched.genes=function(sig, object.genes){
  library(HGNChelper)
  sig.hg=checkGeneSymbols(sig, unmapped.as.na=FALSE)
  object.genes.hg=checkGeneSymbols(object.genes, unmapped.as.na=FALSE)
  matched=object.genes.hg[match(sig.hg$Suggested.Symbol, object.genes.hg$Suggested.Symbol),"x"]
  matched=unique(matched[!is.na(matched)])
  return(matched)
}

na.omit.char=function(x){
  return(x[!is.na(x)])
}

#Borrowed from Stack overflow #https://stackoverflow.com/questions/23018256/printing-p-values-with-0-001
pvalr <- function(pvals, sig.limit = .001, digits = 3, html = FALSE) {

  roundr <- function(x, digits = 1) {
    res <- sprintf(paste0('%.', digits, 'f'), x)
    zzz <- paste0('0.', paste(rep('0', digits), collapse = ''))
    res[res == paste0('-', zzz)] <- zzz
    res
  }

  sapply(pvals, function(x, sig.limit) {
    if (x < sig.limit)
      if (html)
        return(sprintf('&lt; %s', format(sig.limit))) else
          return(sprintf('< %s', format(sig.limit)))
    if (x > .1)
      return(roundr(x, digits = 2)) else
        return(roundr(x, digits = digits))
  }, sig.limit = sig.limit)
}

addmetagenes=function(sig, df, so, fn){
cgenes=sig
cgenes=cgenes[cgenes %in% rownames(ga)]
c_colMeans=as.numeric(scale(colMeans(df[cgenes, ], na.rm = T)))
featurename=fn
so@meta.data[, featurename]=c_colMeans
return(so)
}



## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


###################
#Color palettes 
####################

pal=viridisLite::viridis(n = 10, option = "C", direction = -1)

library(ggplot2)
params.white=ggplot2::theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white", colour = "black",size = 0.5, linetype = "solid"))

gnuplot_colors=c('#E41A1C','#377EB8','#4DAF4A','#984EA3','#FF7F00','#FFFF33','#A65628','#F781BF','#66C2A5','#FC8D62','#8DA0CB','#E78AC3','#A6D854','#FFD92F','#E5C494','#B3B3B3')
cbp=c("#000000", "#E69F00", "#56B4E9", "#009E73","#FD61D1", "#A3A500", "#D55E00", "#CC79A7")
cbp_phase=c("#DDCC77", "#661100", "#117733")

ggplotColours <- function(n = 1, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1)  h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

#ggplot2 to make all cell cluster colors the same. Could make this a generic function
scale_fill_kb <- function(...){
    ggplot2:::manual_scale(
        'fill', 
        values = setNames(c(ggplotColours(length(cellclusts))), cellclusts), 
        ...
    )
}

#get 60 color vector
library(RColorBrewer)
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# replace  "#BF5B17" and "#F0027F" "#FFFF99" as they are not unique
col_vector = as.character(plyr::revalue(as.factor(col_vector), c("#FFFF99"="chartreuse", "#BF5B17"="aquamarine","#F0027F"="firebrick")))


#################################
#stacked vlnplot for Seurat scripts, borrowed from https://divingintogeneticsandgenomics.rbind.io/post/stacked-violin-plot-for-visualizing-single-cell-data-in-seurat/
#################################

library(Seurat)
library(patchwork)
library(ggplot2)

modify_vlnplot<- function(obj,
feature,
pt.size = 0,
plot.margin = margin(0, 0, 0, 0, "cm"),
...) {
p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )+
ylab(feature) +
geom_boxplot(alpha=0.001) +
theme(legend.position = "none",
plot.title= element_blank(),
axis.title.x = element_blank(),
axis.text.x = element_blank(),
axis.ticks.x = element_blank(),
axis.title.y = element_text(size = rel(1), angle = 0),
axis.text.y = element_text(size = rel(1)),
plot.margin = plot.margin )
return(p)
}

#
StackedVlnPlot<- function(obj, features, title1,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(angle = 90, hjust=1), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))

  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1) + plot_annotation(title = title1, theme = theme(plot.title = element_text(size = 20, hjust=0.6)))
  return(p)
}



StackedVlnPlot.kb2=function(obj, features, title1,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_blank(), axis.ticks.x = element_blank())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))

  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1) + plot_annotation(title = title1, theme = theme(plot.title = element_text(size = 20, hjust=0.6)))
  return(p)
}


