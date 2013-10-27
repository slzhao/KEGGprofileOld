### R code from vignette source 'KEGGprofile.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: KEGGprofile.Rnw:15-16
###################################################
  if(exists(".orig.enc")) options(encoding = .orig.enc)


###################################################
### code chunk number 2: load-the-data
###################################################
library(KEGGprofile)
data(pro_pho_expr)
data(pho_sites_count)
ls()
colnames(pro_pho_expr)


###################################################
### code chunk number 3: download_needed_files
###################################################
download_KEGGfile(pathway_id="04110",specis='hsa')


###################################################
### code chunk number 4: find-enriched-pathways
###################################################
genes<-row.names(pho_sites_count)[which(pho_sites_count>=10)]
pho_KEGGresult<-find_enriched_pathway(genes,specis='hsa')
pho_KEGGresult[[1]][,c(1,5)]


###################################################
### code chunk number 5: Visualization-bg
###################################################
## the phosphoproteome data
pho_expr<-pro_pho_expr[,7:12]
temp<-apply(pho_expr,1,function(x) length(which(is.na(x))))
pho_expr<-pho_expr[which(temp==0),]
## transform the expression difference into specific color
col<-col_by_value(pho_expr,col=colorRampPalette(c('green','black','red'))(1024),range=c(-6,6))
## visualization by method 'bg'
temp<-plot_pathway(pho_expr,type="bg",bg_col=col,text_col="white",magnify=1.2,specis='hsa',database_dir=system.file("extdata",package="KEGGprofile"),pathway_id="04110")


###################################################
### code chunk number 6: Visualization-lines
###################################################
## transform the number of phosphorylation sites into specific color
col<-col_by_value(pho_sites_count,col=colorRampPalette(c('white','khaki2'))(4),breaks=c(0,1,4,10,Inf)) ## visualization by method 'lines'
temp<-plot_pathway(pro_pho_expr,type="lines",bg_col=col,line_col=c("brown1","seagreen3"),groups=c(rep("Proteome",6),rep("Phosphoproteome",6)),magnify=1.2,specis='hsa',database_dir=system.file("extdata",package="KEGGprofile"),pathway_id="04110",max_dist=5)


