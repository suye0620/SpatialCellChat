# Audit whether the two existing real fixtures satisfy the merged producer contract.
setwd(local({a<-commandArgs(FALSE);f<-sub('^--file=','',grep('^--file=',a,value=TRUE));if(length(f))dirname(dirname(normalizePath(f)))else getwd()}))
Sys.setenv(RENV_PATHS_LIBRARY='renv/library')
if(!nzchar(Sys.getenv('RENV_PROJECT'))) {if(requireNamespace('renv',quietly=TRUE)) renv::load(getwd()) else source('renv/activate.R')}
suppressPackageStartupMessages({library(methods);library(Matrix)})
source('R/SpatialCellChat_class.R');source('R/database.R');source('R/utilities.R');source('R/merging.R')
a<-readRDS('tests_dev/test_data/visium_pp1_chat.rds');b<-readRDS('tests_dev/test_data/slideseq_full_chat.rds')
result<-data.frame(dataset1='visium_pp1',dataset2='slideseq_mbm05',cells1=ncol(a@assay$norm),cells2=ncol(b@assay$norm),
  gene_axis_identical=identical(rownames(a@assay$norm),rownames(b@assay$norm)),
  shared_genes=length(intersect(rownames(a@assay$norm),rownames(b@assay$norm))),
  group_axis_identical=identical(levels(a@idents),levels(b@idents)),
  shared_groups=length(intersect(levels(a@idents),levels(b@idents))),
  db_identical=identical(a@DB,b@DB),stringsAsFactors=FALSE)
result$producer_result<-tryCatch({mergeSpatialCellChat(list(visium=a,slideseq=b),show.plot=FALSE);'unexpected success'},error=conditionMessage)
print(result)
out<-'tests_dev/figures/merged/real-fixture-contract-2026-11-11.csv'
write.csv(result,out,row.names=FALSE)
if(result$producer_result=='unexpected success')stop('Expected incompatible axes to be rejected')
