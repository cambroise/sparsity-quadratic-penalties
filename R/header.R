## The header function does all the job for pitiful ignorants 
## It loads R and C++ libraries and compile the C++ code when not found
rm(list=ls())
gc()
require(ggplot2)
require(Matrix)
require(RcppArmadillo)
header <- function(path.2.svn = "~/svn/quadrupen") {
  starting.d <- getwd()
  ## source R functions
  setwd(paste(path.2.svn,"/package/quadrupen/R", sep=""))
  source("elastic_net.R")
  source("group_lasso.R")
  source("bootstrap.R")
  source("crossval.R")
  source("pen_fit-class.R")
  ## source C++ libraries
  setwd(paste(path.2.svn,"/package/quadrupen/src", sep=""))
  if (!("enet_ls.so" %in% list.files())) {
    system(paste(path.2.svn,"/bin/Rgcc"," enet_ls",sep="") )
  }
  dyn.load("enet_ls.so")
  if (!("group_l1l2_ls.so" %in% list.files())) {
    system(paste(path.2.svn,"/bin/Rgcc"," group_l1l2_ls",sep="") )
  }
  dyn.load("group_l1l2_ls.so")
  if (!("group_l1linf_ls.so" %in% list.files())) {
    system(paste(path.2.svn,"/bin/Rgcc"," group_l1linf_ls",sep="") )
  }
  dyn.load("group_l1linf_ls.so")
  setwd(starting.d)
}
