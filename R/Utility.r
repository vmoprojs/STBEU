####################################################
### Authors:  Moreno Bevilacqua, Víctor Morales Oñate
### Email:  moreno.bevilacqua@uv.cl, victor.morales@uv.cl
### Instituto de Estadistica
### Universidad de Valparaiso
### File name: Utility.r
### Description:
### This file contains a set of procedures
### for the set up of all the package routines.
### Last change: 29/01/2019.
####################################################

### decomposition of a square  matrix
MatDecomp<-function(mtx,method)    {
  # print(mtx)
  if(method=="cholesky")  {
    mat.decomp <- try(chol(mtx), silent=TRUE)
    if (inherits(mat.decomp , "try-error")) return (FALSE)
  }
  if(method=="svd")      {
    mat.decomp <- svd(mtx, nv = 0)
    cov.logdeth <- try(sum(log(sqrt(mat.decomp$d))), silent=TRUE)
    if (inherits(cov.logdeth, "try-error"))  return (FALSE)
  }
  return(mat.decomp)
}
### inverse a square matrix given a decomposition
MatInv<-function(mat.decomp,method)    {
  
  if(method=="cholesky")  varcov.inv <- chol2inv(mat.decomp)
  if(method=="svd")       
  { 
    tol = sqrt(.Machine$double.eps)
    e <- mat.decomp$d;e[e > tol] <- 1/e[e > tol] 
    varcov.inv<-mat.decomp$u %*% diag(e,nrow=length(e)) %*% t(mat.decomp$u) 
  } 
  return(varcov.inv)
} 

####################################################################################################################
####################################################################################################################


DevOpenCL <- function()
{
    .C("DevOpenCL",PACKAGE='STBEU',DUP = TRUE, NAOK=TRUE)
}

CreateBinary <- function(dev,fname)
{
  # library(STBEU)
  # fname= "Wend"
  # dev = 0
  fname <- paste(fname,".cl",sep="")
  # cat("fname de Composit.r: ",fname,"\n")
  
  path.parent <- getwd()
  
  path <- system.file("CL", fname, package = "STBEU")
  path <- gsub(paste("/",fname,sep = ""),"/",path)
  # .C("create_binary_kernel",  as.integer(GPU),as.character(fname),  PACKAGE='GeoModels',DUP = TRUE, NAOK=TRUE)
  setwd(path)
  # fname = gsub(".cl","",fname)
  .C("create_binary_kernel",as.integer(dev),as.character(fname),PACKAGE='STBEU',DUP = TRUE, NAOK=TRUE)
  setwd(path.parent)
}

#######################################################################
setting_param<-function(cc,theta,fix)
{
  param <- c(theta,fix)
  # if(!(all(names(param) %in% c("nugget","mean","scale_t","scale_s","sill")))) stop("All parameters (fix and theta) must be named")
  if(cc==1){      ### double exponential model
    namespar = c("nugget","mean","scale_t","scale_s","sill")
    npar=length(theta)
    parcor=c(param["scale_s"],param["scale_t"])      #corr parameters
    nuis=c(param["mean"],param["nugget"],param["sill"])   #mean, nugget, sill
    
    flagcor <- c(0,0)
    
    if(any(names(theta)=="scale_s"))
    {
      flagcor[1] <- 1
    }
    if(any(names(fix)=="scale_s"))
    {
      flagcor[1] <- 0
    }
    if(any(names(theta)=="scale_t"))
    {
      flagcor[2] <- 1
    }
    if(any(names(fix)=="scale_t"))
    {
      flagcor[2] <- 0
    }
    # print(flagcor)
    flagnuis <- c(0,0,0)
    if( any(names(fix)=="mean") )
    {
      flagnuis[1] <- 0
    }
    if( any(names(theta)=="mean") )
    {
      flagnuis[1] <- 1
    }
    if( any(names(fix)=="nugget") )
    {
      flagnuis[2] <- 0
    }
    if( any(names(theta)=="nugget") )
    {
      flagnuis[2] <- 1
    }
    if( any(names(fix)=="sill") )
    {
      flagnuis[3] <- 0
    }
    if( any(names(theta)=="sill") )
    {
      flagnuis[3] <- 1
    }
    # print(flagnuis)
    # flagnuis=c(0,0,1)
    nparc=sum(flagcor)  
    nparnuis=sum(flagnuis)
  }
  if(cc==2){       ### gneiting model
    npar=length(theta)
    parcor=c(param["power_s"],param["power_t"],param["scale_s"],param["scale_t"],param["sep"])     #corr parameters
    nuis=c(param["mean"],param["nugget"],param["sill"])   #mean, nugget, sill
    
    flagcor=c(0,0,0,0,0) 
    
    if(any(names(theta)=="power_s"))
    {
      flagcor[1] <- 1
    }
    if(any(names(fix)=="power_s"))
    {
      flagcor[1] <- 0
    }
    if(any(names(theta)=="power_t"))
    {
      flagcor[2] <- 1
    }
    if(any(names(fix)=="power_t"))
    {
      flagcor[2] <- 0
    }
    if(any(names(theta)=="scale_s"))
    {
      flagcor[3] <- 1
    }
    if(any(names(fix)=="scale_s"))
    {
      flagcor[3] <- 0
    }
    if(any(names(theta)=="scale_t"))
    {
      flagcor[4] <- 1
    }
    if(any(names(fix)=="scale_t"))
    {
      flagcor[4] <- 0
    }
    if(any(names(theta)=="sep"))
    {
      flagcor[5] <- 1
    }
    if(any(names(fix)=="sep"))
    {
      flagcor[5] <- 0
    }
    # print(flagcor)
    
    
    flagnuis=c(0,0,0) 
    
    if( any(names(fix)=="mean") )
    {
      flagnuis[1] <- 0
    }
    if( any(names(theta)=="mean") )
    {
      flagnuis[1] <- 1
    }
    if( any(names(fix)=="nugget") )
    {
      flagnuis[2] <- 0
    }
    if( any(names(theta)=="nugget") )
    {
      flagnuis[2] <- 1
    }
    if( any(names(fix)=="sill") )
    {
      flagnuis[3] <- 0
    }
    if( any(names(theta)=="sill") )
    {
      flagnuis[3] <- 1
    }
    
    # flagcor=c(0,0,1,1,0) 
    # flagnuis=c(0,0,1) 
    nparc=sum(flagcor)  
    nparnuis=sum(flagnuis)
  }
  if(cc==3){       ### Wen_time
    npar=length(theta)
    # power2_s,power_s,power2_t,scale_s,scale_t,sep,smooth_t
    parcor=c(param["power2_s"],param["power_s"],param["power2_t"],
             param["scale_s"],param["scale_t"],param["sep"],param["smooth_t"]) #corr parameters
    nuis=c(param["mean"],param["nugget"],param["sill"])   #mean, nugget, sill
    
    flagcor=c(0,0,0,0,0,0,0) 
    
    if(any(names(theta)=="power2_s"))
    {
      flagcor[1] <- 1
    }
    if(any(names(fix)=="power2_s"))
    {
      flagcor[1] <- 0
    }
    if(any(names(theta)=="power_s"))
    {
      flagcor[2] <- 1
    }
    if(any(names(fix)=="power_s"))
    {
      flagcor[2] <- 0
    }
    if(any(names(theta)=="power2_t"))
    {
      flagcor[3] <- 1
    }
    if(any(names(fix)=="power2_t"))
    {
      flagcor[3] <- 0
    }
    if(any(names(theta)=="scale_s"))
    {
      flagcor[4] <- 1
    }
    if(any(names(fix)=="scale_s"))
    {
      flagcor[4] <- 0
    }
    if(any(names(theta)=="scale_t"))
    {
      flagcor[5] <- 1
    }
    if(any(names(fix)=="scale_t"))
    {
      flagcor[5] <- 0
    }
    if(any(names(theta)=="sep"))
    {
      flagcor[6] <- 1
    }
    if(any(names(fix)=="sep"))
    {
      flagcor[6] <- 0
    }
    if(any(names(theta)=="smooth_t"))
    {
      flagcor[7] <- 1
    }
    if(any(names(fix)=="smooth_t"))
    {
      flagcor[7] <- 0
    }
    # cat("flagcor: ",flagcor,"\n")
    
    
    flagnuis=c(0,0,0) 
    
    if( any(names(fix)=="mean") )
    {
      flagnuis[1] <- 0
    }
    if( any(names(theta)=="mean") )
    {
      flagnuis[1] <- 1
    }
    if( any(names(fix)=="nugget") )
    {
      flagnuis[2] <- 0
    }
    if( any(names(theta)=="nugget") )
    {
      flagnuis[2] <- 1
    }
    if( any(names(fix)=="sill") )
    {
      flagnuis[3] <- 0
    }
    if( any(names(theta)=="sill") )
    {
      flagnuis[3] <- 1
    }
    
    # flagcor=c(0,0,1,1,0) 
    # flagnuis=c(0,0,1) 
    nparc=sum(flagcor)  
    nparnuis=sum(flagnuis)
  }
  setup=list(flagcor=flagcor,flagnuis=flagnuis,npar=npar,parcor=parcor,nuis=nuis,nparc=nparc
             ,nparnuis=nparnuis)
  return(setup)
}



checkpar <- function(fix,theta,cc)
{
  if(cc==1)
  {
    namespar = c("nugget","mean","scale_t","scale_s","sill")
    if((length(fix)+length(theta))!=5) stop("Number of parameters don't match the declared model")
    if(is.null(names(fix)) || is.null(names(theta))) 
    {
      warning("fix and theta should be named. If not, take care about the ouput order") 
    }else
    {
      namesfix = names(fix);namestheta = names(theta)
    }
    res = names(theta)
  }
  if(cc==2)
  {
    namespar = c("nugget","mean","scale_t","scale_s","sill",
                 "power_s","power_t","sep")
    if((length(fix)+length(theta))!=8) stop("Number of parameters don't match the declared model")
    if(is.null(names(fix)) || is.null(names(theta))) 
    {
      warning("fix and theta should be named. If not, take care about the ouput order") 
    }else
    {
      namesfix = names(fix);namestheta = names(theta)
    }
    res = names(theta)
  }
  
  if(cc==3)
  {
    namespar = c("nugget","mean","sill",
                 "power2_s","power_s","power2_t","scale_s",
                 "scale_t","sep","smooth_t")
    if((length(fix)+length(theta))!=10) stop("Number of parameters don't match the declared model")
    if(is.null(names(fix)) || is.null(names(theta))) 
    {
      warning("fix and theta should be named. If not, take care about the ouput order") 
    }else
    {
      namesfix = names(fix);namestheta = names(theta)
    }
    res = names(theta)
  }
  return(res)
}
# checkpar(fix,theta)

EvalRho <- function(cc,h,u,parEval)
{
  rr = 0
  sol = .C("CorFct_call",as.integer(cc),as.double(h),as.double(u),
           as.double(parEval),res = as.double(rr),  
           PACKAGE='STBEU',DUP = TRUE, NAOK=TRUE)
  return(sol$res)
}


print.STBEUFit <- function(x,names,GPU,varest, digits = max(3, getOption("digits") - 3), ...)
{
  cat('\n=================================================================')
  if(x$convergence==0) cat('\nSuccessful convergence') else cat('\nConvergence may have failed')
  if(is.null(GPU)) cat('\nCPU device for computation') else cat('\nOpenCL device for computation')
  if(is.null(names)) {cat('\n Parameters (order of theta)\n',x$par)} else
  {
    res = x$par;names(res) = names
    cat('\nParameters:\n')
    print.default(res)
    if(varest)
    {
      cat('\nStandard errors:\n')
      stderr = x$stderr;names(stderr) = names
      print.default(stderr)
      cat('\nVariance-covariance matrix of the estimates:\n')
      varcov = x$varcov;colnames(varcov) = names;rownames(varcov) = names
      print.default(varcov)
    }
  }
  cat('\n=================================================================\n')
}


varestfun<-function(theta,fixed=NULL,coordx,coordy,ncoords,times,ntime,cc,datos,type_dist=1,maxdist=NULL,maxtime=NULL,
                    winc_s=NULL,winstp_s=NULL,winc_t=NULL,winstp_t=NULL,
                    weighted=FALSE,type_sub=NULL,local=c(1,1),GPU=NULL)
  
{
  # cat("\nParametros",theta,"\n")
  setup=setting_param(cc,theta,fixed)
  if(cc == 1)
  {
    if(setup$parcor[1] <0 | setup$parcor[2]<0 | setup$nuis[3]<0 ){
      obj = 1e100
      return (obj)}
  }
  if(cc ==2)
  {
    if(setup$parcor[1] <0 |setup$parcor[1] >1 | setup$parcor[2]<0 | setup$parcor[2]>1| setup$parcor[3]<=0  | setup$parcor[4]<=0 
       | setup$parcor[5]<0 |setup$parcor[5]>1 | setup$nuis[3]<0){
      obj = 1e100
      return (obj)}
  }
  if(cc ==3)
  {
    # setup$parcor
    if(setup$parcor[4]<=0| #scale_s
       setup$parcor[5]<=0| #scale_t
       setup$nuis[3]<0| #sill
       setup$parcor[1]<0| #power2_s R_power
       setup$parcor[2]<0| #power_s R_power_s
       setup$parcor[6]<0| #sep
       setup$parcor[3]<0| #power2_t R_power_t
       setup$parcor[7]<0| #smooth_t
       setup$parcor[7]>4| #smooth_t
       setup$parcor[1]<(2.5+2*setup$parcor[7])
       |setup$parcor[3]<(3.5 + setup$parcor[7])
    ){
      obj = 1e100
      return (obj)}
  }
  
  vari=double(setup$npar*0.5*(setup$npar-1)+setup$npar) ## vector of upper  triangular cov matrix (with diagonal)
  grad = double(setup$npar)
  # cat("B. flagcor:",setup$flagcor,"\n")
  # cat("B. flagnuis:",setup$flagnuis,"\n")
  means=double(setup$npar)                               ## vector of means
  # print(paste(as.character(type_sub),"1",sep = ""))
  p=.C(as.character(type_sub),as.double(coordx), as.double(coordy), as.double(times),as.integer(ncoords),
       as.integer(ntime),as.integer(cc), as.double(datos),as.integer(type_dist),as.double(maxdist),
       as.double(maxtime),as.integer(setup$npar),as.double(setup$parcor),as.integer(setup$nparc),
       as.double(setup$nuis),as.integer(setup$nparnuis),as.integer(setup$flagcor),as.integer(setup$flagnuis),
       vv=as.double(vari),as.double(winc_s), as.double(winstp_s),as.double(winc_t), as.double(winstp_t),
       mm=as.double(means),as.integer(weighted),as.integer(local),as.integer(GPU),
       PACKAGE='STBEU',DUP = TRUE, NAOK=TRUE)
  
  x=p$mm       #means vector
  grad = p$grad
  # cat("B. GRAD:",grad,"\n")
  
  # cat("B. Means vector:",x,"\n")
  # cat("B. Vari :",xpnd(p$vv),"\n")
  
  # varval=solve(diag(grad)%*%solve(xpnd(p$vv))%*%diag(grad)) #cov matrix
  varval=xpnd(p$vv)#cov matrix
  
  return(varval)
}





varestfunHess<-function(theta,fixed=NULL,coordx,coordy,ncoords,times,ntime,cc,datos,type_dist=1,maxdist=NULL,maxtime=NULL,
                        winc_s=NULL,winstp_s=NULL,winc_t=NULL,winstp_t=NULL,
                        weighted=FALSE,type_sub=NULL,local=c(1,1),GPU=NULL)

{
  # cat("\nParametros",theta,"\n")
  setup=setting_param(cc,theta,fixed)
  if(cc == 1)
  {
    if(setup$parcor[1] <0 | setup$parcor[2]<0 | setup$nuis[3]<0 ){
      obj = 1e100
      return (obj)}
  }
  if(cc ==2)
  {
    if(setup$parcor[1] <0 |setup$parcor[1] >1 | setup$parcor[2]<0 | setup$parcor[2]>1| setup$parcor[3]<=0  | setup$parcor[4]<=0
       | setup$parcor[5]<0 |setup$parcor[5]>1 | setup$nuis[3]<0){
      obj = 1e100
      return (obj)}
  }
  if(cc ==3)
  {
    # setup$parcor
    if(setup$parcor[4]<=0| #scale_s
       setup$parcor[5]<=0| #scale_t
       setup$nuis[3]<0| #sill
       setup$parcor[1]<0| #power2_s R_power
       setup$parcor[2]<0| #power_s R_power_s
       setup$parcor[6]<0| #sep
       setup$parcor[3]<0| #power2_t R_power_t
       setup$parcor[7]<0| #smooth_t
       setup$parcor[7]>4| #smooth_t
       setup$parcor[1]<(2.5+2*setup$parcor[7])
       |setup$parcor[3]<(3.5 + setup$parcor[7])
    ){
      obj = 1e100
      return (obj)}
  }

  vari=double(setup$npar*0.5*(setup$npar-1)+setup$npar) ## vector of upper  triangular cov matrix (with diagonal)
  grad = double(setup$npar)
  means=double(setup$npar)                               ## vector of means
  # print(paste(as.character(type_sub),"1",sep = ""))
  p=.C(as.character(type_sub),as.double(coordx), as.double(coordy), as.double(times),as.integer(ncoords),
       as.integer(ntime),as.integer(cc), as.double(datos),as.integer(type_dist),as.double(maxdist),
       as.double(maxtime),as.integer(setup$npar),as.double(setup$parcor),as.integer(setup$nparc),
       as.double(setup$nuis),as.integer(setup$nparnuis),as.integer(setup$flagcor),as.integer(setup$flagnuis),
       vv=as.double(vari),as.double(winc_s), as.double(winstp_s),as.double(winc_t), as.double(winstp_t),
       mm=as.double(means),as.integer(weighted),as.integer(local),as.integer(GPU),
       PACKAGE='STBEU',DUP = TRUE, NAOK=TRUE)

  x=p$mm       #means vector
  return(x)
}



