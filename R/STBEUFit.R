####################################################
### Authors:  Moreno Bevilacqua, Víctor Morales Oñate
### Email:  moreno.bevilacqua@uv.cl, victor.morales@uv.cl
### Instituto de Estadistica
### Universidad de Valparaiso
### File name: STBEUFit.r
### Description:
### This file contains a set of procedures
### for the set up of all the package routines.
### Last change: 13/09/2017.
####################################################
STBEUFit<-function(theta,fix,coords,times,cc,data,type_dist,maxdist,maxtime, winc_s,winstp_s,
                       winc_t,winstp_t,subs,weighted,local=c(1,1),GPU=NULL,varest = FALSE)             
{
  path.parent <- getwd()
  model = cc

  nsub_t=floor((( length(times)-winc_t)/(winc_t*winstp_t)+1))
  if(subs==2 && nsub_t<=2)
  {
    cat("Number of temporal blocks is: ",nsub_t,"\n")
    stop("Number of temporal windows must be greater than 2, check winc_t parameter. START values are printed. We recommend 5 blocks or more")
  }
  if(!is.null(GPU)) 
  {
    if (model ==1) kernel = "DouExp.cl"
    if (model ==2) kernel = "Gneiting.cl"
    path <- system.file("CL", kernel, package = "STBEU")
    path <- gsub(paste("/",kernel,sep=""),"/",path);setwd(path)
    # fname <- paste(fname,"_OCL",sep="")
    .C("create_binary_kernel",  as.integer(GPU),as.character(kernel),  PACKAGE='STBEU',DUP = TRUE, NAOK=TRUE)
    
    if(subs==1 && cc==1) type_sub="SubSamp_space_DE_ocl"
    if(subs==2 && cc==1) type_sub="SubSamp_time_DE_ocl"
    if(subs==3 && cc==1) type_sub="SubSamp_spacetime_DE_ocl"
    
    if(subs==1 && cc==2) type_sub="SubSamp_space_GN_ocl"
    if(subs==2 && cc==2) type_sub="SubSamp_time_GN_ocl"
    if(subs==3 && cc==2) type_sub="SubSamp_spacetime_GN_ocl"
  }else
  {
    if(subs==1) type_sub="SubSamp_space"
    if(subs==2) type_sub="SubSamp_time"
    if(subs==3) type_sub="SubSamp_spacetime"
  }
  coordx=coords[,1];coordy=coords[,2];ncoords=nrow(coords);ntime=length(times)
  npar=length(theta)
  ## optimization
  res=optim(par =theta,fix=fix,fn = eucl_st_ocl,coordx=coordx,coordy=coordy,ncoords=ncoords,
            times=times,ntime=ntime,cc=cc,data=data,type_dist=type_dist,maxdist= maxdist,
            maxtime=maxtime,winc_s=winc_s,winstp_s=winstp_s,winc_t=winc_t,winstp_t=winstp_t,
            weighted=weighted,type_sub=type_sub,local = local, GPU = GPU,control=list(maxit=1000))
  setwd(path.parent)
  names = checkpar(fix=fix,theta=theta,cc=cc)
  
  
  if(varest==TRUE)
  {
    if(subs==1) type_sub="SubSamp_space"
    if(subs==2) type_sub="SubSamp_time"
    if(subs==3) type_sub="SubSamp_spacetime"
    
    varval <- varestfun(theta = res$par,fix = fix,coordx=coordx,coordy=coordy,ncoords=ncoords,
    times=times,ntime=ntime,cc=cc,data=data,type_dist=type_dist,maxdist= maxdist,
    maxtime=maxtime,winc_s=winc_s,winstp_s=winstp_s,winc_t=winc_t,winstp_t=winstp_t,
    weighted=weighted,type_sub=type_sub,local = c(1,1), GPU = NULL)
    
    res$varval <- varval
    res$varcov <- solve(varval)
    res$stderr <- sqrt(diag(solve(varval)))
  }
  
  # return(list(res = res, varval = varval,varcov = solve(varval)))
  print.STBEUFit(x = res, names = names,GPU,varest)# OJO
  return(res)
}

#######################################################################
#######################################################################
eucl_st_ocl<-function(theta,fix,coordx,coordy,ncoords,times,ntime,cc,data,type_dist,maxdist,maxtime,
                      winc_s,winstp_s,winc_t,winstp_t,weighted,type_sub,local,GPU)
  
{
  # print(theta)
  setup=setting_param(cc,theta,fix)
  # cat("length nuis",length(setup$nuis),"length parcor",length(setup$parcor),"\n")
  # print(setup$parcor)
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
  
  vari=double(setup$npar*0.5*(setup$npar-1)+setup$npar) ## vector of upper  triangular cov matrix (with diagonal)
  
  means=double(setup$npar)                               ## vector of means
  
  # if(setup$nuis[1]<0 | setup$nuis[2]<0 | setup$nuis[3]<0) print(setup$nuis)
  # print(type_sub)
  # print(system.time(.C(type_sub,as.double(coordx), as.double(coordy), as.double(times),as.integer(ncoords),
  #                  as.integer(ntime),as.integer(cc), as.double(data),as.integer(type_dist),as.double(maxdist),
  #                  as.double(maxtime),as.integer(setup$npar),as.double(setup$parcor),as.integer(setup$nparc),
  #                  as.double(setup$nuis),as.integer(setup$nparnuis),as.integer(setup$flagcor),as.integer(setup$flagnuis),
  #                  vv=as.double(vari),as.double(winc_s), as.double(winstp_s),as.double(winc_t), as.double(winstp_t),
  #                  mm=as.double(means),as.integer(weighted),as.integer(local),as.integer(GPU))))
  p=.C(type_sub,as.double(coordx), as.double(coordy), as.double(times),as.integer(ncoords),
       as.integer(ntime),as.integer(cc), as.double(data),as.integer(type_dist),as.double(maxdist),
       as.double(maxtime),as.integer(setup$npar),as.double(setup$parcor),as.integer(setup$nparc),
       as.double(setup$nuis),as.integer(setup$nparnuis),as.integer(setup$flagcor),as.integer(setup$flagnuis),
       vv=as.double(vari),as.double(winc_s), as.double(winstp_s),as.double(winc_t), as.double(winstp_t),
       mm=as.double(means),as.integer(weighted),as.integer(local),as.integer(GPU))
  
  x=p$mm       #means vector
  F=xpnd(p$vv) #cov matrix
  Fchol = MatDecomp(F,"cholesky")
  if(length(Fchol)==1)
  {
    Fchol = MatDecomp(F,"svd")
    inv = MatInv(Fchol,"svd")
    obj= crossprod(crossprod(inv, x), x)
  }
  else
  {
    
    inv=chol2inv(Fchol)  #inverse   
    obj= crossprod(crossprod(inv, x), x)# quadratic form
  }
  return(obj)
}



varestfun<-function(theta,fix,coordx,coordy,ncoords,times,ntime,cc,data,type_dist,maxdist,maxtime,
                      winc_s,winstp_s,winc_t,winstp_t,weighted,type_sub,local,GPU)
  
{
  # print(theta)
  setup=setting_param(cc,theta,fix)
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
  
  vari=double(setup$npar*0.5*(setup$npar-1)+setup$npar) ## vector of upper  triangular cov matrix (with diagonal)
  
  means=double(setup$npar)                               ## vector of means
  
  p=.C(type_sub,as.double(coordx), as.double(coordy), as.double(times),as.integer(ncoords),
       as.integer(ntime),as.integer(cc), as.double(data),as.integer(type_dist),as.double(maxdist),
       as.double(maxtime),as.integer(setup$npar),as.double(setup$parcor),as.integer(setup$nparc),
       as.double(setup$nuis),as.integer(setup$nparnuis),as.integer(setup$flagcor),as.integer(setup$flagnuis),
       vv=as.double(vari),as.double(winc_s), as.double(winstp_s),as.double(winc_t), as.double(winstp_t),
       mm=as.double(means),as.integer(weighted),as.integer(local),as.integer(GPU))
  
  x=p$mm       #means vector
  varval=xpnd(p$vv) #cov matrix
  return(varval)
}

