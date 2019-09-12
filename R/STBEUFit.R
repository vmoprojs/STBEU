####################################################
### Authors:  Moreno Bevilacqua, Víctor Morales Oñate
### Email:  moreno.bevilacqua@uv.cl, victor.morales@uv.cl
### Instituto de Estadistica
### Universidad de Valparaiso
### File name: STBEUFit.r
### Description:
### This file contains a set of procedures
### for the set up of all the package routines.
### Last change: 29/01/2019.
####################################################
STBEUFit<-function(theta,fix,coords,times,cc,data,type_dist,maxdist,maxtime, winc_s,winstp_s,
                       winc_t,winstp_t,subs,weighted,local=c(1,1),GPU=NULL,varest = FALSE)             
{
  # param <- c(theta,fix)
  # print(match(names(param) ,c("nugget","mean","scale_t","scale_s","sill")) )
  # print(match(names(param) , c("nugget","mean","scale_t","scale_s","sill")))
  path.parent <- getwd()
  model = cc
  
  if(model == 1)
  {
    param <- c(theta,fix)
    filtro <- match(names(param) ,c("nugget","mean","scale_t","scale_s","sill"))
    if(!filtro || length(filtro)==0) stop("All parameters (fix and theta) must be named")
  }
  
  if(model == 2)
  {
    param <- c(theta,fix)
    filtro <- match(names(param) , c("nugget","mean","scale_t","scale_s","sill","power_s","power_t","sep"))
    # cat("Gnei",filtro,"\n")
    if(!filtro || length(filtro)==0) stop("All parameters (fix and theta) must be named")
  }
  if(model == 3)
  {
    param <- c(theta,fix)
    filtro <- match(names(param) , c("nugget","mean","sill",
                                     "power2_s","power_s","power2_t","scale_s",
                                     "scale_t","sep","smooth_t"))
    # cat("Gnei",filtro,"\n")
    if(!filtro || length(filtro)==0) stop("All parameters (fix and theta) must be named")
  }


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
    if (model ==3) kernel = "Wend.cl"
    path <- system.file("CL", kernel, package = "STBEU")
    path <- gsub(paste("/",kernel,sep=""),"/",path);setwd(path)
    # print(path)
    # fname <- paste(fname,"_OCL",sep="")
    .C("create_binary_kernel",  as.integer(GPU),as.character(kernel),  PACKAGE='STBEU',DUP = TRUE, NAOK=TRUE)
    
    if(subs==1 && cc==1) type_sub="SubSamp_space_DE_ocl"
    if(subs==2 && cc==1) type_sub="SubSamp_time_DE_ocl"
    if(subs==3 && cc==1) type_sub="SubSamp_spacetime_DE_ocl"
    
    if(subs==1 && cc==2) type_sub="SubSamp_space_GN_ocl"
    if(subs==2 && cc==2) type_sub="SubSamp_time_GN_ocl"
    if(subs==3 && cc==2) type_sub="SubSamp_spacetime_GN_ocl"
    
    
    if(subs==1 && cc==3) type_sub="SubSamp_space_WE_ocl"
    if(subs==2 && cc==3) {type_sub="SubSamp_time_WE_ocl"}
    if(subs==3 && cc==3) {type_sub="SubSamp_spacetime_WE_ocl"}
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
  # if(cc ==3)
  # {
  #   print(setup$parcor)
  #   cat("\nsetup$parcor[1]<(2.5+2*setup$parcor[7])",setup$parcor[1]<(2.5+2*setup$parcor[7]),"\n\n\n",
  #       "setup$parcor[2]!=2",setup$parcor[2]!=2,"\n\n\n",
  #       "setup$parcor[2]<(3.5+setup$parcor[7])",setup$parcor[3]<(3.5+setup$parcor[7]),"\n\n\n")
  #   if(setup$parcor[4]<=0|setup$parcor[5]<=0 |
  #      setup$parcor[1]<(2.5+2*setup$parcor[7]) | setup$parcor[2]!=2 |
  #      setup$parcor[6]<0 |setup$parcor[6] >1|
  #      setup$parcor[2]<(3.5+setup$parcor[7]) | setup$nuis[3]<0){
  #     obj = 1e100
  #     return (obj)}
  # }
  
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
  
  means=double(setup$npar)                               ## vector of means
  # cat("MEANS: ",means,"\n",
  #     "parcor: ",setup$parcor,"\n",
  #     "nuis: ",setup$nuis,"\n",
  #     "flagcor: ",setup$flagcor,"\n",
  #     "flagnuis: ",setup$flagnuis,"\n")
  # cat("parcor:",setup$parcor,"\n")
  # cat("parnuis:",setup$parnuis,"\n")
  # cat("flagcor:",setup$flagcor,"\n")
  # cat("flagnuis:",setup$flagnuis,"\n")
  p=.C(type_sub,as.double(coordx), as.double(coordy), as.double(times),as.integer(ncoords),
       as.integer(ntime),as.integer(cc), as.double(data),as.integer(type_dist),as.double(maxdist),
       as.double(maxtime),as.integer(setup$npar),as.double(setup$parcor),as.integer(setup$nparc),
       as.double(setup$nuis),as.integer(setup$nparnuis),as.integer(setup$flagcor),as.integer(setup$flagnuis),
       vv=as.double(vari),as.double(winc_s), as.double(winstp_s),as.double(winc_t), as.double(winstp_t),
       mm=as.double(means),as.integer(weighted),as.integer(local),as.integer(GPU))
  
  x=p$mm       #means vector
  # cat("Ax:",x,"\n")
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
  
  # if(cc ==3)
  # {
  #   if(setup$parcor[4] <0 |setup$parcor[5] <0 |
  #      setup$parcor[6] <0 |setup$parcor[6] >1 |
  #      setup$parcor[7]<0|
  #     setup$nuis[3]<0){
  #     obj = 1e100
  #     return (obj)}
  # }
  
  # if(cc ==3)
  # {
  #   if(setup$parcor[4]<=0||setup$parcor[5]<=0 || 
  #      setup$parcor[1]<(2.5+2*setup$parcor[7]) || setup$parcor[2]>2|| 
  #      setup$parcor[2]<0||setup$parcor[6]<0 ||setup$parcor[6] >1|| 
  #      setup$parcor[2]<(3.5+setup$parcor[7])||setup$parcor[7]<0 || setup$nuis[3]<0){
  #     obj = 1e100
  #     return (obj)}
  # }
  
  if(cc ==3)
  {
    if(setup$parcor[4]<=0|setup$parcor[5]<=0 |
       setup$parcor[1]<(2.5+2*setup$parcor[7]) | setup$parcor[2]!=2 |
       setup$parcor[6]<0 |setup$parcor[6] >1|
       setup$parcor[3]<(3.5+setup$parcor[7]) | setup$nuis[3]<0){
      # print("cc=3!!!!")
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

