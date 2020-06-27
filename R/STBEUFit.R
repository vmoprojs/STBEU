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
STBEUFit<-function(theta,fixed=NULL,coords,times,cc,datos,type_dist=1,
                   maxdist=NULL,maxtime=NULL, winc_s=NULL,winstp_s=NULL,
                   winc_t=NULL,winstp_t=NULL,subs=NULL,weighted=FALSE,
                   local=c(1,1),GPU=NULL,varest = FALSE,optimizer = "Nelder-Mead")             
{
    
  path.parent <- getwd()
  model = cc
  
  if(model == 1)
  {
    param <- c(theta,fixed)
    # cat("Nombres: ",names(param),"\n")
    filtro <- match(names(param) ,c("nugget","mean","scale_t","scale_s","sill"))
    if(!filtro || length(filtro)==0) stop("All parameters (fix and theta) must be named")
  }
  
  if(model == 2)
  {
    param <- c(theta,fixed)
    filtro <- match(names(param) , c("nugget","mean","scale_t","scale_s","sill","power_s","power_t","sep"))
    # cat("Gnei",filtro,"\n")
    if(!filtro || length(filtro)==0) stop("All parameters (fix and theta) must be named")
  }
  if(model == 3)
  {
    param <- c(theta,fixed)
    
    filtro <- match(names(param) , c("nugget","mean","sill",
                                     "power2_s","power_s","power2_t","scale_s",
                                     "scale_t","sep","smooth_t"))
    # cat("Gnei",filtro,"\n")
    if(!filtro || length(filtro)==0) stop("All parameters (fix and theta) must be named")
  }


  nsub_t=floor((( length(times)-winc_t)/(winc_t*winstp_t)+1))
  nsub_s=floor((( (max(coords)-min(coords))-winc_t)/(winc_s*winstp_s)+1))
  
  if(subs==2 && nsub_t<=2)
  {
    cat("Number of temporal blocks is: ",nsub_t,"\n")
    stop("Number of temporal windows must be greater than 2, check winc_t parameter. START values are printed. We recommend 5 blocks or more")
  }
  
  # if(subs==1 && nsub_t<=2)
  # {
  #   cat("Number of spatial blocks is: ",nsub_t,"\n")
  #   stop("Number of spatial windows must be greater than 2, check winc_t parameter. START values are printed. We recommend 5 blocks or more")
  # }
  # 
  # if(subs==1 && nsub_t<=2 && nsub_t<=2)
  # {
  #   cat("Number of spatial blocks is: ",nsub_t,"\n")
  #   stop("Number of spatial windows must be greater than 2, check winc_t parameter. START values are printed. We recommend 5 blocks or more")
  # }
  
  if(!is.null(GPU)) 
  {
   if(GPU == 0)
   {
       if (model ==1) kernel = "DouExp.cl"
       if (model ==2) kernel = "Gneiting.cl"
       if (model ==3) kernel = "Wend.cl"
       path <- system.file("CL", kernel, package = "STBEU")
       path <- gsub(paste("/",kernel,sep=""),"/",path)
       setwd(path)
       .C("create_binary_kernel",  as.integer(GPU),as.character(kernel),  PACKAGE='STBEU',DUP = TRUE, NAOK=TRUE)
   }
   
  if(GPU == 2)
  {
     if (model ==1) kernel = "DouExp.cl"
     if (model ==2) kernel = "Gneiting.cl"
     if (model ==3) kernel = "Wend2.cl"
     path <- system.file("CL", kernel, package = "STBEU")
     path <- gsub(paste("/",kernel,sep=""),"/",path)
     setwd(path)
     .C("create_binary_kernel_GPU",  as.integer(GPU),as.character(kernel),  PACKAGE='STBEU',DUP = TRUE, NAOK=TRUE)
  }
    
    
    if(subs==1 && cc==1) type_sub="SubSamp_space_ocl"
    if(subs==2 && cc==1) type_sub="SubSamp_time_ocl"
    if(subs==3 && cc==1) type_sub="SubSamp_spacetime_ocl"
    
    if(subs==1 && cc==2) type_sub="SubSamp_space_ocl"
    if(subs==2 && cc==2) type_sub="SubSamp_time_ocl"
    if(subs==3 && cc==2) type_sub="SubSamp_spacetime_ocl"
    
    
    if(subs==1 && cc==3) type_sub="SubSamp_space_ocl"
    if(subs==2 && cc==3) {type_sub="SubSamp_time_ocl"}
    if(subs==3 && cc==3) {type_sub="SubSamp_spacetime_ocl"}
    
  }
  if(is.null(GPU))
  {
    if(subs==1) type_sub="SubSamp_space"
    if(subs==2) type_sub="SubSamp_time"
    if(subs==3) type_sub="SubSamp_spacetime"
     kernel = NULL
  }
  
  
  coordx=coords[,1];coordy=coords[,2];ncoords=nrow(coords);ntime=length(times)
  npar=length(theta)
  ## optimization
   print(optimizer)
   if(optimizer=="Nelder-Mead")
   {
     # tiempo = proc.time()
       res=optim(par =theta,fixed=fixed,fn = eucl_st_ocl,coordx=coordx,coordy=coordy,ncoords=ncoords,
       times=times,ntime=ntime,cc=cc,datos=datos,type_dist=type_dist,maxdist= maxdist,
       maxtime=maxtime,winc_s=winc_s,winstp_s=winstp_s,winc_t=winc_t,winstp_t=winstp_t,
       weighted=weighted,type_sub=type_sub,local = local,
       GPU = GPU,control=list(fnscale=1,reltol=1e-14, maxit=10000),kernel = kernel,
       method = optimizer)
      # tiempo = proc.time()-tiempo
      # print(tiempo)
   }else{
       res=optim(par =theta,fixed=fixed,fn = eucl_st_ocl,coordx=coordx,coordy=coordy,ncoords=ncoords,
       times=times,ntime=ntime,cc=cc,datos=datos,type_dist=type_dist,
       maxdist= maxdist,maxtime=maxtime,winc_s=winc_s,winstp_s=winstp_s,
       winc_t=winc_t,winstp_t=winstp_t,
      weighted=weighted,type_sub=type_sub,local = local,
               GPU = GPU,control=list(fnscale=1,factr=1e-10,
               reltol=1e-14, maxit=10000),kernel = kernel,
               method = optimizer)
   }
            
    
# print("TERMINAAAAA")
  # on.exit(rm(res))
  
  names = checkpar(fix=fixed,theta=theta,cc=cc)
  
  
  if(varest==TRUE)
  {
    if(subs==1) type_sub="SubSamp_space"
    if(subs==2) type_sub="SubSamp_time"
    if(subs==3) type_sub="SubSamp_spacetime"

    varval <- varestfun(theta = res$par,fixed = fixed,coordx=coordx,coordy=coordy,ncoords=ncoords,
    times=times,ntime=ntime,cc=cc,datos=datos,type_dist=type_dist,maxdist= maxdist,
    maxtime=maxtime,winc_s=winc_s,winstp_s=winstp_s,winc_t=winc_t,winstp_t=winstp_t,
    weighted=weighted,type_sub=type_sub,local = c(1,1), GPU = NULL)
    # print("I'm BACK")
    
    varvalHess <- varestfunHess(theta = res$par,fixed = fixed,coordx=coordx,coordy=coordy,ncoords=ncoords,
                                times=times,ntime=ntime,cc=cc,datos=datos,type_dist=type_dist,maxdist= maxdist,
                                maxtime=maxtime,winc_s=winc_s,winstp_s=winstp_s,winc_t=winc_t,winstp_t=winstp_t,
                                weighted=weighted,type_sub=type_sub,local = c(1,1), GPU = NULL)
    
    Jaco <- numDeriv::jacobian(varestfunHess,x=res$par,fixed = fixed,coordx=coordx,coordy=coordy,ncoords=ncoords,
                              times=times,ntime=ntime,cc=cc,datos=datos,type_dist=type_dist,maxdist= maxdist,
                              maxtime=maxtime,winc_s=winc_s,winstp_s=winstp_s,winc_t=winc_t,winstp_t=winstp_t,
                              weighted=weighted,type_sub=type_sub,local = c(1,1), GPU = NULL)
    # Jaco <- xpnd(Jaco)
    # Jaco
    
    # cat("\nvarvalHess: ",varvalHess,"\n")
    # cat("\nJacobiano: ",Jaco,"\n")
    # cat("\nvarval: ",solve(varval),"\n")
    
    res$varval <- varval
    res$varcov <- solve(t(Jaco)%*%solve(varval)%*%Jaco)
    res$stderr <- sqrt(diag(res$varcov))
  }
  
  # return(list(res = res, varval = varval,varcov = solve(varval)))
  setwd(path.parent)
  print.STBEUFit(x = res, names = names,GPU,varest)# OJO
  return(res)
}

#######################################################################
#######################################################################
eucl_st_ocl<-function(theta,fixed,coordx,coordy,ncoords,times,ntime,cc,datos,type_dist,maxdist,maxtime,
                      winc_s,winstp_s,winc_t,winstp_t,
                      weighted,type_sub,local,GPU,kernel)
  
{
  setup=setting_param(cc,theta,fixed)
  
  if(cc == 1)
  {
    if(setup$parcor[1] <0 | setup$parcor[2]<0 | setup$nuis[3]<0 ){
      obj = 1e100
      return (obj)}
  }
  
  if(cc ==2)
  {
    # print(setup$parcor)
    if(setup$parcor[1] <0 |setup$parcor[1] >1 | 
       setup$parcor[2]<0 | setup$parcor[2]>1| 
       setup$parcor[3]<=0  | 
       setup$parcor[4]<=0 
       | setup$parcor[5]<0 |setup$parcor[5]>1 | 
       setup$nuis[3]<0){
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
  
  means=double(setup$npar)          ## vector of means
  
  p=.C(as.character(type_sub),as.double(coordx), as.double(coordy), as.double(times),as.integer(ncoords),
       as.integer(ntime),as.integer(cc), as.double(datos),as.integer(type_dist),as.double(maxdist),
       as.double(maxtime),as.integer(setup$npar),vaina = as.double(setup$parcor),as.integer(setup$nparc),
       as.double(setup$nuis),as.integer(setup$nparnuis),as.integer(setup$flagcor),as.integer(setup$flagnuis),
       vv=as.double(vari),as.double(winc_s), as.double(winstp_s),as.double(winc_t), as.double(winstp_t),
       mm=as.double(means),as.integer(weighted),as.integer(local),as.integer(GPU),as.character(kernel),
       PACKAGE='STBEU',DUP = TRUE, NAOK=TRUE)
  # print(str(p))
  # print("OUT++INSIDE")
  # print("Vamo despues de")
  #cat("Parameters:",p$vaina,"\n")
  # cat("FLAGNUIS:",setup$flagnuis,"\n")
  # cat("FLAGCOR:",setup$flagcor,"\n")
  x=p$mm       #means vector
  # cat("A. Means vector:",x,"\n")
  # cat("A. Vari :",p$vv,"\n")
  F.1=xpnd(p$vv) #cov matrix
  Fchol = MatDecomp(F.1,"cholesky")
  if(length(Fchol)==1)
  {
    Fchol = MatDecomp(F.1,"svd")
    inv = MatInv(Fchol,"svd")
    obj= crossprod(crossprod(inv, x), x)
  }
  else
  {
    
    inv=chol2inv(Fchol)  #inverse   
    obj= crossprod(crossprod(inv, x), x)# quadratic form
  }
  # print(obj)
  # if(!is.null(GPU)) gc()
  return(obj)
}



