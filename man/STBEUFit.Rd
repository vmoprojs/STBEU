\name{STBEUFit}
\alias{STBEUFit}
\encoding{UTF-8}
\title{Blowise Eucliden Likelihood Fitting of Gaussian Random Fields.}
\description{
  Blowise Eucliden Likelihood Fitting of Gaussian Random Fields.
  The function returns the model parameters' estimates and the estimates'
  variances. Moreover the function allows to fix any of the parameters in the optimization.
}
\usage{
STBEUFit(theta,fixed=NULL,coords,times,cc,datos,type_dist=1,maxdist=NULL,maxtime=NULL, winc_s=NULL,winstp_s=NULL,
                       winc_t=NULL,winstp_t=NULL,subs=NULL,weighted=FALSE,local=c(1,1),GPU=NULL,varest = FALSE)
}

\arguments{
\item{theta}{A vector of starting parametes for estimation}
\item{fixed}{An optional named vector giving the values of the parameters that
    will be considered as known values. The default is \code{NULL}}
\item{coords}{A numeric (\eqn{d \times 2}{d x 2})-matrix (where
    \code{d} is the number of spatial sites).}
\item{times}{A numeric vector assigning 1-dimension of
    temporal coordinates.}    
\item{cc}{Numerical, the name of a correlation model. Values \eqn{(1,2)} are allowed for Double Exponential and Gneting respectively. See the Section \bold{Details}.}    
\item{datos}{A \eqn{t}{d}-matrix (a single spatial realisation) where \code{d} is the number of spatial sites and \code{t} is the number of temporal coordinates. For the description see the Section \bold{Details}.}
\item{type_dist}{Numerical; \code{1} is Euclinean and \code{2} is Geodesic. The default is \code{1}, the euclidean distance. See the Section \bold{Details}.}
\item{maxdist}{Numeric; an optional positive value indicating the maximum
    spatial distance considered in the composite likelihood
    computation. See the Section \bold{Details} for more information.}
\item{maxtime}{Numeric; an optional positive value indicating the maximum
    temporal separation considered in the euclidean likelihood
    computation (see \bold{Details}).}
\item{winc_s}{Numeric; a positive value for computing the spatial  sub-window
      in the sub-sampling procedure. See \bold{Details} for more information.}
\item{winstp_s}{Numeric; a value in \eqn{(0,1]} for defining the the proportion of 
      overlapping in the spatial  sub-sampling procedure. The case \eqn{1} 
    correspond to no overlapping. See \bold{Details} for more information.}
\item{winc_t}{Numeric; a positive value for computing the temporal sub-window
         in the sub-sampling procedure. See \bold{Details} for more information.}
\item{winstp_t}{Numeric; a value in \eqn{(0,1]} for defining the the proportion of 
      overlapping in the temporal  sub-sampling procedure. The case \eqn{1} 
    correspond to no overlapping. See \bold{Details} for more information.}
\item{subs}{Numeric; a value in \eqn{(1,2,3)} for defining the type of  sub-sampling procedure. Cases are space, time and spacetime respectively. See \bold{Details} for more information.}
\item{weighted}{Logical; if \code{TRUE} the likelihood objects are
    weighted, see the Section \bold{Details}. If \code{FALSE} (the
    default) the euclidean likelihood is not weighted.}
\item{local}{Numeric; number of local work-items of the OpenCL setup. Default is \eqn{c(1,1)}}
\item{GPU}{Numeric; if \code{NULL} (the default) 
      no GPU computation is performed.}
\item{varest}{Logical; if \code{FALSE} (the default) 
      no varcov matrix is returned.}
}

\details{
  The optimization method is \code{Nelder-mead}.

  The \code{cc} corrmodel parameter allows to select a specific correlation
  function for the RF. Options are: \eqn{1} for Double Exponential and \eqn{2} for Gneting.


  The \code{distance} parameter allows to consider differents kinds of spatial distances.
  The settings alternatives are:
  \enumerate{
    \item \code{Eucl}, the euclidean distance (default value);
    \item \code{Geod}, the geodesic distance;
  }
  
  The \code{subs} parameter represents the subsampling
  configurations. The settings alternatives are:
  \enumerate{
    \item \code{Spatial}, preferred option when the number of space sites is greater than temporal coordinates;
    \item \code{Time}, preferred option when the number of temporal coordinates is greater than space sites;
    \item \code{Spacetime}, preferred option when the number of space and temporal are balaced.
  }

   All the nuisance and covariance parameters must be specified
   by the user using  the  \code{start}  and  the \code{fixed} parameter.
   Specifically:

   The \code{start} parameter allows to specify (as starting values for the optimization)  the parameters  to be estimated.
   The \code{fixed} parameter allows to fix  some of the parameters.

    The \code{maxdist} parameter set the maximum
  spatial distance below which pairs of sites with inferior distances
  are considered in the euclidean-likelihood. This can be
  lower of the  maximum spatial distance. \bold{Note} that
  this corresponds to use a weighted composite-likelihood with binary
  weights. Pairs with distance less than \code{maxdist} have weight 1
  and are included in the likelihood computation, instead those with
  greater distance have weight 0 and then excluded.
  The default is \code{NULL}, in this case the effective maximum spatial distance
  between sites is considered.


  The same arguments of \code{maxdist} are valid for \code{maxtime} but
  here the weigthed composite-likelihood regards the case of
  spatial-temporal field. At the moment is
  implemented only for Gaussian RFs. The default
  is \code{NULL}, in this case the effective maximum temporal lag
  between pairs of observations is considered.


  The \code{weighted} parameter specifies if the likelihoods forming the
  composite-likelihood must be weighted. If \code{TRUE} the weights are
  selected by opportune procedures that improve the efficient of the
  maximum composite-likelihood estimator (not implemented yet). If
  \code{FALSE} the efficient improvement procedure is not used.

  For computing the standard errors by the sub-sampling procedure,
  \code{winconst} and \code{winstp} parameters represent respectively a positive constant used to
  determine the sub-window size and the the step with which the
  sub-window moves.

  In the spatio-temporal case the subsampling is meant only in time as
  described by Li et al. (2007). Thus, \code{winconst} represents
  the lenght of the temporal sub-window. By default the size of the
  sub-window is computed following the rule established in Li et al. (2007).
  By default \code{winstp} is the time step.
}

\value{
  Returns an object of class \code{STBEUFit}.
  An object of class \code{STBEUFit} is a list containing
  at most the following components:
  \item{param}{The vector of  parameters' estimates;}
}


\references{
  Composite-likelihood:

  Varin, C., Reid, N. and Firth, D. (2011). An Overview of Composite
  Likelihood Methods. \emph{Statistica Sinica}, \bold{21}, 5--42.

  Varin, C. and Vidoni, P. (2005)
  A Note on Composite Likelihood Inference and Model Selection.
  \emph{Biometrika}, \bold{92}, 519--528.


  Weighted Composite-likelihood for binary RFs:

   Patrick, J. H. and Subhash, R. L. (1998)
  A Composite Likelihood Approach to Binary Spatial Data.
  \emph{Journal of the American Statistical Association, Theory & Methods}, \bold{93}, 1099--1111.


  Weighted Composite-likelihood for Gaussian RFs:

    Bevilacqua, M. Gaetan, C., Mateu, J. and Porcu, E. (2012)
  Estimating space and space-time covariance functions for large data
  sets: a weighted composite likelihood approach.
  \emph{Journal of the American Statistical Association, Theory &
    Methods}, \bold{107}, 268--280.

  Bevilacqua, M. Gaetan, C. (2013)
  On composite likelihood inference based on pairs for
spatial Gaussian RFs
  \emph{Techical Report, Department of Statistics, de Valparaiso University}.


  Sub-sampling estimation:

  Carlstein, E. (1986) The Use of Subseries Values for Estimating the
  Variance.
  \emph{The Annals of Statistics},
  \bold{14}, 1171--1179.

  Heagerty, P. J. and Lumley T. (2000) Window Subsampling of Estimating
  Functions with Application to Regression Models.
  \emph{Journal of the American Statistical Association, Theory &
    Methods}, \bold{95}, 197--211.

  Lee, Y. D. and Lahiri S. N. (2002) Variogram Fitting by Spatial
  Subsampling.
  \emph{Journal of the Royal Statistical Society. Series B},
  \bold{64}, 837--854.

  Li, B., Genton, M. G. and Sherman, M. (2007). A nonparametric assessment of properties of space-time
  covariance functions. \emph{Journal of the American Statistical
    Association}, \bold{102}, 736--744
}



\author{Moreno Bevilacqua, \email{moreno.bevilacqua@uv.cl},\url{https://sites.google.com/a/uv.cl/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/}
}
  
\examples{
#################################********   SPACE   ***####################################### 

#=======================Double Exponential=======================#

rm(list=ls())
########################R package#######################
require(GeoModels)
require(MCMCpack)
require(STBEU)
####location sites ########################################
lambda=8
xx=seq(-lambda,lambda);
coords=as.matrix(expand.grid(xx,xx))   ###regular
# plot(coords)
#set.seed(15)                                                                      ### not regular
#pp<-runifpoint(4*(lambda)^2, win=owin(c(-lambda,lambda),c(-lambda,lambda)))       ### not regular
#coords<-cbind(pp$x,pp$y)  
### not regular
####temporal instants ########################################
times=seq(1,5,1)
##################################################################


type_dist=1                ### type of distance     1=euclidean 2=chordal  3=geodesic 
maxdist=2            ## compact support in weights function for pairwise liklihood
maxtime=2

model=1  #   1=double exponential       2 =gneiting

if (model == 1) {
  # exponential model
  cov.model <-"exp_exp"
  cc=1
  #####
  mean=0
  nugget=0
  scale_s<-1.5/3
  scale_t<-1.5/3
  sill=1
  param=list(nugget=nugget,mean=mean,scale_t=scale_t,scale_s=scale_s,sill=sill)  
  fixed=list(nugget=0)
  fix=c(nugget=nugget)
}
##################################################################
set.seed(276)
#####################################################################################################      
# Simulation of the spatial Gaussian random field:
data <- GeoSim(coordx=coords,coordt=times,corrmodel=cov.model, param=param)$data
mm=mean(c(data))
vv=var(c(data))
###### Composite likelihood based on pairs estimation ###############################################     

start=list(mean=mm,scale_s=scale_s,scale_t=scale_t,sill=vv)
# Maximum composite-likelihood fitting of the random field:
fit <- GeoFit(data=data,coordx=coords,coordt=times,
                    corrmodel=cov.model,maxtime=maxtime,maxdist=maxdist,
                    likelihood="Marginal",type="Pairwise",
                    start=start,fixed=fixed,weighted=TRUE)
# Results:
print(fit$param)
unlist(start)
################################################
#   parameters for the subsampling ####
###############################################
coordx=coords[,1]
coordy=coords[,2]
LX=abs(range(coordx)[1]-range(coordx)[2])
LY=abs(range(coordy)[1]-range(coordy)[2])

lato_fin=3  #changing window size
lx=lato_fin          #lunghezza lato x quadrato subfinestra
ly=lato_fin          #lunghezza lato y quadrato subfinestra
winc=c(lx/sqrt(LX),ly/sqrt(LY))
winstp= 1  ###   1/lato_fin complete overlapping   1 "no" overlapping
###############################################

#####################################################################################################      
theta=start                #starting value  

weighted=0
### eucliden likelihood ################ 
type_subs=1    ### type of subsampling  1=in space    2= in time
tCPU = proc.time()
# names(fix) = c("nugget");names(theta) = c("mean","sill","scale_s","scale_t")
res=STBEUFit(theta,fixed = fix,coords,times,cc,data,type_dist,maxdist ,maxtime,winc,winstp,0,0,type_subs,weighted)
tCPU = proc.time()-tCPU;tCPU
res$par


### OpenCL eucliden likelihood ################ 
local <- c(1,1)
GPU <- 0
tGPU = proc.time()
res1 = STBEUFit(theta,fixed = fix,coords,times,cc,data,type_dist,maxdist ,maxtime,winc,winstp,0,0,type_subs,weighted, local=local, GPU = GPU)
tGPU = proc.time()-tGPU;tGPU
res1$par
tCPU;tGPU




res1 = STBEUFit(theta,fixed = fix,coords,times,cc,data,type_dist,maxdist ,maxtime,winc,winstp,0,0,type_subs,weighted, local=local, GPU = GPU,varest = TRUE)




#=======================Gneiting=======================#

rm(list=ls())


########################R package#######################
require(GeoModels)
require(MCMCpack)
require(STBEU)
####location sites ########################################
lambda=8
xx=seq(-lambda,lambda);
coords=as.matrix(expand.grid(xx,xx))   ###regular
nrow(coords)
# plot(coords)
#set.seed(15)                                                                      ### not regular
#pp<-runifpoint(4*(lambda)^2, win=owin(c(-lambda,lambda),c(-lambda,lambda)))       ### not regular
#coords<-cbind(pp$x,pp$y)  
### not regular
####temporal instants ########################################
times=seq(1,5,1)
##############################################################


type_dist=1                ### type of distance     1=euclidean 2=chordal  3=geodesic 

maxdist=2                  ## compact support in weights function for pairwise liklihood
maxtime=2

model=2  #   1=double exponential       2 =gneiting

if (model == 2) {
  # gneiting model
  cov.model <-"gneiting"
  cc=2
  #####
  mean=0
  nugget=0
  scale_s<-1.5/3
  scale_t<-1.5/20
  sill=4
  power_s=1;power_t=1;sep=0.5
  param=list(nugget=nugget,mean=mean,scale_t=scale_t,scale_s=scale_s,sill=sill,
             power_s=power_s,power_t=power_t,sep=sep)  
  fixed=list(nugget=0,power_s=power_s,power_t=power_t,sep=sep)  
  fix=c(nugget = nugget, power_s = power_s, power_t = power_t,sep = sep) 
}
##################################################################

set.seed(276)
##################################################################
# Simulation of the spatial Gaussian random field:
data <- GeoSim(coordx=coords,coordt=times,corrmodel=cov.model, param=param)$data
mm=mean(c(data))
vv=var(c(data))
###################################################################
###################################################################

###### Composite likelihood based on pairs estimation #############

start=list(mean=mm,scale_s=scale_s,scale_t=scale_t,sill=vv)
# Maximum composite-likelihood fitting of the random field:
fit <- GeoFit(data=data,coordx=coords,coordt=times,
                    corrmodel=cov.model,maxtime=maxtime,maxdist=maxdist,
                    likelihood="Marginal",type="Pairwise",
                    start=start,fixed=fixed,weighted=TRUE)
# Results:
print(fit$param)
unlist(start)
################################################
#   parameters for the subsampling ####
###############################################
coordx=coords[,1]
coordy=coords[,2]
LX=abs(range(coordx)[1]-range(coordx)[2])
LY=abs(range(coordy)[1]-range(coordy)[2])


lato_fin=3  #changing window size
lx=lato_fin          #lunghezza lato x quadrato subfinestra
ly=lato_fin          #lunghezza lato y quadrato subfinestra
winc=c(lx/sqrt(LX),ly/sqrt(LY))
winstp= 1  ###   1/lato_fin complete overlapping   1 "no" overlapping
###############################################
###############################################     
theta=start                #starting value  

weighted=0
### eucliden likelihood ################ 
type_subs=1    ### type of subsampling  1=in space    2= in time
tCPU = proc.time()
res=STBEUFit(theta,fixed = fix,coords,times,cc,data,type_dist,maxdist ,maxtime,winc,winstp,0,0,type_subs,weighted)
tCPU = proc.time()-tCPU;tCPU
res$par


res=STBEUFit(theta,fixed = fix,coords,times,cc,data,type_dist,maxdist ,maxtime,winc,winstp,0,0,type_subs,weighted,varest = TRUE)

### OpenCL eucliden likelihood ################ 
local <- c(1,1)
GPU <- 0
tGPU = proc.time()
res1=STBEUFit(theta,fixed = fix,coords,times,cc,data,type_dist,maxdist ,maxtime,winc,winstp,0,0,type_subs,weighted, local=local, GPU = GPU)
tGPU = proc.time()-tGPU;tGPU
res1$par
tCPU;tGPU



res1=STBEUFit(theta,fixed = fix,coords,times,cc,data,type_dist,maxdist ,maxtime,winc,winstp,0,0,type_subs,weighted, local=local, GPU = GPU,varest = TRUE)





#################################********   TIME   ***####################################### 
#=======================Double Exponential=======================#
rm(list=ls())
########################R package#######################
require(GeoModels)
require(MCMCpack)
require(STBEU)
####location sites ########################################
lambda=2
xx=seq(-lambda,lambda);
coords=as.matrix(expand.grid(xx,xx))   ###regular
nrow(coords)

#set.seed(15)                                                                      ### not regular
#pp<-runifpoint(4*(lambda)^2, win=owin(c(-lambda,lambda),c(-lambda,lambda)))       ### not regular
#coords<-cbind(pp$x,pp$y)  
### not regular
####temporal instants ########################################
times=seq(1,100,1)


type_dist=1                ### type of distance     1=euclidean 2=chordal  3=geodesic 
maxdist=2
maxtime=2
##################################################################
##################################################################
model=1  #   1=double exponential       2 =gneiting

if (model == 1) {
  # exponential model
  cov.model <-"exp_exp"
  cc=1
  #####
  mean=0
  nugget=0
  scale_s<-1.5/3
  scale_t<-1.5/3
  sill=4
  param=list(nugget=nugget,mean=mean,scale_t=scale_t,scale_s=scale_s,sill=sill)  
  fixed=list(nugget=0)
  fix=c(nugget = nugget)
}
##################################################################
##################################################################
# Simulation of the spatial Gaussian random field:
data <- GeoSim(coordx=coords,coordt=times,corrmodel=cov.model, param=param)$data
mm=mean(c(data))
vv=var(c(data))
###################################################################

###### Composite likelihood based on pairs estimation #############
start=list(mean=mm,scale_s=scale_s,scale_t=scale_t,sill=vv)
# Maximum composite-likelihood fitting of the random field:
fit <- GeoFit(data=data,coordx=coords,coordt=times,
                    corrmodel=cov.model,maxtime=maxtime,maxdist=maxdist,
                    likelihood="Marginal",type="Pairwise",
                    start=start,fixed=fixed,weighted=TRUE)


# Results:
print(fit$param)
unlist(start)
################################################
#   parameters for the temporal subsampling ####
###############################################

winc=5    ###  length of temporal window
winstp=1 ###   0.5 half overlapping  1 "no" overlapping
###############################################

#################################################################
theta=start                #starting value 
weighted=0
### eucliden likelihood ################ 
type_subs=2    ### type of subsampling  1=in space    2= in time
tCPU = proc.time()
res=STBEUFit(theta,fixed = fix,coords,times,cc,data,type_dist,maxdist ,maxtime,0,0,winc,winstp,type_subs,weighted)
tCPU = proc.time()-tCPU;tCPU
eu_par=c(res$par[1],res$par[3],res$par[4],res$par[2])
names(eu_par)=names(fit$param)
print(eu_par)

res=STBEUFit(theta,fixed = fix,coords,times,cc,data,type_dist,maxdist ,maxtime,0,0,winc,winstp,type_subs,weighted,varest = TRUE)

### OpenCL eucliden likelihood ################ 
local <- c(1,1)
GPU <- 0
tGPU = proc.time()
res1=STBEUFit(theta,fixed = fix,coords,times,cc,data,type_dist,maxdist ,maxtime,0,0,winc,winstp,type_subs,weighted, local=local, GPU = GPU)
tGPU = proc.time()-tGPU;tGPU
eu_par1=c(res1$par[1],res1$par[3],res1$par[4],res1$par[2])
names(eu_par1)=names(fit$param)
print(eu_par1)
tCPU;tGPU




res1=STBEUFit(theta,fixed = fix,coords,times,cc,data,type_dist,maxdist ,maxtime,0,0,winc,winstp,type_subs,weighted, local=local, GPU = GPU,varest = TRUE)


#=======================Gneiting=======================#

rm(list=ls())
########################R package#######################
require(GeoModels)
require(MCMCpack)
require(STBEU)
####location sites ########################################
lambda=2
xx=seq(-lambda,lambda);
coords=as.matrix(expand.grid(xx,xx))   ###regular
nrow(coords)

#set.seed(15)                                                                      ### not regular
#pp<-runifpoint(4*(lambda)^2, win=owin(c(-lambda,lambda),c(-lambda,lambda)))       ### not regular
#coords<-cbind(pp$x,pp$y)  
### not regular
####temporal instants ########################################
times=seq(1,100,1)


type_dist=1                ### type of distance     1=euclidean 2=chordal  3=geodesic 
maxdist=2
maxtime=2
##################################################################
model=2  #   1=double exponential       2 =gneiting

if (model == 2) {
  # gneiting model
  cov.model <-"gneiting"
  cc=2
  #####
  mean=0
  nugget=0
  scale_s<-1.5/3
  scale_t<-1.5/20
  sill=4
  power_s=1;power_t=1;sep=0.5
  param=list(nugget=nugget,mean=mean,scale_t=scale_t,scale_s=scale_s,sill=sill,
             power_s=power_s,power_t=power_t,sep=sep)  
  fixed=list(nugget=0,power_s=power_s,power_t=power_t,sep=sep)  
  fix=c(nugget = nugget, power_s = power_s, power_t = power_t,sep = sep)  
}
##################################################################
##################################################################
# Simulation of the spatial Gaussian random field:
data <- GeoSim(coordx=coords,coordt=times,corrmodel=cov.model, param=param)$data
mm=mean(c(data))
vv=var(c(data))
###################################################################
###################################################################

###### Composite likelihood based on pairs estimation #############
start=list(mean=mm,scale_s=scale_s,scale_t=scale_t,sill=vv)
# Maximum composite-likelihood fitting of the random field:
fit <- GeoFit(data=data,coordx=coords,coordt=times,
                    corrmodel=cov.model,maxtime=maxtime,maxdist=maxdist,
                    likelihood="Marginal",type="Pairwise",
                    start=start,fixed=fixed,weighted=TRUE)


# Results:
print(fit$param)
unlist(start)
################################################
#   parameters for the temporal subsampling ####
###############################################

winc=4    ###  length of temporal window
winstp=1 ###   0.5 half overlapping  1 "no" overlapping
###############################################
theta=start               #starting value 
weighted=0
### eucliden likelihood ################ 
type_subs=2    ### type of subsampling  1=in space    2= in time
tCPU = proc.time()
res=STBEUFit(theta,fixed = fix,coords,times,cc,data,type_dist,maxdist ,maxtime,0,0,winc,winstp,type_subs,weighted)
tCPU = proc.time()-tCPU;tCPU
eu_par=c(res$par[1],res$par[3],res$par[4],res$par[2])
names(eu_par)=names(fit$param)
print(eu_par)


res=STBEUFit(theta,fixed = fix,coords,times,cc,data,type_dist,maxdist ,maxtime,0,0,winc,winstp,type_subs,weighted,varest = TRUE)

### OpenCL eucliden likelihood ################ 
local <- c(1,1)
GPU <- 0
tGPU = proc.time()
res1=STBEUFit(theta,fixed = fix,coords,times,cc,data,type_dist,maxdist ,maxtime,0,0,winc,winstp,type_subs,weighted, local=local, GPU = GPU)
tGPU = proc.time()-tGPU;tGPU
eu_par1=c(res1$par[1],res1$par[3],res1$par[4],res1$par[2])
names(eu_par1)=names(fit$param)
print(eu_par1)
tCPU;tGPU


res1=STBEUFit(theta,fixed = fix,coords,times,cc,data,type_dist,maxdist ,maxtime,0,0,winc,winstp,type_subs,weighted, local=local, GPU = GPU,varest = TRUE)












#################################********   SPACETIME   ***####################################### 
#=======================Double Exponential=======================#
rm(list=ls())
########################R package#######################
require(GeoModels)
require(MCMCpack)
require(STBEU)
####location sites ########################################
lambda=8
xx=seq(-lambda,lambda);
coords=as.matrix(expand.grid(xx,xx))   ###regular
nrow(coords)
# plot(coords)
#set.seed(15)                                                                      ### not regular
#pp<-runifpoint(4*(lambda)^2, win=owin(c(-lambda,lambda),c(-lambda,lambda)))       ### not regular
#coords<-cbind(pp$x,pp$y)  
### not regular
####temporal instants ########################################
times=seq(1,5,1)

type_dist=1                ### type of distance     1=euclidean 2=chordal  3=geodesic 
maxdist=2
maxtime=2

##################################################################
model=1  #   1=double exponential       2 =gneiting

if (model == 1) {
  # exponential model
  cov.model <-"exp_exp"
  cc=1
  #####
  mean=0
  nugget=0
  scale_s<-1.5/3
  scale_t<-1.5/3
  sill=4
  param=list(nugget=nugget,mean=mean,scale_t=scale_t,scale_s=scale_s,sill=sill)  
  fixed=list(nugget=0)
  fix=c(nugget = nugget)
}

##################################################################
set.seed(276)
##################################################################
# Simulation of the  Gaussian random field:
data <- GeoSim(coordx=coords,coordt=times,corrmodel=cov.model, param=param)$data
##################################################################
mm=mean(c(data))
vv=var(c(data))
###### Composite likelihood based on pairs estimation ############
start=list(mean=mm,scale_s=scale_s,scale_t=scale_t,sill=vv)

# Maximum composite-likelihood fitting of the random field:
fit <- GeoFit(data=data,coordx=coords,coordt=times,
                    corrmodel=cov.model,maxtime=maxtime,maxdist=maxdist,
                    likelihood="Marginal",type="Pairwise",
                    start=start,fixed=fixed,weighted=TRUE)

# Results:
print(fit$param)
unlist(start)
################################################
#   parameters for the subsampling ####
###############################################
coordx=coords[,1]
coordy=coords[,2]
LX=abs(range(coordx)[1]-range(coordx)[2])
LY=abs(range(coordy)[1]-range(coordy)[2])
lato_fin=3  #changing window size
lx=lato_fin          #lunghezza lato x quadrato subfinestra
ly=lato_fin          #lunghezza lato y quadrato subfinestra
winc=c(lx/sqrt(LX),ly/sqrt(LY))
winstp= 1###   1/lato_fin complete overlapping in space  1 "no" overlapping in space
###############################################
winc_t=4   ###  length of temporal window
winstp_t=1 ###   0.5 half overlapping  1 "no" overlapping
###############################################

theta=start                #starting value  
weighted=0                                          #weigthed version  1=yes 0=no
### eucliden likelihood ################ 
type_subs=3   ### type of subsampling  1=in space    2= in time  3 =spacetime
tCPU = proc.time()
res=STBEUFit(theta,fixed = fix,coords,times,cc,data,type_dist,maxdist ,maxtime,winc,winstp,winc_t,winstp_t,type_subs,weighted)
tCPU = proc.time()-tCPU;tCPU
eu_par=c(res$par[1],res$par[3],res$par[4],res$par[2])
names(eu_par)=names(fit$param)
print(eu_par)

res=STBEUFit(theta,fixed = fix,coords,times,cc,data,type_dist,maxdist ,maxtime,winc,winstp,winc_t,winstp_t,type_subs,weighted,varest = TRUE)

### OpenCL eucliden likelihood ################ 
local <- c(1,1)
GPU <- 0
tGPU = proc.time()
res1=STBEUFit(theta,fixed = fix,coords,times,cc,data,type_dist,maxdist ,maxtime,winc,winstp,winc_t,winstp_t,type_subs,weighted, local = local, GPU = GPU)
tGPU = proc.time()-tGPU;tGPU
eu_par1=c(res1$par[1],res1$par[3],res1$par[4],res1$par[2])
names(eu_par1)=names(fit$param)
print(eu_par1)
tCPU;tGPU


res1=STBEUFit(theta,fixed = fix,coords,times,cc,data,type_dist,maxdist ,maxtime,winc,winstp,winc_t,winstp_t,type_subs,weighted, local = local, GPU = GPU,varest = TRUE)




#=======================Gneting=======================#

rm(list=ls())
########################R package#######################
require(GeoModels)
require(MCMCpack)
require(STBEU)
####location sites ########################################
lambda=8
xx=seq(-lambda,lambda);
coords=as.matrix(expand.grid(xx,xx))   ###regular
nrow(coords)
# plot(coords)
#set.seed(15)                                                                      ### not regular
#pp<-runifpoint(4*(lambda)^2, win=owin(c(-lambda,lambda),c(-lambda,lambda)))       ### not regular
#coords<-cbind(pp$x,pp$y)  
### not regular
####temporal instants ########################################
times=seq(1,5,1)

type_dist=1                ### type of distance     1=euclidean 2=chordal  3=geodesic 
maxdist=2
maxtime=2
##################################################################
model=2  #   1=double exponential       2 =gneiting


if (model == 2) {
  # gneiting model
  cov.model <-"gneiting"
  cc=2
  #####
  mean=0
  nugget=0
  scale_s<-1.5/3
  scale_t<-1.5/20
  sill=4
  power_s=1;power_t=1;sep=0.5
  param=list(nugget=nugget,mean=mean,scale_t=scale_t,scale_s=scale_s,sill=sill,
             power_s=power_s,power_t=power_t,sep=sep)  
  fixed=list(nugget=0,power_s=power_s,power_t=power_t,sep=sep)   
  fix=c(nugget = nugget, power_s = power_s, power_t = power_t,sep=sep) 
}
##################################################################
set.seed(276)
##################################################################
# Simulation of the  Gaussian random field:
data <- GeoSim(coordx=coords,coordt=times,corrmodel=cov.model, param=param)$data
##################################################################
##################################################################
mm=mean(c(data))
vv=var(c(data))
###### Composite likelihood based on pairs estimation ############
start=list(mean=mm,scale_s=scale_s,scale_t=scale_t,sill=vv)

# Maximum composite-likelihood fitting of the random field:
fit <- GeoFit(data=data,coordx=coords,coordt=times,
                    corrmodel=cov.model,maxtime=maxtime,maxdist=maxdist,
                    likelihood="Marginal",type="Pairwise",
                    start=start,fixed=fixed,weighted=TRUE)

# Results:
print(fit$param)
unlist(start)
################################################
#   parameters for the subsampling ####
###############################################
coordx=coords[,1]
coordy=coords[,2]
LX=abs(range(coordx)[1]-range(coordx)[2])
LY=abs(range(coordy)[1]-range(coordy)[2])
lato_fin=3  #changing window size
lx=lato_fin          #lunghezza lato x quadrato subfinestra
ly=lato_fin          #lunghezza lato y quadrato subfinestra
winc=c(lx/sqrt(LX),ly/sqrt(LY))
winstp= 1###   1/lato_fin complete overlapping in space  1 "no" overlapping in space
###############################################
winc_t=4   ###  length of temporal window
winstp_t=1 ###   0.5 half overlapping  1 "no" overlapping
###############################################

theta=start                #starting value  
weighted=0                                          #weigthed version  1=yes 0=no
### eucliden likelihood ################ 
type_subs=3   ### type of subsampling  1=in space    2= in time  3 =spacetime
tCPU = proc.time()
res=STBEUFit(theta,fixed = fix,coords,times,cc,data,type_dist,maxdist ,maxtime,winc,winstp,winc_t,winstp_t,type_subs,weighted)
tCPU = proc.time()-tCPU;tCPU
eu_par=c(res$par[1],res$par[3],res$par[4],res$par[2])
names(eu_par)=names(fit$param)
print(eu_par)

res=STBEUFit(theta,fixed = fix,coords,times,cc,data,type_dist,maxdist ,maxtime,winc,winstp,winc_t,winstp_t,type_subs,weighted,varest = TRUE)


### OpenCL eucliden likelihood ################ 
local <- c(1,1)
GPU <- 0
tGPU = proc.time()
res1=STBEUFit(theta,fixed = fix,coords,times,cc,data,type_dist,maxdist ,maxtime,winc,winstp,winc_t,winstp_t,type_subs,weighted, local = local, GPU = GPU)
tGPU = proc.time()-tGPU;tGPU
eu_par1=c(res1$par[1],res1$par[3],res1$par[4],res1$par[2])
names(eu_par1)=names(fit$param)
print(eu_par1)
tCPU;tGPU

res1=STBEUFit(theta,fixed = fix,coords,times,cc,data,type_dist,maxdist ,maxtime,winc,winstp,winc_t,winstp_t,type_subs,weighted, local = local, GPU = GPU,varest = TRUE)




#=======================Wen_time=======================#
###########################################################
require(GeoModels)
library(spam)
require(fields)
require(STBEU)
require(MCMCpack)

corrmodel="Wen_time"; 


####################################################
# Define the spatial-coordinates of the points:
set.seed(8)
dp <- 3
ns <- 7
x <- runif(ns)
y <- runif(ns)
coords=cbind(x,y)


n=nrow(coords)
# Define the temporal-coordinates:
times <- seq(0, 14.75, .3)
# times <- seq(0, 14.75, 4)
tt=length(times)
####################################################


####################################################
smooth_t=2 # k or kappa
scale_t=time_comp_supp=0.75     # compact supportt scale_t
scale_s=0.1
power2_t=3.5+smooth_t #nu
power_s=2
power2_s=2.5+2*smooth_t# tau

sep=1  ## 0 0.5  1
sill=1
mean=0
nugget=0
####################################################

####################################################
param <- list(nugget=nugget,mean=mean,sill = sill,
              power2_s=power2_s,power_s=power_s,power2_t=power2_t,scale_s=scale_s,
              scale_t=scale_t,sep=sep,smooth_t=smooth_t)
fixed=list(nugget=nugget,power2_s=power2_s,power2_t=power2_t,sep=sep,mean = mean,power_s=power_s)  
start=list(scale_s=scale_s,scale_t=time_comp_supp,sill=sill,smooth_t=smooth_t)
####################################################
# Simulation of a spatial Gaussian RF:

set.seed(8519)

data <- GeoSim(coordx=coords, coordt=times, corrmodel=corrmodel, 
               param=param)$data


# Maximum composite-likelihood fitting of the random field:
tCPU0 <- proc.time()
fit <- GeoFit(data=data,coordx=coords,coordt=times,corrmodel=corrmodel,
              likelihood="Full",type="Standard",sparse=TRUE,
              start=start,fixed=fixed)
tCPU0 <- proc.time()-tCPU0
print(tCPU0)


model=3  #   1=double exponential       2 =gneiting 3 = WenTime

# gneiting model
cov.model <-"Wen_time"
cc=3 
#####

################################################
#   parameters for the temporal subsampling ####
###############################################

winc=2    ###  length of temporal window
winstp=1 ###   0.5 half overlapping  1 "no" overlapping
###############################################

#################################################################
theta=start                #starting value 
weighted=0
type_dist=1                ### type of distance     1=euclidean 2=chordal  3=geodesic 
maxdist=.4
maxtime=2
### eucliden likelihood ################ 
type_subs=2    ### type of subsampling  1=in space    2= in time
tCPU = proc.time()
res=STBEUFit(theta,fixed = unlist(fixed),coords,times,cc,data,type_dist,
             maxdist ,maxtime,0,0,winc,winstp,type_subs,weighted)
tCPU = proc.time()-tCPU
print(tCPU)


tGPU = proc.time()
res1=STBEUFit(theta,unlist(fixed),coords,times,cc,data,type_dist,
              maxdist ,maxtime,0,0,winc,
              winstp,type_subs,weighted,
              GPU = 0, local = c(1,1))
tGPU = proc.time()-tGPU
print(tGPU)

}
\keyword{Composite}
