#include "header.h"


/****************************************************************************************************************/
void scalar_space(int *npts, int *ntime,double *coordt, double *maxtime,double *maxdist,int *cormod, double *parcor,int *flagcor,int *flagnuis,int *npar,double *nuis,double *sdata,int *weigthed, double *mom_cond, int *dist, double *scoordx,double *scoordy,double *gradcor,double  *grad, double *ww )
{
    int l=0,m=0,t =0, v =0 ,kk=0;
    double lagt=0.0,lags=0.0, rho=0;
    
    for(l=0;l<npts[0];l++)
    {
        for(t=0;t<*ntime;t++)
        {
            for(m=l;m<npts[0];m++)
            {
                if(l==m)
                {
                    for(v=t+1;v<*ntime;v++)
                    {
                        lagt=fabs(coordt[t]-coordt[v]);
                        if(lagt<=maxtime[0])
                        {
                            
                            //Computing correlation
                            rho=CorFct(cormod,0,lagt,parcor);
                            //Computing the gradient of the corr parameters
                           
                            GradCorrFct(rho,cormod,flagcor,gradcor,0,lagt,parcor);
                            //Compute the gradient of the composite likelihood:
                            Grad_Pair_Gauss(rho,flagnuis,gradcor,grad,npar,nuis,sdata[(t+*ntime*l)],sdata[(v+*ntime*l)]);
                            /*if(weigthed[0])
                             {
                             weights=CorFunBohman(lagt,maxtime[0]);
                             ww[0]=1;ww[1]=1;ww[2]=weights;ww[3]=weights;
                             }
                             else
                             {
                             ww[0]=1;ww[1]=1;ww[2]=1; ww[3]=1;
                             }*/
                            //printf("npar:%d  g1 %f g2 %f g3 %f   gc %f gc %f\n",npar[0],grad[0],grad[1],grad[2],gradcor[0],gradcor[1]);
                            ww[0]=1;ww[1]=1;ww[2]=1; ww[3]=1;
                            for(kk=0;kk<npar[0];kk++)
                            {mom_cond[kk]=mom_cond[kk]+ww[kk]*grad[kk];
                            }
                        }
                    }
                }
                else
                {
                    if(*dist==1) lags=hypot(scoordx[l]-scoordx[m],scoordy[l]-scoordy[m]); // euclidean lag
                    for(v=0;v<*ntime;v++)
                    {
                        lagt=fabs(coordt[t]-coordt[v]);
                        if(lagt<=maxtime[0] && lags<=maxdist[0])
                        {
                            //Computing correlation
                            rho=CorFct(cormod,lags,lagt,parcor);
                            //Computing the gradient of the corr parameters
                            GradCorrFct(rho,cormod,flagcor,gradcor,lags,lagt,parcor);
                            //Compute the gradient of the composite likelihood:
                            Grad_Pair_Gauss(rho,flagnuis,gradcor,grad,npar,nuis,sdata[(t+*ntime*l)],sdata[(v+*ntime*m)]);
                            /*if(weigthed[0])
                             {
                             weights=CorFunBohman(lags,maxdist[0])*CorFunBohman(lagt,maxtime[0]);
                             ww[0]=1;ww[1]=1;ww[2]=weights;ww[3]=weights;
                             }
                             else
                             {
                             ww[0]=1;ww[1]=1;ww[2]=1; ww[3]=1;
                             }*/
                            ww[0]=1;ww[1]=1;ww[2]=1; ww[3]=1;
                            for(kk=0;kk<npar[0];kk++)
                            {mom_cond[kk]=mom_cond[kk]+ww[kk]*grad[kk];
                            }
                        }
                    }
                }
            }
        }
    }
    //printf("final result: %.4f\t%.4f\t%.4f\t%.4f\n",mom_cond[0],mom_cond[1],mom_cond[2],mom_cond[3]);
}

/*******************************************************************************************************************************************/

void SubSamp_space(double *coordx, double *coordy,double *coordt, int *ncoord,int *ntime,int *cormod,double *data,int *dist, double *maxdist,double *maxtime,int *npar,double *parcor,int *nparc,double *nuis,int *nparnuis, int *flagcor, int *flagnuis,double *vari,double *winc, double *winstp,double *a, double *b,double *block_mean,int *weigthed, int *local_wi, int *dev)
{
    double   *rangex, *rangey;
    double *vv,*sdata,*xgrid,*ygrid,*scoordx,*scoordy;
    double *gradcor, *grad, *ww;
    double **tot,*mom_cond,**vector_mean;    //matrix moments conditions
    double delta=0.0, deltax=0.0, deltay=0.0, dimwinx=0.0, dimwiny=0.0;
    double winstx=winstp[0], winsty=winstp[0];
    int *npts, numintx=0, numinty=0;
    int n_win=0,kk=0,h=0,i=0,nsub,j=0,p=0,q=0,r,nmat;
    
    
    
    npts=(int *) R_alloc(1, sizeof(int));          //number of location sites in the subwindow
    rangex=(double *) R_alloc(2, sizeof(double));
    rangey=(double *) R_alloc(2, sizeof(double));
    
    scoordx=(double *) Calloc(*ncoord,double);
    scoordy=(double *) Calloc(*ncoord,double);
    sdata=(double *) Calloc(*ncoord * *ntime,double);
    
    Range(coordx,rangex,ncoord[0]);// range of the x-coordinate (min and max)
    Range(coordy,rangey,ncoord[0]);// range of the y-coordinate  (min and max)
    // set the sub-sampling window based on prototype unit window (R_0)
    // and scaling factor (lambda_n)
    deltax=rangex[1]-rangex[0];//  x window length
    deltay=rangey[1]-rangey[0];//  y window length
    
    //Rprintf("A. MIERDAAA\n");
    
    if(!winc[0])
    {
        delta=fmin(deltax,deltay);  // I set this rule if no constant
        winc[0]=sqrt(delta)/2;
    }
    if(!winstp[0]) winstp[0]=0.5;   //proportion of the overlapping  0< winstp <=1, If winstp=0 then  winstp=0.5
    //Rprintf("%f\n",winstp[0]);
    dimwinx=winc[0] * sqrt(deltax);// sub-window x length depends on a constant: deafault??
    dimwiny=winc[0] * sqrt(deltay);// sub-window y length depends on a constant: deafault??
    
    winstx=*winstp * dimwinx;     // x step is a  proportion of sub-window x length (deafult is 0.5)
    winsty=*winstp * dimwiny;     // y step is a  proportion of sub-window y length (deafult is 0.5)
    numintx=floor((deltax-dimwinx)/(winstx)+1);   //number of overlapping sub-windows is  numintx+1 * numinty+1
    numinty=floor((deltay-dimwiny)/(winsty)+1);
    //printf("--------- numintx %d numinty %d",numintx,numinty);
    //printf("--------- numintx %d numinty %d   deltax %f deltay %f dimwinx: %f dimwiny: %f n_win %d winstx %f winstx %f cc1 %f cc2 %f \n",numintx,numinty,deltax,deltay,dimwinx,dimwiny,n_win,winstx,winsty,((deltax-dimwinx)/(winstx)+1),((deltay-dimwiny)/(winsty)+1));
    
    xgrid=(double *) R_alloc(numintx, sizeof(double));
    ygrid=(double *) R_alloc(numinty, sizeof(double));
    
    SeqStep(rangex,numintx,winstx,xgrid);
    SeqStep(rangey,numinty,winsty,ygrid);
    
    // matrix of means (one for each window)
    
    vector_mean= (double **) Calloc(npar[0],double *);
    for(i=0;i<npar[0];i++)
    {
        vector_mean[i]=(double *) Calloc((numintx+1)*(numinty+1),double);
    }
    
    ww=(double *) Calloc(*npar,double);
    gradcor=(double *) Calloc(*nparc,double);
    grad=(double *) Calloc(*npar,double);
    
    nsub=0;
    n_win=(numintx)*(numinty);   //// number of windows
    
    
    for(i=0;i<=numintx;i++)
    {
        for(j=0;j<=numinty;j++)
        {  // cycle for each block window
            *npts=0;   // number of loc sites in the window
            SetSampling_s(ncoord[0],ntime[0],coordx,coordy,data,npts,scoordx,scoordy,
                          sdata,xgrid[i]+dimwinx,xgrid[i],ygrid[j]+dimwiny,ygrid[j]);//  create data and coordinates of the sub-windows
            //  scoordx    scoordy are the spatial coords of the subwindow
            // sdata the associated data
            // npts number of loc sites in the window
            if( (((xgrid[i]+dimwinx)<rangex[1])||is_equal(xgrid[i]+dimwinx,rangex[1]))&&     //excludin "half" windows
               (((ygrid[j]+dimwiny))<rangey[1]||is_equal(ygrid[j]+dimwiny,rangey[1]))&&       //and windows with few loc sites
               (npts[0]>5))
            {
                mom_cond=(double *) Calloc(*npar,double);
                /******************************************/
                /*computing gradient in the window*/
                /******************************************/
                
                scalar_space(npts,ntime,coordt,maxtime,maxdist,cormod,parcor,flagcor,flagnuis,npar,nuis,sdata,weigthed,mom_cond,dist, scoordx,scoordy,gradcor,grad,ww);
                
                
                /******************************************/
                /******************************************/
                //for(l=0;l<nsens;l++) sensmat[l]=sensmat[l]-sens[l];
                for(kk=0;kk<npar[0];kk++)
                {
                    vector_mean[kk][nsub]=mom_cond[kk];
                }  //vector means for each winndow
                nsub=nsub+1;  // counting number of windows
                Free(mom_cond);
            }
        }
    }
    
    Free(scoordx);
    Free(scoordy);
    Free(sdata);
    Free(grad);
    Free(gradcor);
    Free(ww);
    
    tot= (double **) Calloc(npar[0],double *);
    for(i=0;i<npar[0];i++)
    {
        tot[i]=(double *) Calloc(n_win,double);
    }
    //mean over blocks
    for(q=0;q<n_win ;q++)
    {
        for(kk=0;kk<*npar;kk++)
        {
            block_mean[kk]= block_mean[kk]+ vector_mean[kk][q]/n_win;
        }
    }
    for(q=0;q<n_win ;q++)
    {
        for(kk=0;kk<*npar;kk++)
        {
            tot[kk][q]=  vector_mean[kk][q]+block_mean[kk];
        }
    }
    ////computing variance covarianmce matrix
    nmat=0.5*npar[0]*(npar[0]-1)+npar[0];
    vv=(double *) Calloc(n_win*nmat,double);
    h=0;
    for(r=0;r<n_win ;r++)
    {
        for(p=0;p< *npar;p++)
        {
            for(q=p;q< *npar;q++)
            {
                vv[h]=vv[h]+(tot[p][r])*(tot[q][r])/(n_win);
                h++;
            }
        }
    }
    ////building upper triangular of the variance covarianmce matrix (summing over the matrices )
    for(i=0;i<nmat;i++)
    {
        for(q=0;q<n_win ;q++)
        {
            vari[i]=vari[i]+vv[i+q*nmat];
        }
    }
    Free(vv);
    for(i=0;i<npar[0];i++)
    {
        Free(vector_mean[i]);Free(tot[i]);
    }
    Free(vector_mean);
    Free(tot);
    
    return;
}
/*******************************************************************************************************************************************/

void scalar_time(int *ncoord,int *nstime,double *sublagt,int *cormod,double *parcor,int *flagcor, double *gradcor,int *flagnuis, double *grad,int *npar,double *nuis, double *sdata,int *weigthed,double *maxtime, double *ww, double *mom_cond,int *dist, double *coordx, double *coordy,double *maxdist)
{
    //Rprintf("C: R_power_s: %f R_power: %f R_power_t: %f scale_s: %f scale_t: %f sep: %f smooth: %f \n",parcor[1],parcor[0],parcor[2],parcor[3],parcor[4],parcor[5],parcor[6]);
    // 0,0,0,1,1,0,1
     //printf("\n\nA: C FLAGCOR: %d\t%d\t%d\t%d\t%d\t%d\t%d\n\n",flagcor[0],flagcor[1],flagcor[2],flagcor[3],flagcor[4],flagcor[5],flagcor[6]);
    //printf("CORMOD: %d\n",cormod[0]);
    
    int l= 0, t=0, m=0, v=0,kk=0;
    double lags=0.0,lagt=0.0,rho=0.0;
    
    
    for(l=0;l<ncoord[0];l++)
    {
        for(t=0;t<nstime[0];t++)
        {
            for(m=l;m<ncoord[0];m++)
            {
                if(l==m)
                {
                    for(v=t+1;v<nstime[0];v++)
                    {
                        lagt=fabs(sublagt[t]-sublagt[v]);
                        if(lagt<=maxtime[0])
                        {
                            //Computing correlation
                            rho=CorFct(cormod,0,lagt,parcor);
                            //Computing the gradient of the corr parameters
                            
                            GradCorrFct(rho,cormod,flagcor,gradcor,0,lagt,parcor);
                            
                            //printf("A: C gradcor: %f\t%f\t%f    rho: %f\n\n",gradcor[0],gradcor[1],gradcor[2],rho);
                           
                            //Compute the gradient of the composite likelihood:
                            Grad_Pair_Gauss(rho,flagnuis,gradcor,grad,npar,nuis,sdata[(t+nstime[0]*l)],sdata[(v+nstime[0]*l)]);
                            /*if(*weigthed)
                             {
                             weights=CorFunBohman(lagt,maxtime[0]);
                             ww[0]=1;ww[1]=1;ww[2]=weights;ww[3]=weights;
                             }
                             else
                             {
                             ww[0]=1;ww[1]=1;ww[2]=1; ww[3]=1;
                             }*/
                            ww[0]=1;ww[1]=1;ww[2]=1; ww[3]=1;
                            for(kk=0;kk<npar[0];kk++)
                            {
                                //Rprintf("A: mom_cond[kk]: %f  grad[kk]: %f\n",mom_cond[kk],grad[kk]);
                                mom_cond[kk]=mom_cond[kk]+ww[kk]*grad[kk];
                            }
                        }
                    }
                }
                else
                {
                    if(*dist==1) lags=hypot(coordx[l]-coordx[m],coordy[l]-coordy[m]);    // lag distance in the subwindow
                    // if(*dist==2) lags=Dist_geodesic(coordx[l],coordy[l],coordx[m],coordy[m]);
                    for(v=0;v<nstime[0];v++)
                    {
                        lagt=fabs(sublagt[t]-sublagt[v]);
                        if(lagt<=maxtime[0] && lags<=maxdist[0])
                        {
                            //Computing correlation
                            rho=CorFct(cormod,lags,lagt,parcor);
                            //Computing the gradient of the corr parameters
                            //printf("B: C gradcor: %f\t%f\t%f\n\n",gradcor[0],gradcor[1],gradcor[2]);
                            GradCorrFct(rho,cormod,flagcor,gradcor,lags,lagt,parcor);
                            //Compute the gradient of the composite likelihood:
                            Grad_Pair_Gauss(rho,flagnuis,gradcor,grad,npar,nuis,sdata[(t+nstime[0]*l)],sdata[(v+nstime[0]*m)]);
                            /*if(*weigthed)
                             {
                             weights=CorFunBohman(lags,maxdist[0])*CorFunBohman(lagt,maxtime[0]);
                             ww[0]=1;ww[1]=1;ww[2]=weights;ww[3]=weights;
                             }
                             else
                             {
                             ww[0]=1;ww[1]=1;ww[2]=1; ww[3]=1;
                             }*/
                            ww[0]=1;ww[1]=1;ww[2]=1; ww[3]=1;
                            for(kk=0;kk<npar[0];kk++)
                            {
                                //Rprintf("B: mom_cond[kk]: %f  grad[kk]: %f\n",mom_cond[kk],grad[kk]);
                                mom_cond[kk]=mom_cond[kk]+ww[kk]*grad[kk];
                            }
                        }
                    }
                }
            }
        }
    }
}

/*******************************************************************************************************************************************/


void SubSamp_time(double *coordx, double *coordy,double *coordt, int *ncoord,int *ntime,int *cormod,double *data,int *dist, double *maxdist,double *maxtime,int *npar,double *parcor,int *nparc,double *nuis,int *nparnuis, int *flagcor, int *flagnuis,double *vari,double *a, double *b,double *winc, double *winstp,double *block_mean,int *weigthed,int *local_wi, int *dev)
{
    double beta, *gradcor;
    double *vv,*ww,*sdata, *grad, *sublagt;
    double **tot,*mom_cond,**vector_mean;    //matrix moments conditions
    double step=coordt[0]-coordt[1];  // subsampling works in the regular temporal sampling setting
    int nsub=0, nstime=0;
    int kk=0,h=0,i=0,p=0,q=0,r,nmat;
    //printf("step: %f\n",step);
    //default sub window temporal length
    if(!(*winc))
    {
        beta=CorFct(cormod,0,1,parcor);
        *winc=R_pow(2*beta/(1-R_pow(beta,2)),(double) 2/3)*R_pow(*ntime*1.5,(double)1/3);
        if(*winc<4*step) *winc=2*step;// if the length is too small
        if(*winc>=*ntime) *winc=*ntime-step;
    } // if the length is too big
    //set the spatial-temporal windows:
    int wint = (int) *winc;
    //printf("wint: %d\n",wint);
    if(*winstp==0) *winstp=wint; //defualt for the forward step:the minimum distance
    //else *winstp=*winstp*wint;   //otherwise a proportion of the temporal window
    sublagt=(double *) R_alloc(wint, sizeof(double));
    sublagt[0]=step;
    sdata=(double *)  Calloc(*ncoord*wint,double);
    //set the temporal distances for the sub-sample:
    for(i=0;i<wint;i++)
    {
        sublagt[i+1]=sublagt[i]+step;nstime++;
    }
    nsub=floor((( (*ntime)-wint)/(wint*winstp[0])+1));
    // vector of means for each window
    vector_mean= (double **) Calloc(npar[0],double *);
    for(i=0;i<npar[0];i++) vector_mean[i]=(double *) Calloc(nsub,double);
    //start the sub-sampling procedure:
    
    ww=(double *) Calloc(*npar,double);
    gradcor=(double *) Calloc(*nparc,double);
    grad=(double *) Calloc(*npar,double);
    
    for(i=0;i<nsub;i++)
    {//loop for the number of sub-sampling:
        mom_cond=(double *) Calloc(*npar,double);
        // set the sub-sample of the data:
        SetSampling_t(data,sdata,*ncoord,*ntime,wint,i);
        
        
        /******************************************/
        /*computing gradient in the window*/
        /******************************************/
      // printf("C gradcor: %f\t%f\t%f\n\n",gradcor[0],gradcor[1],gradcor[2]);
        scalar_time(ncoord,&nstime,sublagt,cormod,parcor,flagcor,gradcor,flagnuis,grad,npar,nuis,sdata,weigthed,maxtime,ww,mom_cond,dist,coordx,coordy,maxdist);
        
        /*if(mom_cond[0]<0 || mom_cond[1]<0 ||mom_cond[2]<0)
        {
            Rprintf("%f\t%f\t%f\n",mom_cond[0],mom_cond[1],mom_cond[2]);
        }*/
        
        //printf("C: %f\t%f\t%f\t%f\n",mom_cond[0],mom_cond[1],mom_cond[2],mom_cond[3]);
        /******************************************/
        /******************************************/
        for(kk=0;kk<npar[0];kk++)
        {
            //Rprintf("mom_cond[kk]: %f\n",mom_cond[kk]);
            vector_mean[kk][i]=mom_cond[kk];
        }  //vector means for eac winndow
        Free(mom_cond);
        
    }
    
    Free(sdata);Free(grad);Free(gradcor);Free(ww);
    
    tot= (double **) Calloc(npar[0],double *);
    for(i=0;i<npar[0];i++)tot[i]=(double *) Calloc(nsub,double);
    //mean over blocks
    for(q=0;q<nsub ;q++)
    {
        for(kk=0;kk<*npar;kk++)
        {
            block_mean[kk]= block_mean[kk]+ vector_mean[kk][q]/nsub;
            //printf("block_mean[kk]:%f\n",block_mean[kk]);//problema
            /*if(ISNAN(block_mean[kk]))
            {
                Rprintf("block_mean[kk]: %f vector_mean[kk][q]: %f kk: %d\n",block_mean[kk],vector_mean[kk][q],kk);
            }*/
        }
    }
    //if(block_mean[0]<0 ){return;}
    for(q=0;q<nsub ;q++)
    {
        for(kk=0;kk<*npar;kk++)
        {
            tot[kk][q]=  vector_mean[kk][q]+block_mean[kk];
            //printf("tot[kk][q]: %f\n", tot[kk][q]);//Problema
        }
    }
    ////computing variance covarianmce matrix
    nmat=0.5*npar[0]*(npar[0]-1)+npar[0];
    vv=(double *) Calloc(nsub*nmat,double);
    h=0;
    for(r=0;r<nsub ;r++)
    {
        for(p=0;p< *npar;p++)
        {
            for(q=p;q< *npar;q++)
            {
                vv[h]=vv[h]+(tot[p][r])*(tot[q][r])/(nsub);
                //printf("vv[h]: %f\n", vv[h]); // PROBLEMA
                h++;
            }
        }
    }
    ////building upper triangular of the variance covarianmce matrix (summing over the matrices )
    for(i=0;i<nmat;i++)
    {
        for(q=0;q<nsub ;q++)
        {
            vari[i]=vari[i]+vv[i+q*nmat];
        }
    }
    Free(vv);
    for(i=0;i<npar[0];i++)
    {
        Free(vector_mean[i]);Free(tot[i]);
    }
    Free(vector_mean);
    Free(tot);
    
    
    return;
}

/*void SubSamp_time1(double *coordx, double *coordy,double *coordt, int *ncoord,int *ntime,int *cormod,double *data,int *dist, double *maxdist,double *maxtime,int *npar,double *parcor,int *nparc,double *nuis,int *nparnuis, int *flagcor, int *flagnuis,double *vari,double *a, double *b,double *winc, double *winstp,double *block_mean,int *weigthed,int *local_wi, int *dev,double *grad)
{
    printf("Ok\n");
    double beta, *gradcor;
    double *vv,*ww,*sdata,  *sublagt;
    double **tot,*mom_cond,**vector_mean;    //matrix moments conditions
    double step=coordt[0]-coordt[1];  // subsampling works in the regular temporal sampling setting
    int nsub=0, nstime=0;
    int kk=0,h=0,i=0,p=0,q=0,r,nmat;
    //printf("step: %f\n",step);
    //default sub window temporal length
    if(!(*winc))
    {
        beta=CorFct(cormod,0,1,parcor);
        *winc=R_pow(2*beta/(1-R_pow(beta,2)),(double) 2/3)*R_pow(*ntime*1.5,(double)1/3);
        if(*winc<4*step) *winc=2*step;// if the length is too small
        if(*winc>=*ntime) *winc=*ntime-step;
    } // if the length is too big
    //set the spatial-temporal windows:
    int wint = (int) *winc;
    //printf("wint: %d\n",wint);
    if(*winstp==0) *winstp=wint; //defualt for the forward step:the minimum distance
    //else *winstp=*winstp*wint;   //otherwise a proportion of the temporal window
    sublagt=(double *) R_alloc(wint, sizeof(double));
    sublagt[0]=step;
    sdata=(double *)  Calloc(*ncoord*wint,double);
    //set the temporal distances for the sub-sample:
    for(i=0;i<wint;i++)
    {
        sublagt[i+1]=sublagt[i]+step;nstime++;
    }
    nsub=floor((( (*ntime)-wint)/(wint*winstp[0])+1));
    // vector of means for each window
    vector_mean= (double **) Calloc(npar[0],double *);
    for(i=0;i<npar[0];i++) vector_mean[i]=(double *) Calloc(nsub,double);
    //start the sub-sampling procedure:
    
    ww=(double *) Calloc(*npar,double);
    gradcor=(double *) Calloc(*nparc,double);
    //grad=(double *) Calloc(*npar,double);
    
    for(i=0;i<nsub;i++)
    {//loop for the number of sub-sampling:
        mom_cond=(double *) Calloc(*npar,double);
        // set the sub-sample of the data:
        SetSampling_t(data,sdata,*ncoord,*ntime,wint,i);

      // printf("C gradcor: %f\t%f\t%f\n\n",gradcor[0],gradcor[1],gradcor[2]);
        scalar_time(ncoord,&nstime,sublagt,cormod,parcor,flagcor,gradcor,flagnuis,grad,npar,nuis,sdata,weigthed,maxtime,ww,mom_cond,dist,coordx,coordy,maxdist);
        
  
        //printf("C: %f\t%f\t%f\t%f\n",mom_cond[0],mom_cond[1],mom_cond[2],mom_cond[3]);
       
        for(kk=0;kk<npar[0];kk++)
        {
            //Rprintf("mom_cond[kk]: %f\n",mom_cond[kk]);
            vector_mean[kk][i]=mom_cond[kk];
        }  //vector means for eac winndow
        Free(mom_cond);
        
    }
    
    Free(sdata);Free(gradcor);Free(ww);
    //Free(grad);
    
    tot= (double **) Calloc(npar[0],double *);
    for(i=0;i<npar[0];i++)tot[i]=(double *) Calloc(nsub,double);
    //mean over blocks
    for(q=0;q<nsub ;q++)
    {
        for(kk=0;kk<*npar;kk++)
        {
            block_mean[kk]= block_mean[kk]+ vector_mean[kk][q]/nsub;
            //printf("block_mean[kk]:%f\n",block_mean[kk]);//problema
            
        }
    }
    //if(block_mean[0]<0 ){return;}
    for(q=0;q<nsub ;q++)
    {
        for(kk=0;kk<*npar;kk++)
        {
            tot[kk][q]=  vector_mean[kk][q]+block_mean[kk];
            //printf("tot[kk][q]: %f\n", tot[kk][q]);//Problema
        }
    }
    ////computing variance covarianmce matrix
    nmat=0.5*npar[0]*(npar[0]-1)+npar[0];
    vv=(double *) Calloc(nsub*nmat,double);
    h=0;
    for(r=0;r<nsub ;r++)
    {
        for(p=0;p< *npar;p++)
        {
            for(q=p;q< *npar;q++)
            {
                vv[h]=vv[h]+(tot[p][r])*(tot[q][r])/(nsub);
                //printf("vv[h]: %f\n", vv[h]); // PROBLEMA
                h++;
            }
        }
    }
    ////building upper triangular of the variance covarianmce matrix (summing over the matrices )
    for(i=0;i<nmat;i++)
    {
        for(q=0;q<nsub ;q++)
        {
            vari[i]=vari[i]+vv[i+q*nmat];
        }
    }
    Free(vv);
    for(i=0;i<npar[0];i++)
    {
        Free(vector_mean[i]);Free(tot[i]);
    }
    Free(vector_mean);
    Free(tot);
    
    
    return;
}*/

/*******************************************************************************************************************************************/

void scalar_spacetime(int *npts,int *nstime,double *sublagt,double *maxtime,int *cormod,double *parcor,int *flagcor, double *gradcor,int *flagnuis, double *grad,int *npar,double *nuis,double *s2data,int *weigthed, double *ww, double *mom_cond,double *maxdist,int *dist, double *scoordx, double *scoordy)
{
    int l=0,t=0,m=0,v=0,kk=0;
    double lagt=0.0,rho=0,lags=0.0;
    for(l=0;l<npts[0];l++)
    {
        for(t=0;t<nstime[0];t++)
        {
            for(m=l;m<npts[0];m++)
            {
                if(l==m)
                {
                    for(v=t+1;v<nstime[0];v++)
                    {
                        lagt=fabs(sublagt[t]-sublagt[v]);
                        if(lagt<=maxtime[0])
                        {
                            //Computing correlation
                            rho=CorFct(cormod,0,lagt,parcor);
                            //Computing the gradient of the corr parameters
                            GradCorrFct(rho,cormod,flagcor,gradcor,0,lagt,parcor);
                            //Compute the gradient of the composite likelihood:
                            Grad_Pair_Gauss(rho,flagnuis,gradcor,grad,npar,nuis,s2data[(t+nstime[0]*l)],s2data[(v+nstime[0]*l)]);
                            /*if(*weigthed)
                             {
                             weights=CorFunBohman(lagt,maxtime[0]);
                             ww[0]=1;ww[1]=1;ww[2]=weights;ww[3]=weights;
                             }
                             else
                             {ww[0]=1;ww[1]=1;ww[2]=1; ww[3]=1;
                             }*/
                            ww[0]=1;ww[1]=1;ww[2]=1; ww[3]=1;
                            for(kk=0;kk<npar[0];kk++)
                            {mom_cond[kk]=mom_cond[kk]+ww[kk]*grad[kk];
                            }
                        }
                    }
                }
                else
                {
                    if(*dist==1) lags=hypot(scoordx[l]-scoordx[m],scoordy[l]-scoordy[m]);
                    // if(*dist==2) lags=Dist_geodesic(scoordx[l],scoordy[l],scoordx[m],scoordy[m]);
                    for(v=0;v<nstime[0];v++)
                    {
                        lagt=fabs(sublagt[t]-sublagt[v]);
                        if(lagt<=maxtime[0] && lags<=maxdist[0])
                        {
                            //Computing correlation
                            rho=CorFct(cormod,lags,lagt,parcor);
                            //Computing the gradient of the corr parameters
                            GradCorrFct(rho,cormod,flagcor,gradcor,lags,lagt,parcor);
                            //Compute the gradient of the composite likelihood:
                            Grad_Pair_Gauss(rho,flagnuis,gradcor,grad,npar,nuis,s2data[(t+nstime[0]*l)],s2data[(v+nstime[0]*m)]);
                            /*if(*weigthed)
                             {
                             weights=CorFunBohman(lags,maxdist[0])*CorFunBohman(lagt,maxtime[0]);
                             ww[0]=1;ww[1]=1;ww[2]=weights;ww[3]=weights;
                             }
                             else
                             {
                             ww[0]=1;ww[1]=1;ww[2]=1; ww[3]=1;
                             }*/
                            ww[0]=1;ww[1]=1;ww[2]=1; ww[3]=1;
                            for(kk=0;kk<npar[0];kk++)
                            {
                                mom_cond[kk]=mom_cond[kk]+ww[kk]*grad[kk];
                            }
                        }
                    }
                }
            }
        }
    }
}


/*******************************************************************************************************************************************/


void SubSamp_spacetime(double *coordx, double *coordy,double *coordt, int *ncoord,int *ntime,int *cormod,double *data,int *dist, double *maxdist,double *maxtime,int *npar,double *parcor,int *nparc,double *nuis,int *nparnuis, int *flagcor, int *flagnuis,double *vari,double *winc, double *winstp,double *winc_t, double *winstp_t,double *block_mean,int *weigthed, int *local_wi, int *dev)


{
    //Rprintf("A: %f %f\n",winc[0],winstp[0]);
    
    double  beta, *rangex, *rangey;
    double *vv,*sdata,*s2data,*xgrid,*ygrid,*scoordx,*scoordy,*sublagt;
    double *gradcor,*grad,*ww;
    double **tot,*mom_cond,**vector_mean;    //matrix moments conditions
    double delta=0.0, deltax=0.0, deltay=0.0, dimwinx=0.0, dimwiny=0.0;
    double winstx=winstp[0], winsty=winstp[0];
    double step=coordt[0]-coordt[1];  // subsampling works in the regular temporal sampling setting
    int *npts, numintx=0, numinty=0,nstime=0;
    int n_win=0,kk=0,h=0,i=0,nsub,nsub_t=0,j=0,f=0,p=0,q=0,r,nmat;
    
    npts=(int *) R_alloc(1, sizeof(int));          //number of location sites in the subwindow
    rangex=(double *) R_alloc(2, sizeof(double));
    rangey=(double *) R_alloc(2, sizeof(double));
    
    scoordx=(double *) Calloc(*ncoord,double);
    scoordy=(double *) Calloc(*ncoord,double);
    sdata=(double *) Calloc(*ncoord * *ntime,double);
    
    Range(coordx,rangex,ncoord[0]);// range of the x-coordinate (min and max)
    Range(coordy,rangey,ncoord[0]);// range of the y-coordinate  (min and max)
    // set the sub-sampling window based on prototype unit window (R_0)
    // and scaling factor (lambda_n)
    deltax=rangex[1]-rangex[0];//  x window length
    deltay=rangey[1]-rangey[0];//  y window length
    if(!winc[0])
    {
        delta=fmin(deltax,deltay);  // I set this rule if no constant
        winc[0]=sqrt(delta)/2;
    }
    if(!winstp[0]) winstp[0]=0.5;   //proportion of the overlapping  0< winstp <=1
    //Rprintf("B: %f %f\n",winc[0],winstp[0]);
    dimwinx=winc[0] * sqrt(deltax);// sub-window x length depends on a constant: deafault??
    dimwiny=winc[0] * sqrt(deltay);// sub-window y length depends on a constant: deafault??
    //Rprintf("C:%f %f %f %f %f %f\n",dimwinx,dimwiny,winc[0],winstp[0],deltax,deltay);
    
    winstx=*winstp * dimwinx;     // x step is a  proportion of sub-window x length (deafult is 0.5)
    winsty=*winstp * dimwiny;     // y step is a  proportion of sub-window y length (deafult is 0.5)
    numintx=floor((deltax-dimwinx)/winstx+1);   //number of overlapping sub-windows is  numintx+1 * numinty+1
    numinty=floor((deltay-dimwiny)/winsty+1);
    //Rprintf("C:%f %f %f %f %f %f %f %f\n",dimwinx,dimwiny,winc[0],winstp[0],deltax,deltay,numintx,numinty);
    
    xgrid=(double *) R_alloc(numintx, sizeof(double));
    ygrid=(double *) R_alloc(numinty, sizeof(double));
    
    SeqStep(rangex,numintx,winstx,xgrid);
    SeqStep(rangey,numinty,winsty,ygrid);
    
    //default sub window temporal length
    if(!(*winc_t))
    {
        beta=CorFct(cormod,0,1,parcor);
        *winc_t=R_pow(2*beta/(1-R_pow(beta,2)),(double) 2/3)*R_pow(*ntime*1.5,(double)1/3);
        if(*winc_t<4*step) *winc_t=2*step;// if the length is too small
        if(*winc_t>=*ntime) *winc_t=*ntime-step;} // if the length is too big
    //set the spatial-temporal windows:
    int wint = (int) *winc_t;
    if(*winstp_t==0) *winstp_t=wint; //defualt for the forward step:the minimum distance
    //else *winstp=*winstp*wint;   //otherwise a proportion of the temporal window
    sublagt=(double *) R_alloc(wint, sizeof(double));
    sublagt[0]=step;
    
    s2data=(double *)  Calloc(*ncoord*wint,double);
    //set the temporal distances for the sub-sample:
    for(i=0;i<wint;i++)
    {
        sublagt[i+1]=sublagt[i]+step;nstime++;
    }
    nsub_t=floor(((*ntime-wint)/(wint*winstp_t[0])+1));
    //Rprintf("%d %d\n",nsub_t,floor(((*ntime-wint)/(winstp_t[0])+1)));
    
    
    // vector of means for each window
    vector_mean= (double **) Calloc(npar[0],double *);
    for(i=0;i<npar[0];i++)
    {
        vector_mean[i]=(double *) Calloc((numintx+1)*(numinty+1)*nsub_t,double);
    }
    
    ww=(double *) Calloc(*npar,double);
    gradcor=(double *) Calloc(*nparc,double);
    grad=(double *) Calloc(*npar,double);
    
    nsub=0;
    n_win=numintx*numinty;
    //Rprintf("%d %d %d %f\n",nsub_t,n_win,nsub_t*n_win,winstp_t[0]);
    for(i=0;i<=numintx;i++)
    {
        for(j=0;j<=numinty;j++)
        {
            // cycle for each block∫∫
            *npts=0;   // number of points in the block
            SetSampling_s(ncoord[0],ntime[0],coordx,coordy,data,npts,scoordx,scoordy,sdata,xgrid[i]+dimwinx,xgrid[i],ygrid[j]+dimwiny,ygrid[j]);//  create data and coordinates of the sub-windows
            //  scoordx    scoordy are the spatial coords of the subwindow
            // sdata the associated data
            // npts number of loc sites in the window
            if( (((xgrid[i]+dimwinx)<rangex[1])||is_equal(xgrid[i]+dimwinx,rangex[1]))&&     //le mezze finestre sono escluse.
               (((ygrid[j]+dimwiny))<rangey[1]||is_equal(ygrid[j]+dimwiny,rangey[1]))&&       //Anche le finestre che hanno "pochi" punti(griglia irregolare)
               (npts[0]>5))
            {
                for(f=0;f<nsub_t;f++)
                {//loop for the number of sub-sampling:
                    // set the sub-sample of the data:
                    SetSampling_t(sdata,s2data,npts[0],*ntime,wint,f);
                    mom_cond=(double *) Calloc(*npar,double);
                    /******************************************/
                    /*computing gradient in the window*/
                    /******************************************/
                    
                    scalar_spacetime(npts,&nstime,sublagt,maxtime,cormod,parcor,flagcor,gradcor,flagnuis,grad,npar,nuis,s2data,weigthed,ww,mom_cond,maxdist,dist,scoordx,scoordy);
                    
                    
                    /******************************************/
                    /******************************************/
                    for(kk=0;kk<npar[0];kk++)
                    {
                        vector_mean[kk][nsub]=mom_cond[kk];
                    }  //vector means for eac winndow
                    nsub=nsub+1;
                    Free(mom_cond);
                }
            }
        }
    }
    
    Free(scoordx);Free(scoordy);Free(sdata);Free(s2data);Free(grad);Free(gradcor);Free(ww);
    tot= (double **) Calloc(npar[0],double *);
    for(i=0;i<npar[0];i++)
    {
        tot[i]=(double *) Calloc(nsub,double);
    }
    
    //mean over blocks
    for(q=0;q<nsub;q++)
    {
        for(kk=0;kk<*npar;kk++)
        {
            block_mean[kk]= block_mean[kk]+ vector_mean[kk][q]/nsub;
        }
    }
    for(q=0;q<nsub ;q++)
    {
        for(kk=0;kk<*npar;kk++)
        {
            tot[kk][q]=  vector_mean[kk][q]+block_mean[kk];
        }
    }
    ////computing variance covarianmce matrix
    nmat=0.5*npar[0]*(npar[0]-1)+npar[0];
    vv=(double *) Calloc(nsub*nmat,double);
    h=0;
    for(r=0;r<nsub ;r++)
    {
        for(p=0;p< *npar;p++)
        {
            for(q=p;q< *npar;q++)
            {
                vv[h]=vv[h]+(tot[p][r])*(tot[q][r])/(nsub);
                h++;
            }
        }
    }
    ////building upper triangular of the variance covarianmce matrix (summing over the matrices )
    for(i=0;i<nmat;i++)
    {
        for(q=0;q<nsub ;q++)
        {
            vari[i]=vari[i]+vv[i+q*nmat];
        }
    }
    Free(vv);
    for(i=0;i<npar[0];i++)
    {
        Free(vector_mean[i]);Free(tot[i]);
    }
    Free(vector_mean);
    Free(tot);
    
    return;
}

/*******************************************************************************************************************************************/
