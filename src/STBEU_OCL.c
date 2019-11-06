#include "header.h"

void SubSamp_space_DE_ocl(double *coordx, double *coordy,double *coordt, int *ncoord,int *ntime,int *cormod,double *data,int *dist, double *maxdist,double *maxtime,int *npar,double *parcor,int *nparc,double *nuis,int *nparnuis, int *flagcor, int *flagnuis,double *vari,double *winc, double *winstp,double *a, double *b,double *block_mean,int *weigthed, int *local_wi, int *dev)
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
    
    if(!winc[0])
    {
        delta=fmin(deltax,deltay);  // I set this rule if no constant
        winc[0]=sqrt(delta)/2;
    }
    if(!winstp[0]) winstp[0]=0.5;   //proportion of the overlapping  0< winstp <=1, If winstp=0 then  winstp=0.5
    
    dimwinx=winc[0] * sqrt(deltax);// sub-window x length depends on a constant: deafault??
    dimwiny=winc[0] * sqrt(deltay);// sub-window y length depends on a constant: deafault??
    
    winstx=*winstp * dimwinx;     // x step is a  proportion of sub-window x length (deafult is 0.5)
    winsty=*winstp * dimwiny;     // y step is a  proportion of sub-window y length (deafult is 0.5)
    numintx=floor((deltax-dimwinx)/winstx+1);   //number of overlapping sub-windows is  numintx+1 * numinty+1
    numinty=floor((deltay-dimwiny)/winsty+1);
    
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
    n_win=numintx*numinty;   //// number of windows
    
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
                //-----------------------------------------
                //computing gradient in the window
                //-----------------------------------------
   
                DoubleExpOCL(npts,ntime,coordt,maxtime,maxdist,cormod,parcor,flagcor,flagnuis,npar,nuis,sdata,weigthed,mom_cond,dist,scoordx,scoordy,gradcor,grad,ww,local_wi,dev);
                
                //-----------------------------------------
                //-----------------------------------------
                for(kk=0;kk<npar[0];kk++)
                {
                    vector_mean[kk][nsub]=mom_cond[kk];
                }  //vector means for each winndow
                nsub=nsub+1;  // counting number of windows
                Free(mom_cond);
            }
        }
    }
    Free(scoordx);Free(scoordy);Free(sdata);Free(grad);Free(gradcor);Free(ww);
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



void SubSamp_space_GN_ocl(double *coordx, double *coordy,double *coordt, int *ncoord,int *ntime,int *cormod,double *data,int *dist, double *maxdist,double *maxtime,int *npar,double *parcor,int *nparc,double *nuis,int *nparnuis, int *flagcor, int *flagnuis,double *vari,double *winc, double *winstp,double *a, double *b,double *block_mean,int *weigthed, int *local_wi, int *dev)
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
    
    if(!winc[0])
    {
        delta=fmin(deltax,deltay);  // I set this rule if no constant
        winc[0]=sqrt(delta)/2;
    }
    if(!winstp[0]) winstp[0]=0.5;   //proportion of the overlapping  0< winstp <=1, If winstp=0 then  winstp=0.5

    dimwinx=winc[0] * sqrt(deltax);// sub-window x length depends on a constant: deafault??
    dimwiny=winc[0] * sqrt(deltay);// sub-window y length depends on a constant: deafault??
    
    winstx=*winstp * dimwinx;     // x step is a  proportion of sub-window x length (deafult is 0.5)
    winsty=*winstp * dimwiny;     // y step is a  proportion of sub-window y length (deafult is 0.5)
    numintx=floor((deltax-dimwinx)/winstx+1);   //number of overlapping sub-windows is  numintx+1 * numinty+1
    numinty=floor((deltay-dimwiny)/winsty+1);
    
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
    n_win=numintx*numinty;   //// number of windows
    
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
                //-----------------------------------------
                //computing gradient in the window
                //-----------------------------------------


                GneitingOCL(npts,ntime,coordt,maxtime,maxdist,cormod,parcor,flagcor,flagnuis,npar,nuis,sdata,weigthed,mom_cond,dist,scoordx,scoordy,gradcor,grad,ww,local_wi,dev);
                
                //-----------------------------------------
                //-----------------------------------------
                for(kk=0;kk<npar[0];kk++)
                {
                    vector_mean[kk][nsub]=mom_cond[kk];
                }  //vector means for each winndow
                nsub=nsub+1;  // counting number of windows
                Free(mom_cond);
            }
        }
    }
    Free(scoordx);Free(scoordy);Free(sdata);Free(grad);Free(gradcor);Free(ww);
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

void SubSamp_space_WE_ocl(double *coordx, double *coordy,double *coordt, int *ncoord,int *ntime,int *cormod,double *data,int *dist, double *maxdist,double *maxtime,int *npar,double *parcor,int *nparc,double *nuis,int *nparnuis, int *flagcor, int *flagnuis,double *vari,double *winc, double *winstp,double *a, double *b,double *block_mean,int *weigthed, int *local_wi, int *dev)
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
    
    if(!winc[0])
    {
        delta=fmin(deltax,deltay);  // I set this rule if no constant
        winc[0]=sqrt(delta)/2;
    }
    if(!winstp[0]) winstp[0]=0.5;   //proportion of the overlapping  0< winstp <=1, If winstp=0 then  winstp=0.5

    dimwinx=winc[0] * sqrt(deltax);// sub-window x length depends on a constant: deafault??
    dimwiny=winc[0] * sqrt(deltay);// sub-window y length depends on a constant: deafault??
    
    winstx=*winstp * dimwinx;     // x step is a  proportion of sub-window x length (deafult is 0.5)
    winsty=*winstp * dimwiny;     // y step is a  proportion of sub-window y length (deafult is 0.5)
    numintx=floor((deltax-dimwinx)/winstx+1);   //number of overlapping sub-windows is  numintx+1 * numinty+1
    numinty=floor((deltay-dimwiny)/winsty+1);
    
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
    n_win=numintx*numinty;   //// number of windows
    
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
                WendOCL(npts,ntime,coordt,maxtime,maxdist,cormod,parcor,flagcor,flagnuis,npar,nuis,sdata,weigthed,mom_cond,dist,scoordx,scoordy,gradcor,grad,ww,local_wi,dev);
               
                /******************************************/
                /******************************************/
                for(kk=0;kk<npar[0];kk++)
                {
                    vector_mean[kk][nsub]=mom_cond[kk];
                }  //vector means for each winndow
                nsub=nsub+1;  // counting number of windows
                Free(mom_cond);
            }
        }
    }
    Free(scoordx);Free(scoordy);Free(sdata);Free(grad);Free(gradcor);Free(ww);
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

void SubSamp_time_DE_ocl(double *coordx, double *coordy,double *coordt, int *ncoord,int *ntime,int *cormod,double *data,int *dist, double *maxdist,double *maxtime,int *npar,double *parcor,int *nparc,double *nuis,int *nparnuis, int *flagcor, int *flagnuis,double *vari,double *a, double *b,double *winc, double *winstp,double *block_mean,int *weigthed, int *local_wi, int *dev)
{
    double beta, *gradcor;
    double *vv,*ww,*sdata, *grad, *sublagt;
    double **tot,*mom_cond,**vector_mean;    //matrix moments conditions
    double step=coordt[0]-coordt[1];  // subsampling works in the regular temporal sampling setting
    int nsub=0, nstime=0;
    int kk=0,h=0,i=0,p=0,q=0,r,nmat;
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
    //SetSampling_t(data,sdata,*ncoord,*ntime,wint,0);
    //for(i=0;i<*ncoord* *ntime;i++){printf("sdata[%d]: %f\n",i,sdata[i]);}
    
    for(i=0;i<nsub;i++)
    {//loop for the number of sub-sampling:
        mom_cond=(double *) Calloc(*npar,double);
        // set the sub-sample of the data:
        SetSampling_t(data,sdata,*ncoord,*ntime,wint,i);
        
        
        /******************************************/
        /*computing gradient in the window*/
        /******************************************/
        //scalar_time(ncoord,&nstime,sublagt,cormod,parcor,flagcor,gradcor,flagnuis,grad,npar,nuis,sdata,weigthed,maxtime,ww,mom_cond,dist,coordx,coordy,maxdist);
        DoubleExpOCL(ncoord,&nstime,sublagt,maxtime,maxdist,cormod,parcor,flagcor,flagnuis,npar,nuis,sdata,weigthed,mom_cond,dist,coordx,coordy,gradcor,grad,ww,local_wi,dev);
        
        /******************************************/
        /******************************************/
        for(kk=0;kk<npar[0];kk++)
        {
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
            block_mean[kk]= block_mean[kk]+ vector_mean[kk][q]/nsub;;
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
                vv[h]=vv[h]+(tot[p][r])*(tot[q][r])/(nsub);h++;
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

void SubSamp_time_GN_ocl(double *coordx, double *coordy,double *coordt, int *ncoord,int *ntime,int *cormod,double *data,int *dist, double *maxdist,double *maxtime,int *npar,double *parcor,int *nparc,double *nuis,int *nparnuis, int *flagcor, int *flagnuis,double *vari,double *a, double *b,double *winc, double *winstp,double *block_mean,int *weigthed, int *local_wi, int *dev)
{
    double beta, *gradcor;
    double *vv,*ww,*sdata, *grad, *sublagt;
    double **tot,*mom_cond,**vector_mean;    //matrix moments conditions
    double step=coordt[0]-coordt[1];  // subsampling works in the regular temporal sampling setting
    int nsub=0, nstime=0;
    int kk=0,h=0,i=0,p=0,q=0,r,nmat;
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
    //SetSampling_t(data,sdata,*ncoord,*ntime,wint,0);
    //for(i=0;i<*ncoord* *ntime;i++){printf("sdata[%d]: %f\n",i,sdata[i]);}

    for(i=0;i<nsub;i++)
    {//loop for the number of sub-sampling:
        mom_cond=(double *) Calloc(*npar,double);
        // set the sub-sample of the data:
        SetSampling_t(data,sdata,*ncoord,*ntime,wint,i);
        
        
        /******************************************/
        /*computing gradient in the window*/
        /******************************************/
        //scalar_time(ncoord,&nstime,sublagt,cormod,parcor,flagcor,gradcor,flagnuis,grad,npar,nuis,sdata,weigthed,maxtime,ww,mom_cond,dist,coordx,coordy,maxdist);
        GneitingOCL(ncoord,&nstime,sublagt,maxtime,maxdist,cormod,parcor,flagcor,flagnuis,npar,nuis,sdata,weigthed,mom_cond,dist,coordx,coordy,gradcor,grad,ww,local_wi,dev);
        
        /******************************************/
        /******************************************/
        for(kk=0;kk<npar[0];kk++)
        {
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
            block_mean[kk]= block_mean[kk]+ vector_mean[kk][q]/nsub;;
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
                vv[h]=vv[h]+(tot[p][r])*(tot[q][r])/(nsub);h++;
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

void SubSamp_time_WE_ocl(double *coordx, double *coordy,double *coordt, int *ncoord,int *ntime,int *cormod,double *data,int *dist, double *maxdist,double *maxtime,int *npar,double *parcor,int *nparc,double *nuis,int *nparnuis, int *flagcor, int *flagnuis,double *vari,double *a, double *b,double *winc, double *winstp,double *block_mean,int *weigthed, int *local_wi, int *dev)
{
    double beta, *gradcor;
    double *vv,*ww,*sdata, *grad, *sublagt;
    double **tot,*mom_cond,**vector_mean;    //matrix moments conditions
    double step=coordt[0]-coordt[1];  // subsampling works in the regular temporal sampling setting
    int nsub=0, nstime=0;
    int kk=0,h=0,i=0,p=0,q=0,r,nmat;
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
        //scalar_time(ncoord,&nstime,sublagt,cormod,parcor,flagcor,gradcor,flagnuis,grad,npar,nuis,sdata,weigthed,maxtime,ww,mom_cond,dist,coordx,coordy,maxdist);
        WendOCL(ncoord,&nstime,sublagt,maxtime,maxdist,cormod,parcor,flagcor,flagnuis,npar,nuis,sdata,weigthed,mom_cond,dist,coordx,coordy,gradcor,grad,ww,local_wi,dev);
        
        /******************************************/
        /******************************************/
        for(kk=0;kk<npar[0];kk++)
        {
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
                vv[h]=vv[h]+(tot[p][r])*(tot[q][r])/(nsub);h++;
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

void SubSamp_spacetime_DE_ocl(double *coordx, double *coordy,double *coordt, int *ncoord,int *ntime,int *cormod,double *data,int *dist, double *maxdist,double *maxtime,int *npar,double *parcor,int *nparc,double *nuis,int *nparnuis, int *flagcor, int *flagnuis,double *vari,double *winc, double *winstp,double *winc_t, double *winstp_t,double *block_mean,int *weigthed, int *local_wi, int *dev)


{
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
    dimwinx=winc[0] * sqrt(deltax);// sub-window x length depends on a constant: deafault??
    dimwiny=winc[0] * sqrt(deltay);// sub-window y length depends on a constant: deafault??
    
    winstx=*winstp * dimwinx;     // x step is a  proportion of sub-window x length (deafult is 0.5)
    winsty=*winstp * dimwiny;     // y step is a  proportion of sub-window y length (deafult is 0.5)
    numintx=floor((deltax-dimwinx)/winstx+1);   //number of overlapping sub-windows is  numintx+1 * numinty+1
    numinty=floor((deltay-dimwiny)/winsty+1);
    
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
                    
                    DoubleExpOCL(npts,&nstime,sublagt,maxtime,maxdist,cormod,parcor,flagcor,flagnuis,npar,nuis,s2data,weigthed,mom_cond,dist,scoordx,scoordy,gradcor,grad,ww,local_wi,dev);
                    
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

void SubSamp_spacetime_GN_ocl(double *coordx, double *coordy,double *coordt, int *ncoord,int *ntime,int *cormod,double *data,int *dist, double *maxdist,double *maxtime,int *npar,double *parcor,int *nparc,double *nuis,int *nparnuis, int *flagcor, int *flagnuis,double *vari,double *winc, double *winstp,double *winc_t, double *winstp_t,double *block_mean,int *weigthed, int *local_wi, int *dev)


{
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
    
    dimwinx=winc[0] * sqrt(deltax);// sub-window x length depends on a constant: deafault??
    dimwiny=winc[0] * sqrt(deltay);// sub-window y length depends on a constant: deafault??
    
    winstx=*winstp * dimwinx;     // x step is a  proportion of sub-window x length (deafult is 0.5)
    winsty=*winstp * dimwiny;     // y step is a  proportion of sub-window y length (deafult is 0.5)
    numintx=floor((deltax-dimwinx)/winstx+1);   //number of overlapping sub-windows is  numintx+1 * numinty+1
    numinty=floor((deltay-dimwiny)/winsty+1);
    
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
                    
                    GneitingOCL(npts,&nstime,sublagt,maxtime,maxdist,cormod,parcor,flagcor,flagnuis,npar,nuis,s2data,weigthed,mom_cond,dist,scoordx,scoordy,gradcor,grad,ww,local_wi,dev);
                    
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

void SubSamp_spacetime_WE_ocl(double *coordx, double *coordy,double *coordt, int *ncoord,int *ntime,int *cormod,double *data,int *dist, double *maxdist,double *maxtime,int *npar,double *parcor,int *nparc,double *nuis,int *nparnuis, int *flagcor, int *flagnuis,double *vari,double *winc, double *winstp,double *winc_t, double *winstp_t,double *block_mean,int *weigthed, int *local_wi, int *dev)


{
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
    
    dimwinx=winc[0] * sqrt(deltax);// sub-window x length depends on a constant: deafault??
    dimwiny=winc[0] * sqrt(deltay);// sub-window y length depends on a constant: deafault??
    
    winstx=*winstp * dimwinx;     // x step is a  proportion of sub-window x length (deafult is 0.5)
    winsty=*winstp * dimwiny;     // y step is a  proportion of sub-window y length (deafult is 0.5)
    numintx=floor((deltax-dimwinx)/winstx+1);   //number of overlapping sub-windows is  numintx+1 * numinty+1
    numinty=floor((deltay-dimwiny)/winsty+1);
    
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
                    
                    WendOCL(npts,&nstime,sublagt,maxtime,maxdist,cormod,parcor,flagcor,flagnuis,npar,nuis,s2data,weigthed,mom_cond,dist,scoordx,scoordy,gradcor,grad,ww,local_wi,dev);
                    
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




void DoubleExpOCL(int *npts, int *ntime,double *coordt, double *maxtime,double *maxdist,int *cormod, double *parcor,int *flagcor,int *flagnuis,int *npar,double *nuis,double *data,int *weigthed, double *mom_cond, int *dist, double *coordx,double *coordy,double *gradcor,double  *grad, double *ww,int *local_wi, int *dev)
{
    
    
    
    // Setting data and params for sclar_space function
    double *h_x,*h_y,*h_data, *h_t;
    
    h_x = (double *) R_alloc(npts[0], sizeof(double));
    h_y = (double *) R_alloc(npts[0], sizeof(double));
    h_data = (double *) R_alloc(npts[0]*ntime[0], sizeof(double));
    h_t = (double *) R_alloc(ntime[0], sizeof(double));
    
    h_x = coordx;
    h_y = coordy;
    h_data = data;
    h_t = coordt;
    int i;
    int nparc=0, nparnuis=0;
    int *int_par;
    double *dou_par; // objects to pass int and double parametes to the kernel function
    //int_par = (int*)calloc(1+1+2+3+1+1+1+1+1, sizeof(int));//length sum of: npts+ntime+flagcor+flagnuis+npar+weigthed+dist+nparc+nparnuis+cormod
    //dou_par = (double*)calloc(1+1+2+3, sizeof(int));//length sum of: maxtime+maxdist+parcor+nuis
    
    
    
    for(i =0; i<2;i++){nparc += flagcor[i];}
    for(i =0; i<3;i++){nparnuis += flagnuis[i];}
    int_par = (int*)calloc((13), sizeof(int));
    dou_par = (double*)calloc((7), sizeof(double));
    
    int_par[0] = npts[0];
    int_par[1] = ntime[0];
    int_par[2] = flagcor[0];
    int_par[3] = flagcor[1];
    int_par[4] = flagnuis[0];
    int_par[5] = flagnuis[1];
    int_par[6] = flagnuis[2];
    int_par[7] = npar[0];//
    int_par[8] = weigthed[0];
    int_par[9] = dist[0];
    int_par[10] = nparc;
    int_par[11] = nparnuis;
    int_par[12] = cormod[0];
    
    
    
    dou_par[0] = maxtime[0];
    dou_par[1] = maxdist[0];
    dou_par[2] = parcor[0];
    dou_par[3] = parcor[1];
    dou_par[4] = nuis[0];
    dou_par[5] = nuis[1];
    dou_par[6] = nuis[2];
    
    
    
    cl_int err;
    cl_device_id        device;     // compute device id
    cl_context       context;       // compute context
    cl_command_queue commands;      // compute command queue
    cl_program       program;       // compute program
    cl_kernel        kernel;     // compute kernel
    
    
    // Vars for querying Device Info:
    
    char* value;
    size_t valueSize;
    cl_uint maxComputeUnits;
    cl_ulong long_entries;
    size_t p_size;
    
    // Set up OpenCL context, queue, kernel, etc.
    
    cl_uint deviceIndex = dev[0];
    
    // Get list of devices
    cl_device_id devices[MAX_DEVICES];
    unsigned numDevices = getDeviceList(devices);
    
    // Check device index in range
    if (deviceIndex >= numDevices)
    {
        //printf("Invalid device index (try '--list')\n");
        //return EXIT_FAILURE;
    }
    
    device = devices[deviceIndex];
    
    char name[MAX_INFO_STRING];
    getDeviceName(device, name);
    
    
    // Get Device Info for Execution Model:
    
    // print hardware device version
    clGetDeviceInfo(device, CL_DEVICE_VERSION, 0, NULL, &valueSize);
    value = (char*) malloc(valueSize);
    clGetDeviceInfo(device, CL_DEVICE_VERSION, valueSize, value, NULL);
    //printf(" %d.%d Hardware version: %s\n", i+1, 1, value);
    free(value);
    
    // print software driver version
    clGetDeviceInfo(device, CL_DRIVER_VERSION, 0, NULL, &valueSize);
    value = (char*) malloc(valueSize);
    clGetDeviceInfo(device, CL_DRIVER_VERSION, valueSize, value, NULL);
    //printf(" %d.%d Software version: %s\n", i+1, 2, value);
    free(value);
    
    // print c version supported by compiler for device
    clGetDeviceInfo(device, CL_DEVICE_OPENCL_C_VERSION, 0, NULL, &valueSize);
    value = (char*) malloc(valueSize);
    clGetDeviceInfo(device, CL_DEVICE_OPENCL_C_VERSION, valueSize, value, NULL);
    //printf(" %d.%d OpenCL C version: %s\n", i+1, 3, value);
    free(value);
    
    // print parallel compute units
    clGetDeviceInfo(device, CL_DEVICE_MAX_COMPUTE_UNITS,
                    sizeof(maxComputeUnits), &maxComputeUnits, NULL);
    //printf(" %d.%d Parallel compute units: %d\n", i+1, 4, maxComputeUnits);
    
    
    clGetDeviceInfo(device,CL_DEVICE_GLOBAL_MEM_SIZE,sizeof(cl_ulong),&long_entries,NULL);
    //printf(" %d.%d Global Memory (MB):\t%llu\n",i+1, 5,long_entries/1024/1024);
    clGetDeviceInfo(device,CL_DEVICE_GLOBAL_MEM_CACHE_SIZE,sizeof(cl_ulong),&long_entries,NULL);
    //printf(" %d.%d Global Memory Cache (MB):\t%llu\n",i+1, 6,long_entries/1024/1024);
    clGetDeviceInfo(device,CL_DEVICE_LOCAL_MEM_SIZE,sizeof(cl_ulong),&long_entries,NULL);
    //printf(" %d.%d Local Memory (KB):\t%llu\n",i+1, 7,long_entries/1024);
    clGetDeviceInfo(device,CL_DEVICE_MAX_CLOCK_FREQUENCY,sizeof(cl_ulong),&long_entries,NULL);
    //printf(" %d.%d Max clock (MHz) :\t%llu\n",i+1, 8,long_entries);
    clGetDeviceInfo(device,CL_DEVICE_MAX_WORK_GROUP_SIZE,sizeof(size_t),&p_size,NULL);
    //printf(" %d.%d Max Work Group Size:\t%zu\n",i+1, 9,p_size);
    clGetDeviceInfo(device,CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS,sizeof(size_t),&p_size,NULL);
    //printf(" %d.%d Maximum dimensions:\t%zu\n",i+1, 10,p_size);
    clGetDeviceInfo(device,CL_DEVICE_MAX_MEM_ALLOC_SIZE,sizeof(size_t),&p_size,NULL);
    //printf(" %d.%d Max device buffer size (MB):\t%zu\n",i+1, 11,p_size/1024/1024);
    
    
    
    // Create a compute context
    context = clCreateContext(0, 1, &device, NULL, NULL, &err);
    checkError(err, "Creating context");
    
    // Create a command queue
    commands = clCreateCommandQueue(context, device, 0, &err);
    checkError(err, "Creating command queue");
    // Create the compute program from the source buffer
    //char *kernelsource = getKernelSource("DouExp.cl");
    //program = clCreateProgramWithSource(context, 1, (const char **) & kernelsource, NULL, &err);
    //checkError(err, "Creating program");
    
    
    FILE *fp;
    size_t binary_size;
    char *binary_buf;
    cl_int binary_status;
    
    fp = fopen("DouExp.cl.bin", "r");
    if (!fp) {
        Rprintf("Failed to load kernel.\n");
    }
    binary_buf = (char *)malloc(MAX_BINARY_SIZE);
    binary_size = fread(binary_buf, 1, MAX_BINARY_SIZE, fp);
    fclose(fp);
    
    program = clCreateProgramWithBinary(
                                        context, 1, &device, (const size_t *)&binary_size,
                                        (const unsigned char **)&binary_buf, &binary_status, &err
                                        );
    free(binary_buf);
    
    
    
    
    // Build the program
    err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
    if (err != CL_SUCCESS)
    {
        size_t len;
        char buffer[2048*4];
        
        Rprintf("Error: Failed to build program executable!\n%s\n", err_code(err));
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
        SEP;
        Rprintf("Build Log:\n%s\n", buffer);
        SEP;
        //return EXIT_FAILURE;
    }
    // Create the compute kernel from the program
    
    kernel = clCreateKernel(program, "scalarspaceocl", &err);
    checkError(err, "Creating kernel");
    
    // Find kernel work-group size
    size_t work_group_size;
    err = clGetKernelWorkGroupInfo (kernel, device, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &work_group_size, NULL);
    //checkError(err, "Getting kernel work group info");
    //printf("Maximum Kernel work-group size: %lu\n", work_group_size);
    
    
    // Creating buffers
    size_t coords_buff = sizeof(double) * npts[0];
    size_t coordt_buff = sizeof(double) * ntime[0];
    size_t data_buff = sizeof(double) * npts[0]*ntime[0];
    size_t int_par_buff = sizeof(int) * (13);
    size_t dou_par_buff = sizeof(double) * (7);
    
    size_t g1,g2;
    const int ll1 =local_wi[0];
    const int ll2 =local_wi[1];
    g1 = npts[0] + (ll1 - (npts[0] & (ll1-1))); // SPACE
    g2 = ntime[0] + (ll2 - (ntime[0] & (ll2-1))); //TIME
    
    
    //printf("GLOBAL:\t%zu\t%zu\n",g1,g2);
    size_t local[2] = {ll1,ll2};
    size_t global[2] = { g1,g2};
    
    int length = g1*g2;
    //printf("LENGTH: %d\n",length);
    size_t length_buff = sizeof(double)* (length);
    
    double *h_mom_cond0,*h_mom_cond1,*h_mom_cond2,*h_mom_cond3;
    h_mom_cond0= (double*)calloc(length, sizeof(double));
    h_mom_cond1= (double*)calloc(length, sizeof(double));
    h_mom_cond2= (double*)calloc(length, sizeof(double));
    h_mom_cond3= (double*)calloc(length, sizeof(double));
    
   
    
    //st_alloc = wtime();
    cl_mem d_x = clCreateBuffer(context, CL_MEM_READ_ONLY, coords_buff, NULL, &err);
    checkError(err, "Creating buffer device X");
    err = clEnqueueWriteBuffer(commands, d_x, CL_TRUE, 0, coords_buff, (void*)h_x, 0, NULL, NULL);
    checkError(err, "Writing buffer device X");
    
    
    cl_mem d_y = clCreateBuffer(context, CL_MEM_READ_ONLY, coords_buff, NULL, &err);
    err = clEnqueueWriteBuffer(commands, d_y, CL_TRUE, 0, coords_buff, (void*)h_y, 0, NULL, NULL);
    checkError(err, "Creating buffer device Y");
    
    
    cl_mem d_data = clCreateBuffer(context, CL_MEM_READ_ONLY, data_buff, NULL, &err);
    err = clEnqueueWriteBuffer(commands, d_data, CL_TRUE, 0, data_buff, (void*)h_data, 0, NULL, NULL);
    checkError(err, "Creating buffer device data");
    
    cl_mem d_t = clCreateBuffer(context, CL_MEM_READ_ONLY, coordt_buff, NULL, &err);
    err = clEnqueueWriteBuffer(commands, d_t, CL_TRUE, 0, coordt_buff, (void*)h_t, 0, NULL, NULL);
    checkError(err, "Creating buffer device time");
    
    cl_mem m0_sol = clCreateBuffer(context, CL_MEM_WRITE_ONLY, length_buff, NULL, &err);
    checkError(err, "Creating buffer device m0_sol");
    err = clEnqueueWriteBuffer(commands, m0_sol, CL_TRUE, 0, length_buff, (void*)h_mom_cond0, 0, NULL, NULL);
    checkError(err, "Writing buffer device m0_sol");
    
    
    cl_mem m1_sol = clCreateBuffer(context, CL_MEM_WRITE_ONLY, length_buff, NULL, &err);
    checkError(err, "Creating buffer device m1_sol");
    err = clEnqueueWriteBuffer(commands, m0_sol, CL_TRUE, 0, length_buff, (void*)h_mom_cond1, 0, NULL, NULL);
    checkError(err, "Writing buffer device m1_sol");
    
    cl_mem m2_sol = clCreateBuffer(context, CL_MEM_WRITE_ONLY, length_buff, NULL, &err);
    checkError(err, "Creating buffer device m2_sol");
    err = clEnqueueWriteBuffer(commands, m0_sol, CL_TRUE, 0, length_buff, (void*)h_mom_cond2, 0, NULL, NULL);
    checkError(err, "Writing buffer device m2_sol");
    
    cl_mem m3_sol = clCreateBuffer(context, CL_MEM_WRITE_ONLY, length_buff, NULL, &err);
    checkError(err, "Creating buffer device m3_sol");
    err = clEnqueueWriteBuffer(commands, m0_sol, CL_TRUE, 0, length_buff, (void*)h_mom_cond3, 0, NULL, NULL);
    checkError(err, "Writing buffer device m3_sol");
    
    cl_mem d_int_par = clCreateBuffer(context, CL_MEM_READ_ONLY, int_par_buff, NULL, &err);
    err = clEnqueueWriteBuffer(commands, d_int_par, CL_TRUE, 0, int_par_buff, (void*)int_par, 0, NULL, NULL);
    checkError(err, "Creating buffer device int params");
    
    
    cl_mem d_dou_par = clCreateBuffer(context, CL_MEM_READ_ONLY, dou_par_buff, NULL, &err);
    err = clEnqueueWriteBuffer(commands, d_dou_par, CL_TRUE, 0, dou_par_buff, (void*)dou_par, 0, NULL, NULL);
    checkError(err, "Creating buffer device double params");
    
    
    // Push the data out to device
    clFinish(commands);
    
    err = clSetKernelArg(kernel,  0, sizeof(cl_mem), &d_t); //coordt
    err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &d_x); //coordx
    err |= clSetKernelArg(kernel, 2, sizeof(cl_mem), &d_y); //coordx
    err |= clSetKernelArg(kernel, 3, sizeof(cl_mem), &d_data); //data
    err |= clSetKernelArg(kernel, 4, sizeof(cl_mem), &m0_sol); //mom_cond0
    err |= clSetKernelArg(kernel, 5, sizeof(cl_mem), &m1_sol); //mom_cond1
    err |= clSetKernelArg(kernel, 6, sizeof(cl_mem), &m2_sol); //mom_cond2
    err |= clSetKernelArg(kernel, 7, sizeof(cl_mem), &m3_sol); //mom_cond3
    err |= clSetKernelArg(kernel, 8, sizeof(cl_mem), &d_int_par); //mom_cond3
    err |= clSetKernelArg(kernel, 9, sizeof(cl_mem), &d_dou_par); //mom_cond3
    checkError(err, "Setting kernel args length");
    
    
    // Execute the kernel
    
    err = clEnqueueNDRangeKernel(commands, kernel, 2, NULL, global,local, 0, NULL, NULL);
    checkError(err,"clEnqueueNDRangeKernel\n");
    clFinish(commands);
    
    
    err = clEnqueueReadBuffer(commands, m0_sol, CL_TRUE, 0,length_buff, h_mom_cond0, 0, NULL, NULL);
    err = clEnqueueReadBuffer(commands, m1_sol, CL_TRUE, 0,length_buff, h_mom_cond1, 0, NULL, NULL);
    err = clEnqueueReadBuffer(commands, m2_sol, CL_TRUE, 0,length_buff, h_mom_cond2, 0, NULL, NULL);
    err = clEnqueueReadBuffer(commands, m3_sol, CL_TRUE, 0,length_buff, h_mom_cond3, 0, NULL, NULL);
    clFinish(commands);
    
    
    
    double m0=0,m1=0,m2=0,m3=0;
    for(i=0; i<(npts[0]*ntime[0]); i++)
    {
        m0 += h_mom_cond0[i];
        m1 += h_mom_cond1[i];
        m2 += h_mom_cond2[i];
        m3 += h_mom_cond3[i];
    }
    mom_cond[0] =m0;mom_cond[1] =m1;mom_cond[2] =m2;mom_cond[3] =m3;
    //printf("final result: %.4f\t%.4f\t%.4f\t%.4f\n", m0,m1,m2,m3);
    
    
    // clean up inside kernels
    
    clReleaseProgram(program);
    clReleaseKernel(kernel);
    clReleaseCommandQueue(commands);
    clReleaseContext(context);
    
    
    clReleaseMemObject(d_x);
    clReleaseMemObject(d_y);
    clReleaseMemObject(d_data);
    clReleaseMemObject(d_t);
    clReleaseMemObject(m0_sol);
    clReleaseMemObject(m1_sol);
    clReleaseMemObject(m2_sol);
    clReleaseMemObject(m3_sol);
    clReleaseMemObject(d_int_par);
    clReleaseMemObject(d_dou_par);
    
    
    free(h_mom_cond0);
    free(h_mom_cond1);
    free(h_mom_cond2);
    free(h_mom_cond3);
    
    
}

//------------------------------------------------------------------------------

void GneitingOCL(int *npts, int *ntime,double *coordt, double *maxtime,double *maxdist,int *cormod, double *parcor,int *flagcor,int *flagnuis,int *npar,double *nuis,double *data,int *weigthed, double *mom_cond, int *dist, double *coordx,double *coordy,double *gradcor,double  *grad, double *ww,int *local_wi, int *dev)
{
    
    
    
    // Setting data and params for sclar_space function
    double *h_x,*h_y,*h_data, *h_t;
    h_x = (double *) R_alloc(npts[0], sizeof(double));
    h_y = (double *) R_alloc(npts[0], sizeof(double));
    h_data = (double *) R_alloc(npts[0]*ntime[0], sizeof(double));
    h_t = (double *) R_alloc(ntime[0], sizeof(double));
    
    h_x = coordx;
    h_y = coordy;
    h_data = data;
    h_t = coordt;
    
    
    int i;
    int nparc=0, nparnuis=0;
    int *int_par;
    double *dou_par; // objects to pass int and double parametes to the kernel function
    //int_par = (int*)calloc(1+1+2+3+1+1+1+1+1, sizeof(int));//length sum of: npts+ntime+flagcor+flagnuis+npar+weigthed+dist+nparc+nparnuis+cormod
    //dou_par = (double*)calloc(1+1+2+3, sizeof(int));//length sum of: maxtime+maxdist+parcor+nuis
    
    
    
    for(i =0; i<5;i++){nparc += flagcor[i];}
    for(i =0; i<3;i++){nparnuis += flagnuis[i];}
    int_par = (int*)calloc((16), sizeof(int));
    dou_par = (double*)calloc((10), sizeof(double));
    
    //printf("NPARC:%d\tNPARNUIS:%d\n",nparc,nparnuis);
    
    int_par[0] = npts[0];
    int_par[1] = ntime[0];
    int_par[2] = flagcor[0];
    int_par[3] = flagcor[1];
    int_par[4] = flagcor[2];
    int_par[5] = flagcor[3];
    int_par[6] = flagcor[4];
    int_par[7] = flagnuis[0];
    int_par[8] = flagnuis[1];
    int_par[9] = flagnuis[2];
    int_par[10] = npar[0];//npar
    int_par[11] = weigthed[0];
    int_par[12] = dist[0];
    int_par[13] = nparc;
    int_par[14] = nparnuis;
    int_par[15] = cormod[0];
    
    dou_par[0] = maxtime[0];
    dou_par[1] = maxdist[0];
    dou_par[2] = parcor[0];
    dou_par[3] = parcor[1];//parcor
    dou_par[4] = parcor[2];//parcor
    dou_par[5] = parcor[3];//parcor
    dou_par[6] = parcor[4];//parcor
    dou_par[7] = nuis[0];
    dou_par[8] = nuis[1];
    dou_par[9] = nuis[2];
    
    
    
    cl_int err;
    cl_device_id        device;     // compute device id
    cl_context       context;       // compute context
    cl_command_queue commands;      // compute command queue
    cl_program       program;       // compute program
    cl_kernel        kernel;     // compute kernel
    
    
    // Vars for querying Device Info:
    
    char* value;
    size_t valueSize;
    cl_uint maxComputeUnits;
    cl_ulong long_entries;
    size_t p_size;
    
    // Set up OpenCL context, queue, kernel, etc.
    
    cl_uint deviceIndex = dev[0];
    
    // Get list of devices
    cl_device_id devices[MAX_DEVICES];
    unsigned numDevices = getDeviceList(devices);
    
    // Check device index in range
    if (deviceIndex >= numDevices)
    {
        //printf("Invalid device index (try '--list')\n");
        //return EXIT_FAILURE;
    }
    
    device = devices[deviceIndex];
    
    char name[MAX_INFO_STRING];
    getDeviceName(device, name);
    //printf("\nUsing OpenCL device: %s\n", name);
    
    
    // Get Device Info for Execution Model:
    
    // print hardware device version
    clGetDeviceInfo(device, CL_DEVICE_VERSION, 0, NULL, &valueSize);
    value = (char*) malloc(valueSize);
    clGetDeviceInfo(device, CL_DEVICE_VERSION, valueSize, value, NULL);
    //printf(" %d.%d Hardware version: %s\n", i+1, 1, value);
    free(value);
    
    // print software driver version
    clGetDeviceInfo(device, CL_DRIVER_VERSION, 0, NULL, &valueSize);
    value = (char*) malloc(valueSize);
    clGetDeviceInfo(device, CL_DRIVER_VERSION, valueSize, value, NULL);
    //printf(" %d.%d Software version: %s\n", i+1, 2, value);
    free(value);
    
    // print c version supported by compiler for device
    clGetDeviceInfo(device, CL_DEVICE_OPENCL_C_VERSION, 0, NULL, &valueSize);
    value = (char*) malloc(valueSize);
    clGetDeviceInfo(device, CL_DEVICE_OPENCL_C_VERSION, valueSize, value, NULL);
    //printf(" %d.%d OpenCL C version: %s\n", i+1, 3, value);
    free(value);
    
    // print parallel compute units
    clGetDeviceInfo(device, CL_DEVICE_MAX_COMPUTE_UNITS,
                    sizeof(maxComputeUnits), &maxComputeUnits, NULL);
    //printf(" %d.%d Parallel compute units: %d\n", i+1, 4, maxComputeUnits);
    
    
    clGetDeviceInfo(device,CL_DEVICE_GLOBAL_MEM_SIZE,sizeof(cl_ulong),&long_entries,NULL);
    //printf(" %d.%d Global Memory (MB):\t%llu\n",i+1, 5,long_entries/1024/1024);
    clGetDeviceInfo(device,CL_DEVICE_GLOBAL_MEM_CACHE_SIZE,sizeof(cl_ulong),&long_entries,NULL);
    //printf(" %d.%d Global Memory Cache (MB):\t%llu\n",i+1, 6,long_entries/1024/1024);
    clGetDeviceInfo(device,CL_DEVICE_LOCAL_MEM_SIZE,sizeof(cl_ulong),&long_entries,NULL);
    //printf(" %d.%d Local Memory (KB):\t%llu\n",i+1, 7,long_entries/1024);
    clGetDeviceInfo(device,CL_DEVICE_MAX_CLOCK_FREQUENCY,sizeof(cl_ulong),&long_entries,NULL);
    //printf(" %d.%d Max clock (MHz) :\t%llu\n",i+1, 8,long_entries);
    clGetDeviceInfo(device,CL_DEVICE_MAX_WORK_GROUP_SIZE,sizeof(size_t),&p_size,NULL);
    //printf(" %d.%d Max Work Group Size:\t%zu\n",i+1, 9,p_size);
    clGetDeviceInfo(device,CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS,sizeof(size_t),&p_size,NULL);
    //printf(" %d.%d Maximum dimensions:\t%zu\n",i+1, 10,p_size);
    clGetDeviceInfo(device,CL_DEVICE_MAX_MEM_ALLOC_SIZE,sizeof(size_t),&p_size,NULL);
    //printf(" %d.%d Max device buffer size (MB):\t%zu\n",i+1, 11,p_size/1024/1024);
    
    
    
    // Create a compute context
    context = clCreateContext(0, 1, &device, NULL, NULL, &err);
    checkError(err, "Creating context");
    
    // Create a command queue
    commands = clCreateCommandQueue(context, device, 0, &err);
    checkError(err, "Creating command queue");
    // Create the compute program from the source buffer
    //char *kernelsource = getKernelSource("Gneiting.cl");
    //program = clCreateProgramWithSource(context, 1, (const char **) & kernelsource, NULL, &err);
    //checkError(err, "Creating program");
    
    FILE *fp;
    size_t binary_size;
    char *binary_buf;
    cl_int binary_status;
    
    fp = fopen("Gneiting.cl.bin", "r");
    if (!fp) {
        Rprintf("Failed to load kernel.\n");
    }
    binary_buf = (char *)malloc(MAX_BINARY_SIZE);
    binary_size = fread(binary_buf, 1, MAX_BINARY_SIZE, fp);
    fclose(fp);
    
    program = clCreateProgramWithBinary(
                                        context, 1, &device, (const size_t *)&binary_size,
                                        (const unsigned char **)&binary_buf, &binary_status, &err
                                        );
    free(binary_buf);
    
    
    
    // Build the program
    err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
    if (err != CL_SUCCESS)
    {
        size_t len;
        char buffer[2048*4];
        
        Rprintf("Error: Failed to build program executable!\n%s\n", err_code(err));
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
        SEP;
        Rprintf("Build Log:\n%s\n", buffer);
        SEP;
        //return EXIT_FAILURE;
    }
    // Create the compute kernel from the program
    
    kernel = clCreateKernel(program, "scalarspaceocl", &err);
    checkError(err, "Creating kernel");
    
    // Find kernel work-group size
    size_t work_group_size;
    err = clGetKernelWorkGroupInfo (kernel, device, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &work_group_size, NULL);
    checkError(err, "Getting kernel work group info");
    //printf("Maximum Kernel work-group size: %lu\n", work_group_size);
    
    
    // Creating buffers
    size_t coords_buff = sizeof(double) * npts[0];
    size_t coordt_buff = sizeof(double) * ntime[0];
    size_t data_buff = sizeof(double) * npts[0]*ntime[0];
    size_t int_par_buff = sizeof(int) * (16);
    size_t dou_par_buff = sizeof(double) * (16);
    
    size_t g1,g2;
    const int ll1 =local_wi[0];
    const int ll2 =local_wi[1];
    g1 = npts[0] + (ll1 - (npts[0] & (ll1-1))); // SPACE
    g2 = ntime[0] + (ll2 - (ntime[0] & (ll2-1))); //TIME
    
    
    //printf("GLOBAL:\t%zu\t%zu\n",g1,g2);
    size_t local[2] = {ll1,ll2};
    size_t global[2] = { g1,g2};
    
    int length = g1*g2;
    //printf("LENGTH: %d\n",length);
    size_t length_buff = sizeof(double)* (length);
    
    double *h_mom_cond0,*h_mom_cond1,*h_mom_cond2,*h_mom_cond3;
    h_mom_cond0= (double*)calloc(length, sizeof(double));
    h_mom_cond1= (double*)calloc(length, sizeof(double));
    h_mom_cond2= (double*)calloc(length, sizeof(double));
    h_mom_cond3= (double*)calloc(length, sizeof(double));
    
    
    
    //st_alloc = wtime();
    cl_mem d_x = clCreateBuffer(context, CL_MEM_READ_ONLY, coords_buff, NULL, &err);
    checkError(err, "Creating buffer device X");
    err = clEnqueueWriteBuffer(commands, d_x, CL_TRUE, 0, coords_buff, (void*)h_x, 0, NULL, NULL);
    checkError(err, "Writing buffer device X");
    
    
    cl_mem d_y = clCreateBuffer(context, CL_MEM_READ_ONLY, coords_buff, NULL, &err);
    err = clEnqueueWriteBuffer(commands, d_y, CL_TRUE, 0, coords_buff, (void*)h_y, 0, NULL, NULL);
    checkError(err, "Creating buffer device Y");
    
    
    cl_mem d_data = clCreateBuffer(context, CL_MEM_READ_ONLY, data_buff, NULL, &err);
    err = clEnqueueWriteBuffer(commands, d_data, CL_TRUE, 0, data_buff, (void*)h_data, 0, NULL, NULL);
    checkError(err, "Creating buffer device data");
    
    cl_mem d_t = clCreateBuffer(context, CL_MEM_READ_ONLY, coordt_buff, NULL, &err);
    err = clEnqueueWriteBuffer(commands, d_t, CL_TRUE, 0, coordt_buff, (void*)h_t, 0, NULL, NULL);
    checkError(err, "Creating buffer device time");
    
    cl_mem m0_sol = clCreateBuffer(context, CL_MEM_WRITE_ONLY, length_buff, NULL, &err);
    checkError(err, "Creating buffer device m0_sol");
    err = clEnqueueWriteBuffer(commands, m0_sol, CL_TRUE, 0, length_buff, (void*)h_mom_cond0, 0, NULL, NULL);
    checkError(err, "Writing buffer device m0_sol");
    
    
    cl_mem m1_sol = clCreateBuffer(context, CL_MEM_WRITE_ONLY, length_buff, NULL, &err);
    checkError(err, "Creating buffer device m1_sol");
    err = clEnqueueWriteBuffer(commands, m0_sol, CL_TRUE, 0, length_buff, (void*)h_mom_cond1, 0, NULL, NULL);
    checkError(err, "Writing buffer device m1_sol");
    
    cl_mem m2_sol = clCreateBuffer(context, CL_MEM_WRITE_ONLY, length_buff, NULL, &err);
    checkError(err, "Creating buffer device m2_sol");
    err = clEnqueueWriteBuffer(commands, m0_sol, CL_TRUE, 0, length_buff, (void*)h_mom_cond2, 0, NULL, NULL);
    checkError(err, "Writing buffer device m2_sol");
    
    cl_mem m3_sol = clCreateBuffer(context, CL_MEM_WRITE_ONLY, length_buff, NULL, &err);
    checkError(err, "Creating buffer device m3_sol");
    err = clEnqueueWriteBuffer(commands, m0_sol, CL_TRUE, 0, length_buff, (void*)h_mom_cond3, 0, NULL, NULL);
    checkError(err, "Writing buffer device m3_sol");
    
    cl_mem d_int_par = clCreateBuffer(context, CL_MEM_READ_ONLY, int_par_buff, NULL, &err);
    err = clEnqueueWriteBuffer(commands, d_int_par, CL_TRUE, 0, int_par_buff, (void*)int_par, 0, NULL, NULL);
    checkError(err, "Creating buffer device int params");
    
    
    cl_mem d_dou_par = clCreateBuffer(context, CL_MEM_READ_ONLY, dou_par_buff, NULL, &err);
    err = clEnqueueWriteBuffer(commands, d_dou_par, CL_TRUE, 0, dou_par_buff, (void*)dou_par, 0, NULL, NULL);
    checkError(err, "Creating buffer device double params");
    
    
    // Push the data out to device
    clFinish(commands);
    
    err = clSetKernelArg(kernel,  0, sizeof(cl_mem), &d_t); //coordt
    err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &d_x); //coordx
    err |= clSetKernelArg(kernel, 2, sizeof(cl_mem), &d_y); //coordx
    err |= clSetKernelArg(kernel, 3, sizeof(cl_mem), &d_data); //data
    err |= clSetKernelArg(kernel, 4, sizeof(cl_mem), &m0_sol); //mom_cond0
    err |= clSetKernelArg(kernel, 5, sizeof(cl_mem), &m1_sol); //mom_cond1
    err |= clSetKernelArg(kernel, 6, sizeof(cl_mem), &m2_sol); //mom_cond2
    err |= clSetKernelArg(kernel, 7, sizeof(cl_mem), &m3_sol); //mom_cond3
    err |= clSetKernelArg(kernel, 8, sizeof(cl_mem), &d_int_par); //mom_cond3
    err |= clSetKernelArg(kernel, 9, sizeof(cl_mem), &d_dou_par); //mom_cond3
    checkError(err, "Setting kernel args length");
    
    
    // Execute the kernel
    
    err = clEnqueueNDRangeKernel(commands, kernel, 2, NULL, global,local, 0, NULL, NULL);
    checkError(err,"clEnqueueNDRangeKernel\n");
    clFinish(commands);
    
    
    err = clEnqueueReadBuffer(commands, m0_sol, CL_TRUE, 0,length_buff, h_mom_cond0, 0, NULL, NULL);
    err = clEnqueueReadBuffer(commands, m1_sol, CL_TRUE, 0,length_buff, h_mom_cond1, 0, NULL, NULL);
    err = clEnqueueReadBuffer(commands, m2_sol, CL_TRUE, 0,length_buff, h_mom_cond2, 0, NULL, NULL);
    err = clEnqueueReadBuffer(commands, m3_sol, CL_TRUE, 0,length_buff, h_mom_cond3, 0, NULL, NULL);
    clFinish(commands);
    
    
    
    double m0=0,m1=0,m2=0,m3=0;
    for(i=0; i<(npts[0]*ntime[0]); i++)
    {
        m0 += h_mom_cond0[i];
        m1 += h_mom_cond1[i];
        m2 += h_mom_cond2[i];
        m3 += h_mom_cond3[i];
    }
    mom_cond[0] =m0;mom_cond[1] =m1;mom_cond[2] =m2;mom_cond[3] =m3;
    //printf("final result: %.4f\t%.4f\t%.4f\t%.4f\n", m0,m1,m2,m3);
    
    
    // clean up inside kernels
    
    clReleaseProgram(program);
    clReleaseKernel(kernel);
    clReleaseCommandQueue(commands);
    clReleaseContext(context);
    
    
    clReleaseMemObject(d_x);
    clReleaseMemObject(d_y);
    clReleaseMemObject(d_data);
    clReleaseMemObject(d_t);
    clReleaseMemObject(m0_sol);
    clReleaseMemObject(m1_sol);
    clReleaseMemObject(m2_sol);
    clReleaseMemObject(m3_sol);
    clReleaseMemObject(d_int_par);
    clReleaseMemObject(d_dou_par);
    
    
    free(h_mom_cond0);
    free(h_mom_cond1);
    free(h_mom_cond2);
    free(h_mom_cond3);
    
    
}


void WendOCL(int *npts, int *ntime,double *coordt, double *maxtime,double *maxdist,int *cormod, double *parcor,int *flagcor,int *flagnuis,int *npar,double *nuis,double *data,int *weigthed, double *mom_cond, int *dist, double *coordx,double *coordy,double *gradcor,double  *grad, double *ww,int *local_wi, int *dev)
{
    
    
    
    // Setting data and params for sclar_space function
    double *h_x,*h_y,*h_data, *h_t;
    h_x = (double *) R_alloc(npts[0], sizeof(double));
    h_y = (double *) R_alloc(npts[0], sizeof(double));
    h_data = (double *) R_alloc(npts[0]*ntime[0], sizeof(double));
    h_t = (double *) R_alloc(ntime[0], sizeof(double));
    
    h_x = coordx;
    h_y = coordy;
    h_data = data;
    h_t = coordt;
    
    
    int i;
    int nparc=0, nparnuis=0;
    int *int_par;
    double *dou_par; // objects to pass int and double parametes to the kernel function
    //int_par = (int*)calloc(1+1+2+3+1+1+1+1+1, sizeof(int));//length sum of: npts+ntime+flagcor+flagnuis+npar+weigthed+dist+nparc+nparnuis+cormod
    //dou_par = (double*)calloc(1+1+2+3, sizeof(int));//length sum of: maxtime+maxdist+parcor+nuis
    
    
    
    for(i =0; i<7;i++){nparc += flagcor[i];}
    for(i =0; i<3;i++){nparnuis += flagnuis[i];} // pilas con el numero de parametros
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    
    //printf("NPARC:%d\tNPARNUIS:%d\n",nparc,nparnuis);
    
    int_par[0] = npts[0];
    int_par[1] = ntime[0];
    int_par[2] = flagcor[0];
    int_par[3] = flagcor[1];
    int_par[4] = flagcor[2];
    int_par[5] = flagcor[3];
    int_par[6] = flagcor[4];
    int_par[7] = flagnuis[0];
    int_par[8] = flagnuis[1];
    int_par[9] = flagnuis[2];
    int_par[10] = npar[0];//npar
    int_par[11] = weigthed[0];
    int_par[12] = dist[0];
    int_par[13] = nparc;
    int_par[14] = nparnuis;
    int_par[15] = cormod[0];
    
    int_par[16] = flagcor[5];
    int_par[17] = flagcor[6];
    
    dou_par[0] = maxtime[0];
    dou_par[1] = maxdist[0];
    dou_par[2] = parcor[0];
    dou_par[3] = parcor[1];//parcor
    dou_par[4] = parcor[2];//parcor
    dou_par[5] = parcor[3];//parcor
    dou_par[6] = parcor[4];//parcor
    dou_par[7] = nuis[0];
    dou_par[8] = nuis[1];
    dou_par[9] = nuis[2];
    
    dou_par[10] = parcor[5];//parcor
    dou_par[11] = parcor[6];//parcor
    
    double *h_grad,*h_gradcor;
    h_gradcor = (double*)calloc(nparc, sizeof(double));
    h_grad = (double*)calloc(npar[0], sizeof(double));
    
    //printf("CL wend: %f %f %f %f %f %f %f\n",parcor[0],parcor[1],parcor[2],parcor[3],parcor[4],dou_par[10],dou_par[11]);
    
    cl_int err;
    cl_device_id        device;     // compute device id
    cl_context       context;       // compute context
    cl_command_queue commands;      // compute command queue
    cl_program       program;       // compute program
    cl_kernel        kernel;     // compute kernel
    
    
    // Vars for querying Device Info:
    
    char* value;
    size_t valueSize;
    cl_uint maxComputeUnits;
    cl_ulong long_entries;
    size_t p_size;
    
    // Set up OpenCL context, queue, kernel, etc.
    
    cl_uint deviceIndex = dev[0];
    
    // Get list of devices
    cl_device_id devices[MAX_DEVICES];
    unsigned numDevices = getDeviceList(devices);
    
    // Check device index in range
    if (deviceIndex >= numDevices)
    {
        Rprintf("Invalid device index (try '--list')\n");
        //return EXIT_FAILURE;
    }
    
    device = devices[deviceIndex];
    
    char name[MAX_INFO_STRING];
    getDeviceName(device, name);
    //printf("\nUsing OpenCL device: %s\n", name);
    
    
    // Get Device Info for Execution Model:
    
    // print hardware device version
    clGetDeviceInfo(device, CL_DEVICE_VERSION, 0, NULL, &valueSize);
    value = (char*) malloc(valueSize);
    clGetDeviceInfo(device, CL_DEVICE_VERSION, valueSize, value, NULL);
    //printf(" %d.%d Hardware version: %s\n", i+1, 1, value);
    free(value);
    
    
    
    // print software driver version
    clGetDeviceInfo(device, CL_DRIVER_VERSION, 0, NULL, &valueSize);
    value = (char*) malloc(valueSize);
    clGetDeviceInfo(device, CL_DRIVER_VERSION, valueSize, value, NULL);
    //printf(" %d.%d Software version: %s\n", i+1, 2, value);
    free(value);
    
    
    
    // print c version supported by compiler for device
    clGetDeviceInfo(device, CL_DEVICE_OPENCL_C_VERSION, 0, NULL, &valueSize);
    value = (char*) malloc(valueSize);
    clGetDeviceInfo(device, CL_DEVICE_OPENCL_C_VERSION, valueSize, value, NULL);
    //printf(" %d.%d OpenCL C version: %s\n", i+1, 3, value);
    free(value);
    
    
    
    // print parallel compute units
    clGetDeviceInfo(device, CL_DEVICE_MAX_COMPUTE_UNITS,
                    sizeof(maxComputeUnits), &maxComputeUnits, NULL);
    //printf(" %d.%d Parallel compute units: %d\n", i+1, 4, maxComputeUnits);
    
    
    clGetDeviceInfo(device,CL_DEVICE_GLOBAL_MEM_SIZE,sizeof(cl_ulong),&long_entries,NULL);
    //printf(" %d.%d Global Memory (MB):\t%llu\n",i+1, 5,long_entries/1024/1024);
    clGetDeviceInfo(device,CL_DEVICE_GLOBAL_MEM_CACHE_SIZE,sizeof(cl_ulong),&long_entries,NULL);
    //printf(" %d.%d Global Memory Cache (MB):\t%llu\n",i+1, 6,long_entries/1024/1024);
    clGetDeviceInfo(device,CL_DEVICE_LOCAL_MEM_SIZE,sizeof(cl_ulong),&long_entries,NULL);
    //printf(" %d.%d Local Memory (KB):\t%llu\n",i+1, 7,long_entries/1024);
    clGetDeviceInfo(device,CL_DEVICE_MAX_CLOCK_FREQUENCY,sizeof(cl_ulong),&long_entries,NULL);
    //printf(" %d.%d Max clock (MHz) :\t%llu\n",i+1, 8,long_entries);
    clGetDeviceInfo(device,CL_DEVICE_MAX_WORK_GROUP_SIZE,sizeof(size_t),&p_size,NULL);
    //printf(" %d.%d Max Work Group Size:\t%zu\n",i+1, 9,p_size);
    clGetDeviceInfo(device,CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS,sizeof(size_t),&p_size,NULL);
    //printf(" %d.%d Maximum dimensions:\t%zu\n",i+1, 10,p_size);
    clGetDeviceInfo(device,CL_DEVICE_MAX_MEM_ALLOC_SIZE,sizeof(size_t),&p_size,NULL);
    //printf(" %d.%d Max device buffer size (MB):\t%zu\n",i+1, 11,p_size/1024/1024);
    
    
    
    // Create a compute context
    context = clCreateContext(0, 1, &device, NULL, NULL, &err);
    checkError(err, "Creating context");
    
    
    
    // Create a command queue
    commands = clCreateCommandQueue(context, device, 0, &err);
    checkError(err, "Creating command queue");
    // Create the compute program from the source buffer
    
    
    //char *kernelsource = getKernelSource("Wend.cl");
    //program = clCreateProgramWithSource(context, 1, (const char **) & kernelsource, NULL, &err);
    //checkError(err, "Creating program");
    
    FILE *fp;
    size_t binary_size;
    char *binary_buf;
    cl_int binary_status;
    
    fp = fopen("Wend.cl.bin", "r");
    if (!fp) {
        Rprintf("Failed to load kernel.\n");
    }
    binary_buf = (char *)malloc(MAX_BINARY_SIZE);
    binary_size = fread(binary_buf, 1, MAX_BINARY_SIZE, fp);
    fclose(fp);
    
    program = clCreateProgramWithBinary(
                                        context, 1, &device, (const size_t *)&binary_size,
                                        (const unsigned char **)&binary_buf, &binary_status, &err
                                        );
    free(binary_buf);
    
    
    
    // Build the program
    err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
    if (err != CL_SUCCESS)
    {
        size_t len;
        char buffer[2048*100];
        
        Rprintf("Error: Failed to build program executable!\n%s\n", err_code(err));
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
        SEP;
        Rprintf("Build Log:\n%s\n", buffer);
        SEP;
        //return EXIT_FAILURE;
    }
    // Create the compute kernel from the program
    
    kernel = clCreateKernel(program, "scalarspaceocl", &err);
    checkError(err, "Creating kernel");
    
    // Find kernel work-group size
    size_t work_group_size;
    err = clGetKernelWorkGroupInfo (kernel, device, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &work_group_size, NULL);
    checkError(err, "Getting kernel work group info");
    //printf("Maximum Kernel work-group size: %lu\n", work_group_size);
    
    
    // Creating buffers
    size_t coords_buff = sizeof(double) * npts[0];
    size_t coordt_buff = sizeof(double) * ntime[0];
    size_t data_buff = sizeof(double) * npts[0]*ntime[0];
    size_t int_par_buff = sizeof(int) * (int)(32);
    size_t dou_par_buff = sizeof(double) * (double)(32);
    
    size_t g1,g2;
    const int ll1 =local_wi[0];
    const int ll2 =local_wi[1];
    g1 = npts[0] + (ll1 - (npts[0] & (ll1-1))); // SPACE
    g2 = ntime[0] + (ll2 - (ntime[0] & (ll2-1))); //TIME
    
    
    //printf("GLOBAL:\t%zu\t%zu\n",g1,g2);
    size_t local[2] = {ll1,ll2};
    size_t global[2] = { g1,g2};
    
    int length = g1*g2;
    //printf("LENGTH: %d\n",length);
    size_t length_buff = sizeof(double)* (length);
    //size_t n_grad = sizeof(double)* (npar[0]);
    //size_t n_gradcor = sizeof(double)* (nparc);
    
    double *h_mom_cond0,*h_mom_cond1,*h_mom_cond2,*h_mom_cond3;
    h_mom_cond0= (double*)calloc(length, sizeof(double));
    h_mom_cond1= (double*)calloc(length, sizeof(double));
    h_mom_cond2= (double*)calloc(length, sizeof(double));
    h_mom_cond3= (double*)calloc(length, sizeof(double));
    
    
    
    //st_alloc = wtime();
    cl_mem d_x = clCreateBuffer(context, CL_MEM_READ_ONLY, coords_buff, NULL, &err);
    checkError(err, "Creating buffer device X");
    err = clEnqueueWriteBuffer(commands, d_x, CL_TRUE, 0, coords_buff, (void*)h_x, 0, NULL, NULL);
    checkError(err, "Writing buffer device X");
    
    
    cl_mem d_y = clCreateBuffer(context, CL_MEM_READ_ONLY, coords_buff, NULL, &err);
    err = clEnqueueWriteBuffer(commands, d_y, CL_TRUE, 0, coords_buff, (void*)h_y, 0, NULL, NULL);
    checkError(err, "Creating buffer device Y");
    
    
    cl_mem d_data = clCreateBuffer(context, CL_MEM_READ_ONLY, data_buff, NULL, &err);
    err = clEnqueueWriteBuffer(commands, d_data, CL_TRUE, 0, data_buff, (void*)h_data, 0, NULL, NULL);
    checkError(err, "Creating buffer device data");
    
    cl_mem d_t = clCreateBuffer(context, CL_MEM_READ_ONLY, coordt_buff, NULL, &err);
    err = clEnqueueWriteBuffer(commands, d_t, CL_TRUE, 0, coordt_buff, (void*)h_t, 0, NULL, NULL);
    checkError(err, "Creating buffer device time");
    
    cl_mem m0_sol = clCreateBuffer(context, CL_MEM_WRITE_ONLY, length_buff, NULL, &err);
    checkError(err, "Creating buffer device m0_sol");
    err = clEnqueueWriteBuffer(commands, m0_sol, CL_TRUE, 0, length_buff, (void*)h_mom_cond0, 0, NULL, NULL);
    checkError(err, "Writing buffer device m0_sol");
    
    
    cl_mem m1_sol = clCreateBuffer(context, CL_MEM_WRITE_ONLY, length_buff, NULL, &err);
    checkError(err, "Creating buffer device m1_sol");
    err = clEnqueueWriteBuffer(commands, m0_sol, CL_TRUE, 0, length_buff, (void*)h_mom_cond1, 0, NULL, NULL);
    checkError(err, "Writing buffer device m1_sol");
    
    cl_mem m2_sol = clCreateBuffer(context, CL_MEM_WRITE_ONLY, length_buff, NULL, &err);
    checkError(err, "Creating buffer device m2_sol");
    err = clEnqueueWriteBuffer(commands, m0_sol, CL_TRUE, 0, length_buff, (void*)h_mom_cond2, 0, NULL, NULL);
    checkError(err, "Writing buffer device m2_sol");
    
    cl_mem m3_sol = clCreateBuffer(context, CL_MEM_WRITE_ONLY, length_buff, NULL, &err);
    checkError(err, "Creating buffer device m3_sol");
    err = clEnqueueWriteBuffer(commands, m0_sol, CL_TRUE, 0, length_buff, (void*)h_mom_cond3, 0, NULL, NULL);
    checkError(err, "Writing buffer device m3_sol");
    
    /*cl_mem m1_grad = clCreateBuffer(context, CL_MEM_WRITE_ONLY, n_grad, NULL, &err);
    checkError(err, "Creating buffer device m1_grad");
    err = clEnqueueWriteBuffer(commands, m1_grad, CL_TRUE, 0, n_grad, (void*)h_grad, 0, NULL, NULL);
    checkError(err, "Writing buffer device m1_grad");
    
    
    cl_mem m1_gradcor = clCreateBuffer(context, CL_MEM_READ_WRITE, n_gradcor, NULL, &err);
    checkError(err, "Creating buffer device m1_gradcor");
    err = clEnqueueWriteBuffer(commands, m1_gradcor, CL_TRUE, 0, n_gradcor, (void*)h_gradcor, 0, NULL, NULL);
    checkError(err, "Writing buffer device m1_grad");*/
    
    cl_mem d_int_par = clCreateBuffer(context, CL_MEM_READ_WRITE, int_par_buff, NULL, &err);
    err = clEnqueueWriteBuffer(commands, d_int_par, CL_TRUE, 0, int_par_buff, (void*)int_par, 0, NULL, NULL);
    checkError(err, "Creating buffer device int params");
    
    
    cl_mem d_dou_par = clCreateBuffer(context, CL_MEM_READ_ONLY, dou_par_buff, NULL, &err);
    err = clEnqueueWriteBuffer(commands, d_dou_par, CL_TRUE, 0, dou_par_buff, (void*)dou_par, 0, NULL, NULL);
    checkError(err, "Creating buffer device double params");
    
    
    // Push the data out to device
    clFinish(commands);
    
    err = clSetKernelArg(kernel,  0, sizeof(cl_mem), &d_t); //coordt
    err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &d_x); //coordx
    err |= clSetKernelArg(kernel, 2, sizeof(cl_mem), &d_y); //coordx
    err |= clSetKernelArg(kernel, 3, sizeof(cl_mem), &d_data); //data
    err |= clSetKernelArg(kernel, 4, sizeof(cl_mem), &m0_sol); //mom_cond0
    err |= clSetKernelArg(kernel, 5, sizeof(cl_mem), &m1_sol); //mom_cond1
    err |= clSetKernelArg(kernel, 6, sizeof(cl_mem), &m2_sol); //mom_cond2
    err |= clSetKernelArg(kernel, 7, sizeof(cl_mem), &m3_sol); //mom_cond3
    err |= clSetKernelArg(kernel, 8, sizeof(cl_mem), &d_int_par); //mom_cond3
    err |= clSetKernelArg(kernel, 9, sizeof(cl_mem), &d_dou_par); //mom_cond3
   // err |= clSetKernelArg(kernel, 10, sizeof(cl_mem), &m1_grad); //grad
   // err |= clSetKernelArg(kernel, 11, sizeof(cl_mem), &m1_gradcor); //gradcor
    checkError(err, "Setting kernel args length");
    
    
    // Execute the kernel
    
    err = clEnqueueNDRangeKernel(commands, kernel, 2, NULL, global,local, 0, NULL, NULL);
    checkError(err,"clEnqueueNDRangeKernel\n");
    clFinish(commands);
    
    
    err = clEnqueueReadBuffer(commands, m0_sol, CL_TRUE, 0,length_buff, h_mom_cond0, 0, NULL, NULL);
    err = clEnqueueReadBuffer(commands, m1_sol, CL_TRUE, 0,length_buff, h_mom_cond1, 0, NULL, NULL);
    err = clEnqueueReadBuffer(commands, m2_sol, CL_TRUE, 0,length_buff, h_mom_cond2, 0, NULL, NULL);
    err = clEnqueueReadBuffer(commands, m3_sol, CL_TRUE, 0,length_buff, h_mom_cond3, 0, NULL, NULL);
    clFinish(commands);
    
    
    
    double m0=0,m1=0,m2=0,m3=0;
    for(i=0; i<(npts[0]*ntime[0]); i++)
    {
        m0 += h_mom_cond0[i];
        m1 += h_mom_cond1[i];
        m2 += h_mom_cond2[i];
        //printf("h_mom_cond1[i] %f h_mom_cond2[i] %f\n",h_mom_cond1[i],h_mom_cond1[i]);
        m3 += h_mom_cond3[i];
    }
    mom_cond[0] =m0;mom_cond[1] =m1;mom_cond[2] =m2;mom_cond[3] =m3;
    //printf("CL: final result: %.4f\t%.4f\t%.4f\t%.4f\n", m0,m1,m2,m3);
    
    
    // clean up inside kernels
    
    clReleaseProgram(program);
    clReleaseKernel(kernel);
    clReleaseCommandQueue(commands);
    clReleaseContext(context);
    
    
    clReleaseMemObject(d_x);
    clReleaseMemObject(d_y);
    clReleaseMemObject(d_data);
    clReleaseMemObject(d_t);
    clReleaseMemObject(m0_sol);
    clReleaseMemObject(m1_sol);
    clReleaseMemObject(m2_sol);
    clReleaseMemObject(m3_sol);
    clReleaseMemObject(d_int_par);
    clReleaseMemObject(d_dou_par);
    //clReleaseMemObject(m1_grad);
    //clReleaseMemObject(m1_gradcor);
    
    
    free(h_mom_cond0);
    free(h_mom_cond1);
    free(h_mom_cond2);
    free(h_mom_cond3);
    
    free(int_par);
    free(dou_par);
    
    free(h_x);
    free(h_y);
    free(h_data);
    free(h_t);
    
    free(h_gradcor);
    free(h_grad);
    
    
}







