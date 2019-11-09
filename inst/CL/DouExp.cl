void Grad_Pair_Gauss(double rho, int flag0, int flag1, int flag2, double *gradcor,double *grad,int npar, double par0,double par1,double par2, double u, double v);
double DStabSc(double lag, double power, double scale, double rho);
void GradCorrFct(double rho, int flag0, int flag1,
                 double *grad, double h, double u,double par0,double par1);
double CorFunBohman(double lag,double scale);
double CorFunStable(double lag, double power, double scale);
double CorFct(double h, double u, double par0,double par1);



void Grad_Pair_Gauss(double rho, int flag0, int flag1, int flag2,double *gradcor, double *grad,
                     int npar, double par0,double par1,double par2, double u, double v)
{

    //printf("CL grad: %f\t%f\t%f\t%f\n\n",grad[0],grad[1],grad[2],grad[3]);
    //printf("CL gracor: %f\t%f\t%f\t%f\n",gradcor[0],gradcor[1],gradcor[2],gradcor[3]);
    // Initialization variables:
    double mean=par0,nugget=par1,sill=par2;
    //printf("CL %f %f %f\n",mean,nugget,sill);
    double a=nugget+sill,b=sill*rho,pa=a*a,pb=b*b;
    double c=-pa+pb,d=pa+pb,k=1/(c*c);
    double C=0.0,L=0.0,R=0.0;
    double pn=nugget*nugget,ps=sill*sill,pu=0.0, pv=0.0;
    int h=0, i=0, j=0;
    //defines useful quantities:
    u=u-mean;
    v=v-mean;
    pu=pow(u,2);
    pv=pow(v,2);
    R=pu+pv;L=u*v;
    // Derivatives  respect with the mean
    if(flag0==1){grad[i]=(u+v)/(a+b);i++;}
    // Derivative  respect with the nugget
    if(flag1==1)
    {
        grad[i]=0.5*k*(R*d-L*4*b*a-2*a*(pa-pb));i++;
    }
    // Derivative respect with the sill
    if(flag2==1)
    {
        grad[i]=-0.5*k*(2*(pa*a-pb*(2*sill+3*nugget)+rho*b*(pb-pn))+R*(c+2*nugget*b*rho)+2*L*rho*(ps-pn-pb));
        i++;
    }

    // Derivatives with respect to the correlation parameters
    h=0;
    C=-k*sill*(R*a*b-L*d+b*c);
    //printf("%f\n",C);
    //printf("%d %d\n",i , npar);
    for(j=i;j<npar;j++){grad[j]=C*gradcor[h];h++;}
    return;
}

// Derivatives with respect to scale of the Stable correlation model:
double DStabSc(double lag, double power, double scale, double rho)
{
    if(lag) return rho*power*pow(lag/scale,power)/scale;
    else return 0.0;
}

//  Derivatives with respect ot the correlations parameters:
void GradCorrFct(double rho, int flag0, int flag1,
                 double *grad, double h, double u,double par0,double par1)
{
    int i=0;
    double power_s=0, power_t=0;
    double scale_s=0, scale_t=0, sep=0;
    //spatial gradients of correlations:
    //Double Exponentil
    scale_s=par0;
    scale_t=par1;
    if(flag0==1)
    {//spatial-scale parameter
        grad[i]=DStabSc(h,1,scale_s,rho);
        i++;
    }
    //temporal-scale parameter
    if(flag1==1) grad[i]=DStabSc(u,1,scale_t,rho);
    return;
}
/********************************************************************/
/************** correlation  models **********************************/
/********************************************************************/
double CorFunBohman(double lag,double scale)
{
    double rho=0.0,x=0;
    x=lag/scale;
    if(x<=1)
    {
        if (x>0) rho=(1-x)*(sin(2*M_PI*x)/(2*M_PI*x))+(1-cos(2*M_PI*x))/(2*M_PI*M_PI*x);
        else   rho=1;
    }
    else rho=0;
    return rho;
}
// Stable class of correlation models:
double CorFunStable(double lag, double power, double scale)
{
    double rho=0.0;
    // Computes the correlation:
    rho=exp(-pow(lag/scale,power));
    return rho;
}
double CorFct(double h, double u, double par0,double par1)
{
    double arg=0.0,  power_s=0.0, power_t=0.0;
    double rho=0.0, sep=0,  scale_s=0.0, scale_t=0;
    
    scale_s=par0;
    scale_t=par1;
    rho=CorFunStable(h,1,scale_s)*CorFunStable(u,1,scale_t);
    return rho;
}


//*********************************************************************
__kernel void scalarspaceocl(__global const double *coordt,__global const double *scoordx,__global const double *scoordy, __global const double *sdata, __global double *mom_cond0,__global double *mom_cond1,__global double *mom_cond2,__global double *mom_cond3,__global const int *int_par,__global const double *dou_par)
{
    
    int npts	=	int_par[0];
    int ntime	=	int_par[1];
    int flagcor0	=	int_par[2];
    int flagcor1	=	int_par[3];
    int flagnuis0	=	int_par[7];
    int flagnuis1	=	int_par[8];
    int flagnuis2	=	int_par[9];
    int npar	=	int_par[10];
    int weigthed	=	int_par[11];
    int dist	=	int_par[12];
    int nparc	=	int_par[13];
    int nparnuis	=	int_par[14];
    int cormod	=	int_par[15];
    
    double maxtime	=	dou_par[0];
    double maxdist	=	dou_par[1];
    double parcor0	=	dou_par[2];
    double parcor1	=	dou_par[3];
    double nuis0	=	dou_par[7];
    double nuis1	=	dou_par[8];
    double nuis2	=	dou_par[9];
    
    
    double grad[4];
    //grad[0]=0;grad[1]=0;grad[2]=0;grad[3]=0;
    double gradcor[2];
    //printf("%d\n",nparc);
    //gradcor[0]=0;gradcor[1]=0;
    
    //int l=0,t=0,m=0,v=0,gid = get_global_id(0);
    
    int l = get_global_id(0);
    int t = get_global_id(1);
    
    int ls = get_global_size(0);
    int ms = get_global_size(1);
    
    int m=0,v =0;
    double lagt=0.0f,lags=0.0, rho=0.0,weights=0.0;
    double ww0 =0,ww1 =0,ww2 =0,ww3 =0;
    
    int gid = (npts*t+l);
    //int gid = (t+ntime*l);
    
    int m1 = get_local_id(0);
    int v1 = get_local_id(1);
    
    int lsize_m = get_local_size(0);
    int lsize_v = get_local_size(1);
    
    int wx = (l-m1)/lsize_m;
    int wy = (t-v1)/lsize_v;
    
    //int gidx = (npts*t+l);
    //int gidy = (ntime*l+t);
    
    int gidx = (wx*lsize_m+m1);
    int gidy = (wy*lsize_v+v1);
    
    int i = (npts*gidy+gidx);
    int j = (ntime*gidx+gidy);
    
    
    
    bool isValid = true;
    
    if(l >= npts) isValid = false;
    
    if(t >= ntime) isValid = false;
    
    if(isValid)
        
    {
        mom_cond0[i] = 0;mom_cond1[i] = 0;mom_cond2[i] = 0;mom_cond3[i] = 0;
        //printf("%d\t%d\t%d\n",gidx,l,t);
        for(m = l;m<npts;m++)
        {
            if(l==m)
            {
                for(v = t+1;v<ntime;v++)
                {
                    lagt=fabs(coordt[t]-coordt[v]);
                    
                    if(lagt<=maxtime)
                    {
                        if(!isnan(sdata[(t+ntime*l)])&&!isnan(sdata[(v+ntime*l)]) ){
                        //Computing correlation
                        
                        rho=CorFct(0,lagt,parcor0,parcor1);
                        //Computing the gradient of the corr parameters
                        
                        GradCorrFct(rho,flagcor0,flagcor1,gradcor,0,lagt,parcor0,parcor1);
                        //Compute the gradient of the composite likelihood:
                        
                    Grad_Pair_Gauss(rho,flagnuis0,flagnuis1,flagnuis2,gradcor,grad,npar,nuis0,nuis1,nuis2,sdata[(t+ntime*l)],sdata[(v+ntime*l)]);
                
                        if(weigthed)
                        {
                            weights=CorFunBohman(lagt,maxtime);
                            ww0 =1,ww1 =1,ww2 =weights,ww3 =weights;
                        }
                        else
                        {
                            ww0 =1,ww1 =1,ww2 =1,ww3 =1;
                        }
                        
                        mom_cond0[i]+=ww0*grad[0];
                        mom_cond1[i]+=ww1*grad[1];
                        mom_cond2[i]+=ww2*grad[2];
                        mom_cond3[i]+=ww3*grad[3];
                    }
                    }
                }
            }
            else
            {
                if(dist==1) lags=hypot(scoordx[l]-scoordx[m],scoordy[l]-scoordy[m]);
                for(v=0;v<ntime;v++)
                {
                    lagt=fabs(coordt[t]-coordt[v]);
                    if(lagt<=maxtime && lags<=maxdist)
                    {
                        if(!isnan(sdata[(t+ntime*l)])&&!isnan(sdata[(v+ntime*m)]) ){
                        rho=CorFct(lags,lagt,parcor0,parcor1);
                        //Computing the gradient of the corr parameters
                        GradCorrFct(rho,flagcor0,flagcor1,gradcor,lags,lagt,parcor0,parcor1);
                        //Compute the gradient of the composite likelihood:
                        
                        Grad_Pair_Gauss(rho,flagnuis0,flagnuis1,flagnuis2,gradcor,grad,npar,nuis0,nuis1,nuis2,sdata[(t+ntime*l)],sdata[(v+ntime*m)]);

                        if(weigthed)
                        {
                            weights=CorFunBohman(lags,maxdist)*CorFunBohman(lagt,maxtime);
                            ww0 =1,ww1 =1,ww2 =weights,ww3 =weights;
                        }
                        else
                        {
                            ww0 =1,ww1 =1,ww2 =1,ww3 =1;
                        }
                        mom_cond0[i]+=ww0*grad[0];
                        mom_cond1[i]+=ww1*grad[1];
                        mom_cond2[i]+=ww2*grad[2];
                        mom_cond3[i]+=ww3*grad[3];
                    }

                    }
                }
            }
        }
    }
}
