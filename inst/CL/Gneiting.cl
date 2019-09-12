/********************************************************************/
/************ gradient of the pairwise CL ****************/
/********************************************************************/
double *Grad_Pair_Gauss(double rho, int flag0, int flag1, int flag2, double *gradcor, double *grad,
                     int npar, double par0,double par1,double par2, double u, double v)
{
    // Initialization variables:
    double gcor[2];
    gcor[0] =gradcor[0];gcor[1] =gradcor[1];
    double mean=par0,nugget=par1,sill=par2;
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
    if(flag1==1){grad[i]=0.5*k*(R*d-L*4*b*a-2*a*(pa-pb));i++;}
    
    // Derivative respect with the sill
    if(flag2==1)
    {
        grad[i]=-0.5*k*(2*(pa*a-pb*(2*sill+3*nugget)+rho*b*(pb-pn))+R*(c+2*nugget*b*rho)+2*L*rho*(ps-pn-pb));
        i++;
    }
    
    // Derivatives with respect to the correlation parameters
    h=0;
    C=-k*sill*(R*a*b-L*d+b*c);
    
    for(j=i;j<npar;j++) // CON I+1 LA PARTE A SE ARREGLA
    {
        grad[j]=(gcor[h])*C;
        h++;
        
    }
    return grad;
}
// Derivatives with respect to the separable parameter of the Gneiting correlation model:
double DGneiting_sep(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep)
{
    double a=0,arg=0,rho=0;
    arg=1+pow(u/scale_t, power_t);
    rho=exp(-pow(h/scale_s, power_s)/(pow(arg, 0.5*sep*power_s)))/arg;
    if(arg) a=0.5*pow(h/scale_s, power_s)*pow(arg,-0.5*sep*power_s)*power_s*rho*log(arg);
    return(a);
}
// Derivatives with respect to spatial power of the Gneiting correlation model:
double DGneiting_pw_s(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep)
{
    double arg=0.0,rho=0.0,a=0.0;
    arg=1+pow(u/scale_t, power_t);
    rho=exp(-pow(h/scale_s, power_s)/(pow(arg, 0.5*sep*power_s)))/arg;
    if(h && arg)
    {
        a=(-pow(h/scale_s, power_s)*pow(arg,-0.5*sep*power_s)*log(h/scale_s) + 0.5*pow(h/scale_s, power_s)*pow(arg,-0.5*sep*power_s)*sep*log(arg) )*rho;
    }
    return(a);
}

// Derivatives with respect to temporal power of the Gneiting correlation model:
double DGneiting_pw_t(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep)
{
    double arg=0.0,rho=0.0,a=0.0;
    arg=1+pow(u/scale_t, power_t);
    rho=exp(-pow(h/scale_s, power_s)/(pow(arg, 0.5*sep*power_s)))/arg;
    if(u)
    {
        a=( -pow(u/scale_t, power_t)*log(u/scale_t)*rho+0.5*rho*pow(h/scale_s, power_s)*power_s*sep*log(u/scale_t)*pow(arg,-0.5*sep*power_s))/arg;
    }
    return(a) ;
}

// Derivatives with respect to the temporal scale parameter of the Gneiting correlation model:
double DGneiting_sc_t(double h,double u, double power_s,double power_t,
                      double scale_s,double scale_t,double sep)
{
    double arg=0.0,rho=0.0,a=0.0;
    arg=1+pow(u/scale_t, power_t);
    rho=exp(-pow(h/scale_s, power_s)/(pow(arg, 0.5*sep*power_s)))/arg;
    if (arg) a= (power_t*pow(u/scale_t, power_t)*rho)/(scale_t*pow(arg,2))-( 0.5*power_t*power_s*sep*rho*pow(h/scale_s, power_s)*pow(u/scale_t, power_t)*pow(arg,-0.5*sep*power_s-2))/scale_t;
    return(a);
}
// Derivatives with respect to the spatial scale parameter of the Gneiting correlation model:
double DGneiting_sc_s(double h,double u,double power_s,double power_t,
                      double scale_s,double scale_t,double sep)
{
    double arg=0.0,rho=0.0,a=0.0;
    arg=1+pow(u/scale_t, power_t);
    rho=exp(-pow(h/scale_s, power_s)/(pow(arg, 0.5*sep*power_s)))/arg;
    a=(pow(h/scale_s, power_s)*power_s*rho*pow(arg,-0.5*power_s*sep-1))/scale_s;
    
    return (a);
}

//  Derivatives with respect ot the correlations parameters:
void GradCorrFct(double rho, int flag0, int flag1, int flag2, int flag3, int flag4,
                 double *grad, double h, double u,double par0,double par1,double par2,double par3,double par4)
{
    int i=0;
    double power_s=0, power_t=0;
    double scale_s=0, scale_t=0, sep=0;
    //spatial gradients of correlations:
    //Gneiting spatio-temporal correlation
    power_s=par0;
    power_t=par1;
    scale_s=par2;
    scale_t=par3;
    sep=par4;
    
    if(flag0==1){//spatial-power parameter
        grad[i]=DGneiting_pw_s(h,u,power_s,power_t,scale_s,scale_t,sep);i++;}
    if(flag1==1){//temporal-power parameter
        grad[i]=DGneiting_pw_t(h,u,power_s,power_t,scale_s,scale_t,sep);i++;}
    if(flag2==1){//spatial-scale parameter
        grad[i]=DGneiting_sc_s(h,u,power_s,power_t,scale_s,scale_t,sep);i++;}
    if(flag3==1){//temporal-scale parameter
        grad[i]=DGneiting_sc_t(h,u,power_s,power_t,scale_s,scale_t,sep);i++;}
    //separable parameter
    if(flag4==1)
        grad[i]=DGneiting_sep(h,u,power_s,power_t,scale_s,scale_t,sep);
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
double CorFct(double h, double u, double par0,double par1,double par2,double par3,double par4)
{
    double arg=0.0,  power_s=0.0, power_t=0.0;
    double rho=0.0, sep=0,  scale_s=0.0, scale_t=0;
    
    power_s =   par0;
    power_t =   par1;
    scale_s =   par2;
    scale_t =   par3;
    sep     =   par4;
    arg=1+pow(u/scale_t, power_t);
    rho=exp(-(pow(h/scale_s, power_s))*pow(arg, -0.5*sep*power_s))/pow(arg,1);
    
    return rho;
}

__kernel void scalarspaceocl(__global const double *coordt,__global const double *scoordx,__global const double *scoordy, __global const double *sdata, __global double *mom_cond0,__global double *mom_cond1,__global double *mom_cond2,__global double *mom_cond3,__global const int *int_par,__global const double *dou_par)
{
    
    int npts	=   int_par[0];
    int ntime	=   int_par[1];
    int flagcor0	=   int_par[2];
    int flagcor1	=   int_par[3];
    int flagcor2	=   int_par[4];
    int flagcor3	=   int_par[5];
    int flagcor4	=   int_par[6];
    int flagnuis0	=   int_par[7];
    int flagnuis1   =   int_par[8];
    int flagnuis2	=   int_par[9];
    int npar        =   int_par[10];
    int weigthed	=   int_par[11];
    int dist        =   int_par[12];
    int nparc       =   int_par[13];
    int nparnuis	=   int_par[14];
    int cormod      =   int_par[15];
    
    double maxtime	=   dou_par[0];
    double maxdist	=   dou_par[1];
    double parcor0	=   dou_par[2];
    double parcor1	=   dou_par[3];
    double parcor2	=   dou_par[4];
    double parcor3	=   dou_par[5];
    double parcor4	=   dou_par[6];
    double nuis0	=   dou_par[7];
    double nuis1	=   dou_par[8];
    double nuis2	=   dou_par[9];
    
    double grad[4];
    double gradcor[2];
    
    int l = get_global_id(0);
    int t = get_global_id(1);
    
    int ls = get_global_size(0);
    int ms = get_global_size(1);
    
    int m=0,v =0;
    double lagt=0.0f,lags=0.0, rho=0.0,weights=0.0;
    double ww0 =0,ww1 =0,ww2 =0,ww3 =0;
    
    int gid = (npts*t+l);
    
    int m1 = get_local_id(0);
    int v1 = get_local_id(1);
    
    int lsize_m = get_local_size(0);
    int lsize_v = get_local_size(1);
    
    int wx = (l-m1)/lsize_m;
    int wy = (t-v1)/lsize_v;
    
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
        for(m = l;m<npts;m++)
        {
            if(l==m)
            {
                for(v = t+1;v<ntime;v++)
                {
                    lagt=fabs(coordt[t]-coordt[v]);
                    
                    if(lagt<=maxtime)
                    {
                        //Computing correlation
                        
                        rho=CorFct(0,lagt,parcor0,parcor1,parcor2,parcor3,parcor4);
                        //Computing the gradient of the corr parameters
                        
                        GradCorrFct(rho,flagcor0,flagcor1,flagcor2,flagcor3,flagcor4,gradcor,0,lagt,parcor0,parcor1,parcor2,parcor3,parcor4);
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
            else
            {
                if(dist==1) lags=hypot(scoordx[l]-scoordx[m],scoordy[l]-scoordy[m]);
                for(v=0;v<ntime;v++)
                {
                    lagt=fabs(coordt[t]-coordt[v]);
                    if(lagt<=maxtime && lags<=maxdist)
                    {
                        rho=CorFct(lags,lagt,parcor0,parcor1,parcor2,parcor3,parcor4);
                        //Computing the gradient of the corr parameters
                        GradCorrFct(rho,flagcor0,flagcor1,flagcor2,flagcor3,flagcor4,gradcor,lags,lagt,parcor0,parcor1,parcor2,parcor3,parcor4);
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
