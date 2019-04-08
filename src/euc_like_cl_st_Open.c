#include "header.h"

#define EPS DBL_EPSILON
//#define EPS 1.0e-60

//---------START WENDLAND FUNCTIONS-----------

/* START Wendland covariance */

/* integrand  in  generalized wendland function*/
double int_gen(double x,double mu, double alpha,double lag,double supp)
{
    double res=0.0,y;
    y=lag/supp;
    res=R_pow(1-x,mu-1)*R_pow(x*x-y*y,alpha)/beta(2*alpha+1,mu);
    return (res);///(R_pow(2,alpha-1)*gamma(alpha)*R_pow(supp,2*alpha)));
}
void integr_gen(double *x, int n, void *ex){
    int i;double mu,alpha,beta,y;
    mu =    ((double*)ex)[0];  //mu
    alpha = ((double*)ex)[1];  //alpha
    beta =     ((double*)ex)[2];  //csupp
    y =     ((double*)ex)[3];  //h
    for (i=0;i<n;i++) {x[i]=int_gen(x[i],mu,alpha,y,beta);}
    return;
}
// function computing generalized wendland
double wendintegral(double x, double *param) {
    
    double ex[4], lower, upper, epsabs, epsrel, result, abserr, *work;
    int neval, ier, subdiv, lenw, last, *iwork;
    subdiv = 100;
    epsabs = R_pow(DOUBLE_EPS, 0.25);
    epsrel = epsabs;
    lenw = 4 * subdiv;           /* as instructed in WRE */
    iwork =   (int *) Calloc(subdiv, int);  /* idem */
    work = (double *) Calloc(lenw, double); /* idem */
    ex[0] = param[0]; ex[1] = param[1]; ex[2] = param[2];ex[3]=x;
    lower=x/param[2];
    upper=1;
    // Compute the integral
    
    if(x<=param[2]) {
        Rdqags(integr_gen, (void *) &ex,
               &lower, &upper, &epsabs, &epsrel, &result,
               &abserr, &neval, &ier, &subdiv, &lenw, &last, iwork, work);
        
    }else   {result=0;}
    Free(iwork);Free(work);
    return(result);
}
/* generalized wendland function*/
double CorFunW_gen(double lag,double R_power1,double smooth,double scale)  // mu alpha beta
{
    
    //Rprintf("PARAM A: %f %f %f %f\n",lag,R_power1,smooth,scale);
    double rho=0.0,x=0;
    if(smooth==0) {
        
        x=lag/scale;
        if(x<=1) rho=R_pow(1-x,R_power1);
        else rho=0;
        return(rho);
    }
    if(smooth==1) {
        
        x=lag/scale;
        if(x<=1) rho=R_pow(1-x,R_power1+1)*(1+x*(R_power1+1));
        else rho=0;
        return(rho);
    }
    if(smooth==2) {
        
        x=lag/scale;
        if(x<=1) rho=R_pow(1-x,R_power1+2)*(1+x*(R_power1+2)+x*x*(R_power1*R_power1 +4*R_power1 +3 )/3  );
        else rho=0;
        return(rho);
    }
    
    x=lag;
    double *param;
    param=(double *) Calloc(3,double);
    param[0]=R_power1;param[1]=smooth;param[2]=scale;  //mu,alpha //beta
    
    rho=wendintegral(x,param);
    Free(param);
    return(rho);
}


double wen_time(double *par, double h,double u)
{
    //Rprintf("%f %f %f %f %f %f\n",par[0],par[1],par[2],par[3],par[4],par[5]);
   // Rprintf("smooth: %f\n",par[6]);
    //0:power2_s 1:power_s  2:power2_t   3:scale_s     4:scale_t   5:sep 6:smooth
    double R_power_s=2;
    double R_power=par[0];
    double R_power_t=par[2];
    double scale_s=par[3];
    double scale_t=par[4];
    double sep=par[5];
    double smooth=par[6];
    //Rprintf("R_power_s: %f R_power: %f R_power_t: %f scale_s: %f scale_t: %f sep: %f smooth: %f \n",R_power_s,R_power,R_power_t,scale_s,scale_t,sep,smooth);
    double arg=R_pow(1+R_pow(h/scale_s,R_power_s/2),-1/(R_power_s/2));
    double rho=R_pow(arg,R_power)*CorFunW_gen(u,R_power_t,smooth,scale_t*R_pow(arg,sep));
    return(rho);
}



/* END Wendland covariance */


/* START DERIVATIVES Wendland covariance */

// SCALE_S:
double deri_scale_s_wen_time(double *par, double h,double u)
{
    double R_power=par[0];
    double R_power_t=par[2];
    double scale_s=par[3];
    double scale_t=par[4];
    double sep=par[5];
    double smooth=par[6];
    
    //double EPS = 1.0e-10;
    
    double delta=sqrt(EPS)*scale_s;
    
    double *par1;
    par1 = (double *) Calloc(7, double); /* idem */
    // define vector par1
    par1[0]=R_power;
    par1[1]=2;
    par1[2]=R_power_t;
    par1[3]=scale_s + delta;
    par1[4]=scale_t;
    par1[5]=sep;
    par1[6]=smooth;
    double grad=(wen_time(par1,h,u)-wen_time(par,h,u))/delta;
    //Rprintf("%f %f %f %f\n",wen_time(par1,h,u),wen_time(par,h,u),delta,EPS);
    return(grad);
}



// SCALE_T:
double deri_scale_t_wen_time(double *par, double h,double u)
{
    double R_power=par[0];
    double R_power_t=par[2];
    double scale_s=par[3];
    double scale_t=par[4];
    double sep=par[5];
    double smooth=par[6];
    
    //double EPS = 1.0e-10;
    
    double delta=sqrt(EPS)*scale_t;
    
    double *par1;
    par1 = (double *) Calloc(7, double); /* idem */
    // define vector par1
    par1[0]=R_power;
    par1[1]=2;
    par1[2]=R_power_t;
    par1[3]=scale_s;
    par1[4]=scale_t+ delta;
    par1[5]=sep;
    par1[6]=smooth;
    double grad=(wen_time(par1,h,u)-wen_time(par,h,u))/delta;
    //Rprintf("%f %f %f %f\n",wen_time(par1,h,u),wen_time(par,h,u),delta,EPS);
    return(grad);
}



// SMOOTH:
double deri_smooth_wen_time(double *par, double h,double u)
{
    double R_power=par[0];
    double R_power_t=par[2];
    double scale_s=par[3];
    double scale_t=par[4];
    double sep=par[5];
    double smooth=par[6];
    
    //double EPS = 1.0e-10;
    
    double delta=sqrt(EPS)*smooth;
    
    double *par1;
    par1 = (double *) Calloc(7, double); /* idem */
    // define vector par1
    par1[0]=R_power;
    par1[1]=2;
    par1[2]=R_power_t;
    par1[3]=scale_s;
    par1[4]=scale_t;
    par1[5]=sep;
    par1[6]=smooth+ delta;
    double grad=(wen_time(par1,h,u)-wen_time(par,h,u))/delta;
    //Rprintf("%f %f %f %f\n",wen_time(par1,h,u),wen_time(par,h,u),delta,EPS);
    return(grad);
}



// SILL:
double deri_sill_wen_time(double *par, double h,double u)
{
    double grad=wen_time(par,h,u);
    return(grad);
}



// SEP:
double deri_sep_wen_time(double *par, double h,double u)
{
    double R_power=par[0];
    double R_power_t=par[2];
    double scale_s=par[3];
    double scale_t=par[4];
    double sep=par[5];
    double smooth=par[6];
    
    //double EPS = 1.0e-10;
    
    double delta=sqrt(EPS)*sep;
    
    double *par1;
    par1 = (double *) Calloc(7, double); /* idem */
    // define vector par1
    par1[0]=R_power;
    par1[1]=2;
    par1[2]=R_power_t;
    par1[3]=scale_s;
    par1[4]=scale_t;
    par1[5]=sep+ delta;
    par1[6]=smooth;
    double grad=(wen_time(par1,h,u)-wen_time(par,h,u))/delta;
    //Rprintf("%f %f %f %f\n",wen_time(par1,h,u),wen_time(par,h,u),delta,EPS);
    return(grad);
}



// R_power:
double deri_R_power_wen_time(double *par, double h,double u)
{
    double R_power=par[0];
    double R_power_t=par[2];
    double scale_s=par[3];
    double scale_t=par[4];
    double sep=par[5];
    double smooth=par[6];
    
    //double EPS = 1.0e-10;
    
    double delta=sqrt(EPS)*R_power;
    
    double *par1;
    par1 = (double *) Calloc(7, double); /* idem */
    // define vector par1
    par1[0]=R_power+ delta;
    par1[1]=2;
    par1[2]=R_power_t;
    par1[3]=scale_s;
    par1[4]=scale_t;
    par1[5]=sep;
    par1[6]=smooth;
    double grad=(wen_time(par1,h,u)-wen_time(par,h,u))/delta;
    //Rprintf("%f %f %f %f\n",wen_time(par1,h,u),wen_time(par,h,u),delta,EPS);
    return(grad);
}



// R_power_t:
double deri_R_power_t_wen_time(double *par, double h,double u)
{
    double R_power=par[0];
    double R_power_t=par[2];
    double scale_s=par[3];
    double scale_t=par[4];
    double sep=par[5];
    double smooth=par[6];
    
    //double EPS = 1.0e-10;
    
    double delta=sqrt(EPS)*R_power_t;
    
    double *par1;
    par1 = (double *) Calloc(7, double); /* idem */
    // define vector par1
    par1[0]=R_power;
    par1[1]=2;
    par1[2]=R_power_t+ delta;
    par1[3]=scale_s;
    par1[4]=scale_t;
    par1[5]=sep;
    par1[6]=smooth;
    double grad=(wen_time(par1,h,u)-wen_time(par,h,u))/delta;
    //Rprintf("%f %f %f %f\n",wen_time(par1,h,u),wen_time(par,h,u),delta,EPS);
    return(grad);
}


/* END DERIVATIVES Wendland covariance */


//---------END WENDLAND FUNCTIONS-----------












/********************************************************************/
/************ gradient of the pairwise CL ****************/
/********************************************************************/
void Grad_Pair_Gauss(double rho, int *flag, double *gradcor, double *grad,
                     int *npar, double *par, double u, double v)
{
    // Initialization variables:
    //printf("%f\t%f\t%f\t%f\t%f\t%f\n",gradcor[0],gradcor[1],grad[0],grad[1],grad[2],grad[3]);
    double mean=par[0],nugget=par[1],sill=par[2];
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
    if(flag[0]==1){grad[i]=(u+v)/(a+b);i++;}
    // Derivative  respect with the nugget
    if(flag[1]==1)
    {
        grad[i]=0.5*k*(R*d-L*4*b*a-2*a*(pa-pb));i++;
    }
    // Derivative respect with the sill
    if(flag[2]==1)
    {
        grad[i]=-0.5*k*(2*(pa*a-pb*(2*sill+3*nugget)+rho*b*(pb-pn))+R*(c+2*nugget*b*rho)+2*L*rho*(ps-pn-pb));
        // Rprintf("%f \n",grad[i],mean);
        i++;
    }
    //printf("%f\t%f\t%f\t%f\t%f\t%f\n",gradcor[0],gradcor[1],grad[0],grad[1],grad[2],grad[3]);
    //printf("I: %d\n",i);
    // Derivatives with respect to the correlation parameters
    h=0;
    C=-k*sill*(R*a*b-L*d+b*c);
    //printf("%f\n",C);
    for(j=i;j<*npar;j++){grad[j]=C*gradcor[h];h++;}
    //printf("grad: %f\t%d\n",grad[j],j);
    //printf("H: %d\n",h);h++;}
    //printf("%f\t%f\t%f\t%f\t%f\t%f\n",gradcor[0],gradcor[1],grad[0],grad[1],grad[2],grad[3]);
    return;
}




/********************************************************************/
/************ derivative of correlation function ****************/
/********************************************************************/

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

// Derivatives with respect to scale of the Stable correlation model:
double DStabSc(double lag, double power, double scale, double rho)
{
    if(lag) return rho*power*pow(lag/scale,power)/scale;
    else return 0.0;
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


double CorFct(int *cormod, double h, double u, double *par)
{
    //Rprintf("%f %f %f %f %f %f %f\n",par[0],par[1],par[2],par[3],par[4],par[5],par[6]);
    double arg=0.0,  power_s=0.0, power_t=0.0;
    double rho=0.0, sep=0,  scale_s=0.0, scale_t=0;
    
    switch(*cormod) // Correlation functions are in alphabetical order
    {
            // START non-separable correlation functions:
        case 1:// Double exp:
            scale_s=par[0];
            scale_t=par[1];
            rho=CorFunStable(h,1,scale_s)*CorFunStable(u,1,scale_t);
            break;
        case 2:   //Gneiting correlation model as in (14) Gneitint (2002) with tau=1
            power_s=par[0];
            power_t=par[1];
            scale_s=par[2];
            scale_t=par[3];
            sep=par[4];
            arg=1+pow(u/scale_t, power_t);
            rho=exp(-(pow(h/scale_s, power_s))*pow(arg, -0.5*sep*power_s))/pow(arg,1);
            break;
            
        case 3:   //WendLand Time correlation model
            rho=wen_time(par,h,u);
            //Rprintf("rho: %f\n",rho);
            break;
    }
    return rho;
}
/********************************************************************/
/********************************************************************/
/********************************************************************/

//  Derivatives with respect ot the correlations parameters:
void GradCorrFct(double rho, int *cormod, int *flag,
                 double *grad, double h, double u,double *par)
{
    //Rprintf("FLAGS: %d %d %d %d %d %d %d \n",flag[0],flag[1],flag[2],flag[3],flag[4],flag[5],flag[6]);
    //Rprintf("PAR: %f %f %f %f %f %f %f\n",par[0],par[1],par[2],par[3],par[4],par[5],par[6]);
    int i=0;
    double power_s=0, power_t=0;
    double scale_s=0, scale_t=0, sep=0;
    
    switch(*cormod)// Correlation functions are in alphabetical order
    {
            //spatial gradients of correlations:
        case 1://Double Exponentil
            scale_s=par[0];
            scale_t=par[1];
            if(flag[0]==1){//spatial-scale parameter
                grad[i]=DStabSc(h,1,scale_s,rho);i++;}
            //temporal-scale parameter
            if(flag[1]==1) grad[i]=DStabSc(u,1,scale_t,rho);
            break;
        case 2://Gneiting spatio-temporal correlation
            power_s=par[0];
            power_t=par[1];
            scale_s=par[2];
            scale_t=par[3];
            sep=par[4];
            
            if(flag[0]==1){//spatial-power parameter
                grad[i]=DGneiting_pw_s(h,u,power_s,power_t,scale_s,scale_t,sep);i++;}
            if(flag[1]==1){//temporal-power parameter
                grad[i]=DGneiting_pw_t(h,u,power_s,power_t,scale_s,scale_t,sep);i++;}
            if(flag[2]==1){//spatial-scale parameter
                grad[i]=DGneiting_sc_s(h,u,power_s,power_t,scale_s,scale_t,sep);i++;}
            if(flag[3]==1){//temporal-scale parameter
                grad[i]=DGneiting_sc_t(h,u,power_s,power_t,scale_s,scale_t,sep);i++;}
            //separable parameter
            if(flag[4]==1)
                grad[i]=DGneiting_sep(h,u,power_s,power_t,scale_s,scale_t,sep);
            break;
        case 3://WendLand time
            //0:power2_s 1:power_s  2:power2_t   3:scale_s     4:scale_t   5:sep 6:smooth
            if(flag[0]==1){//power2_s parameter
                grad[i]=deri_R_power_wen_time(par, h,u);i++;}
            if(flag[1]==1){//power_s parameter
                grad[i]=0;i++;}
            if(flag[2]==1){//power2_t parameter
                grad[i]=deri_R_power_t_wen_time(par, h,u);i++;}
            if(flag[3]==1){//scale_s parameter
                grad[i]=deri_scale_s_wen_time(par, h,u);i++;}
            if(flag[4]==1){//scale_t parameter
                grad[i]=deri_scale_t_wen_time(par, h,u);i++;}
            if(flag[5]==1){//sep parameter
                grad[i]=deri_sep_wen_time(par, h,u);i++;}
            //smooth parameter
            if(flag[6]==1)
                grad[i]=deri_smooth_wen_time(par, h,u);
            break;
    }
    
}

// Computes the Geodesic distance between to coordinates:
double Dist_geodesic(double loni, double lati, double lonj, double latj)
{
    double ai, bi, aj, bj, val;
    val = 0.0;
    if (loni == lonj && lati == latj) return val;
    ai = (lati)*M_PI/180;
    bi = (loni)*M_PI/180;
    aj = (latj)*M_PI/180;
    bj = (lonj)*M_PI/180;
    val = sin(ai) * sin(aj) + cos(ai) * cos(aj) * cos(bi - bj);
    val = acos(val) *  REARTH;
    return val;
}


/**************** return a bivariate vector with maximum and minimum of a vector *********/
void Range(double *x, double *ran, int size)
{
    int i=0;
    
    ran[0] = x[0];
    ran[1] = x[0];
    
    for(i = 1; i < size; i++)
    {
        ran[0] = fmin(ran[0], x[i]);
        ran[1] = fmax(ran[1], x[i]);
    }
    return;
}
/**********************************************************************************/
int is_equal(double val1, double val2)
{
    return fabs(val1-val2)<MAXERR;
}


void SeqStep(double *x, int len, double step, double *res)
{
    int i=0;
    res[0]=x[0];
    for(i=0;i<len;i++) {
        res[i+1]=res[i]+step;
    //    Rprintf("res[i]+step %f\n",res[i]+step);
    }
    return;
}
/****************************************************************************************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/
// subsampling in space
void SetSampling_s(int ncoord,int ntime,double *coordx, double *coordy, double *data, int *npts,double *scoordx, double *scoordy, double *sdata, double xmax,double xmin, double ymax,double ymin)
{
    int i=0, j=0, f=0,h=0;
    for(i=0;i<ncoord;i++)
    {
        if((xmin<coordx[i]||is_equal(xmin,coordx[i]))&&
           (xmax>coordx[i]||is_equal(xmax,coordx[i]))&&
           (ymin<coordy[i]||is_equal(ymin,coordy[i]))&&
           (ymax>coordy[i]||is_equal(ymax,coordy[i])))
        {
            scoordx[j]=coordx[i];
            scoordy[j]=coordy[i];
            //Rprintf("%f %f\n",scoordx[j],scoordy[j]);
            for(h=0;h<ntime;h++)
            {
                sdata[f]=data[h+ntime*i];f++;
            }
            j++;
        }
    }
    *npts = j;
    return;
}
// subsampling in time
void SetSampling_t(double *data,double *sdata,int ncoord,int ntime,int wint,int k)
{
    int i=0,j=0,p=0;
    for(i=0;i<(ncoord);i++)
        for(j=(k+(ntime*i));j<(k+wint+(ntime*i));j++)
        {
            sdata[p]=data[j];p++;
            //printf("%d\t%d\t%d\t\n",p,i,j);
        }
    
    return;
}

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
    //double *gradcor, *grad, *ww;
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
    
    //Rprintf("X: %d  %f  %f  %f  %f  %f  %f\n",numintx,deltax,rangex[1],rangex[0],dimwinx,winstx,winstp[0]);
    //Rprintf("Y: %d  %f  %f  %f  %f  %f\n",numinty,deltay,rangey[1],rangey[0],dimwiny,winsty);
    //Rprintf("B. MIERDAAA\n");
    
    xgrid=(double *) R_alloc(numintx, sizeof(double));
    ygrid=(double *) R_alloc(numinty, sizeof(double));
    
    //Rprintf("C. MIERDAAA\n");
    
    //Rprintf("Y: %d  %d  %f  %f  %f  %f %f %f\n",numintx,numinty,rangex[1],rangex[0],rangey[1],rangey[0],winstx,winsty);
    SeqStep(rangex,numintx,winstx,xgrid);
    SeqStep(rangey,numinty,winsty,ygrid);
    
    //Rprintf("D. MIERDAAA\n");
    
    // matrix of means (one for each window)
    
    vector_mean= (double **) Calloc(npar[0],double *);
    for(i=0;i<npar[0];i++)
    {
        vector_mean[i]=(double *) Calloc((numintx+1)*(numinty+1),double);
    }
    //Rprintf("E. MIERDAAA\n");
    
    ww=(double *) Calloc(*npar,double);
    gradcor=(double *) Calloc(*nparc,double);
    grad=(double *) Calloc(*npar,double);
    
    nsub=0;
    n_win=numintx*numinty;   //// number of windows
    //Rprintf("%d\n",n_win);
    
    
    /*int nsens=0;
    double *sens,*sensmat;
    nsens=*npar * (*npar+1)/2;
    sens=(double *) Calloc(nsens,double);// One sensitive contribute
    sensmat=(double *) Calloc(nsens,double);// sensitivity matrix*/
    
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
    
    
    //printf("npar:%d, sill %f alfa_s %f alfa_s %f\n",npar[0],grad[2],gradcor[0],gradcor[1]);
    //printf("npar:%d, g1 %f g2 %f g3 %f  gc %f gc %f\n",npar[0],grad[0],grad[1],grad[2],gradcor[0],gradcor[1]);
    //printf("nsens:%d, 1 %f 2 %f 3 %f 4 %f 5 %f 6 %f\n",nsens,sensmat[0],sensmat[1],sensmat[2],sensmat[3],sensmat[4],sensmat[5]);
    
    //printf("nparc:%d npar:%d\n",*nparc,*npar);
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
    //Rprintf("%f\n",winstp[0]);
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
   
                DoubleExpOCL(npts,ntime,coordt,maxtime,maxdist,cormod,parcor,flagcor,flagnuis,npar,nuis,sdata,weigthed,mom_cond,dist,scoordx,scoordy,gradcor,grad,ww,local_wi,dev);
                //printf("%f\t%f\t%f\t%f\n",mom_cond[0],mom_cond[1],mom_cond[2],mom_cond[3]);
                
                
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
    //Rprintf("%f\n",winstp[0]);
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


                GneitingOCL(npts,ntime,coordt,maxtime,maxdist,cormod,parcor,flagcor,flagnuis,npar,nuis,sdata,weigthed,mom_cond,dist,scoordx,scoordy,gradcor,grad,ww,local_wi,dev);
                
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
    //Rprintf("%f\n",winstp[0]);
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


void scalar_time(int *ncoord,int *nstime,double *sublagt,int *cormod,double *parcor,int *flagcor, double *gradcor,int *flagnuis, double *grad,int *npar,double *nuis, double *sdata,int *weigthed,double *maxtime, double *ww, double *mom_cond,int *dist, double *coordx, double *coordy,double *maxdist)
{
    //Rprintf("C: R_power_s: %f R_power: %f R_power_t: %f scale_s: %f scale_t: %f sep: %f smooth: %f \n",parcor[1],parcor[0],parcor[2],parcor[3],parcor[4],parcor[5],parcor[6]);
    
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
    //Rprintf("mean: %f nugget: %f sill: %f scale_s: %f scale_t: %f\n",nuis[0],nuis[1],nuis[2],parcor[0],parcor[1]);
    //if(nuis[1]<0||nuis[2]<0 || parcor[0]<0 ||parcor[1]<0) {*res=1e100;  return;}
    //if(nuis[2]<0 ){return;}
    
    
    //Rprintf("C: R_power_s: %f R_power: %f R_power_t: %f scale_s: %f scale_t: %f sep: %f smooth: %f \n",parcor[1],parcor[0],parcor[2],parcor[3],parcor[4],parcor[5],parcor[6]);
    
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
        //Rprintf("wint[%d]: %f\tnstime: %d\n",i,sublagt[i],nstime);
    }
    nsub=floor((( (*ntime)-wint)/(wint*winstp[0])+1));
    /*if(nsub<=2)
    {
        Rprintf("Number of temporal windows must be greater than 2, check winc_t parameter. START values are printed\n");
        return;
    }*/
    //Rprintf("nsub: %d %d %f\n",nsub,wint,winstp[0]);
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
        scalar_time(ncoord,&nstime,sublagt,cormod,parcor,flagcor,gradcor,flagnuis,grad,npar,nuis,sdata,weigthed,maxtime,ww,mom_cond,dist,coordx,coordy,maxdist);
        
        /*if(mom_cond[0]<0 || mom_cond[1]<0 ||mom_cond[2]<0)
        {
            Rprintf("%f\t%f\t%f\n",mom_cond[0],mom_cond[1],mom_cond[2]);
        }*/
        
        
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
            //Rprintf("block_mean[kk]:%f\n",block_mean[kk]);
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

void SubSamp_time_DE_ocl(double *coordx, double *coordy,double *coordt, int *ncoord,int *ntime,int *cormod,double *data,int *dist, double *maxdist,double *maxtime,int *npar,double *parcor,int *nparc,double *nuis,int *nparnuis, int *flagcor, int *flagnuis,double *vari,double *a, double *b,double *winc, double *winstp,double *block_mean,int *weigthed, int *local_wi, int *dev)
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
        //printf("wint[%d]: %f\n",i,sublagt[i]);
    }
    nsub=floor((( (*ntime)-wint)/(wint*winstp[0])+1));
    //printf("nsub: %d\n",nsub);
    // vector of means for each window
    vector_mean= (double **) Calloc(npar[0],double *);
    for(i=0;i<npar[0];i++) vector_mean[i]=(double *) Calloc(nsub,double);
    //start the sub-sampling procedure:
    
    ww=(double *) Calloc(*npar,double);
    gradcor=(double *) Calloc(*nparc,double);
    grad=(double *) Calloc(*npar,double);
    //SetSampling_t(data,sdata,*ncoord,*ntime,wint,0);
    //for(i=0;i<*ncoord* *ntime;i++){printf("sdata[%d]: %f\n",i,sdata[i]);}
    
    //printf("HASTA ");
    for(i=0;i<nsub;i++)
    {//loop for the number of sub-sampling:
        mom_cond=(double *) Calloc(*npar,double);
        // set the sub-sample of the data:
        SetSampling_t(data,sdata,*ncoord,*ntime,wint,i);
        
        
        /******************************************/
        /*computing gradient in the window*/
        /******************************************/
        //scalar_time(ncoord,&nstime,sublagt,cormod,parcor,flagcor,gradcor,flagnuis,grad,npar,nuis,sdata,weigthed,maxtime,ww,mom_cond,dist,coordx,coordy,maxdist);
        //printf("nstime\n",&nstime);
        DoubleExpOCL(ncoord,&nstime,sublagt,maxtime,maxdist,cormod,parcor,flagcor,flagnuis,npar,nuis,sdata,weigthed,mom_cond,dist,coordx,coordy,gradcor,grad,ww,local_wi,dev);
        //printf("%f\t%f\t%f\t%f\n",mom_cond[0],mom_cond[1],mom_cond[2],mom_cond[3]);
        
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
        //printf("wint[%d]: %f\n",i,sublagt[i]);
    }
    nsub=floor((( (*ntime)-wint)/(wint*winstp[0])+1));
    //printf("nsub: %d\n",nsub);
    // vector of means for each window
    vector_mean= (double **) Calloc(npar[0],double *);
    for(i=0;i<npar[0];i++) vector_mean[i]=(double *) Calloc(nsub,double);
    //start the sub-sampling procedure:
    
    ww=(double *) Calloc(*npar,double);
    gradcor=(double *) Calloc(*nparc,double);
    grad=(double *) Calloc(*npar,double);
    //SetSampling_t(data,sdata,*ncoord,*ntime,wint,0);
    //for(i=0;i<*ncoord* *ntime;i++){printf("sdata[%d]: %f\n",i,sdata[i]);}
    
    //printf("HASTA ");
    for(i=0;i<nsub;i++)
    {//loop for the number of sub-sampling:
        mom_cond=(double *) Calloc(*npar,double);
        // set the sub-sample of the data:
        SetSampling_t(data,sdata,*ncoord,*ntime,wint,i);
        
        
        /******************************************/
        /*computing gradient in the window*/
        /******************************************/
        //scalar_time(ncoord,&nstime,sublagt,cormod,parcor,flagcor,gradcor,flagnuis,grad,npar,nuis,sdata,weigthed,maxtime,ww,mom_cond,dist,coordx,coordy,maxdist);
        //printf("nstime\n",&nstime);
        GneitingOCL(ncoord,&nstime,sublagt,maxtime,maxdist,cormod,parcor,flagcor,flagnuis,npar,nuis,sdata,weigthed,mom_cond,dist,coordx,coordy,gradcor,grad,ww,local_wi,dev);
        //printf("%f\t%f\t%f\t%f\n",mom_cond[0],mom_cond[1],mom_cond[2],mom_cond[3]);
        
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
        //printf("wint[%d]: %f\n",i,sublagt[i]);
    }
    nsub=floor((( (*ntime)-wint)/(wint*winstp[0])+1));
    //printf("nsub: %d\n",nsub);
    // vector of means for each window
    vector_mean= (double **) Calloc(npar[0],double *);
    for(i=0;i<npar[0];i++) vector_mean[i]=(double *) Calloc(nsub,double);
    //start the sub-sampling procedure:
    
    ww=(double *) Calloc(*npar,double);
    gradcor=(double *) Calloc(*nparc,double);
    grad=(double *) Calloc(*npar,double);
    //SetSampling_t(data,sdata,*ncoord,*ntime,wint,0);
    //for(i=0;i<*ncoord* *ntime;i++){printf("sdata[%d]: %f\n",i,sdata[i]);}
    
    //printf("HASTA ");
    for(i=0;i<nsub;i++)
    {//loop for the number of sub-sampling:
        mom_cond=(double *) Calloc(*npar,double);
        // set the sub-sample of the data:
        SetSampling_t(data,sdata,*ncoord,*ntime,wint,i);
        
        
        /******************************************/
        /*computing gradient in the window*/
        /******************************************/
        //scalar_time(ncoord,&nstime,sublagt,cormod,parcor,flagcor,gradcor,flagnuis,grad,npar,nuis,sdata,weigthed,maxtime,ww,mom_cond,dist,coordx,coordy,maxdist);
        //printf("nstime\n",&nstime);
       //Rprintf("A.MIERDAAAAA\n");
        WendOCL(ncoord,&nstime,sublagt,maxtime,maxdist,cormod,parcor,flagcor,flagnuis,npar,nuis,sdata,weigthed,mom_cond,dist,coordx,coordy,gradcor,grad,ww,local_wi,dev);
        
        //printf("%f\t%f\t%f\t%f\n",mom_cond[0],mom_cond[1],mom_cond[2],mom_cond[3]);
        
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
    //Rprintf("%f\n",winstp[0]);
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
    //Rprintf("%f\n",winstp[0]);
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
    //Rprintf("%f\n",winstp[0]);
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


// OPEN CL ///

char * getKernelSource(char *filename)
{
    FILE *file = fopen(filename, "r");
    if (!file)
    {
        //fprintf(stderr, "Error: Could not open kernel source file\n");
        printf("Error: Could not open kernel source file\n");
        //exit(EXIT_FAILURE);
    }
    fseek(file, 0, SEEK_END);
    int len = ftell(file) + 1;
    rewind(file);
    
    char *source = (char *)calloc(sizeof(char), len);
    if (!source)
    {
        //fprintf(stderr, "Error: Could not allocate memory for source string\n");
        //exit(EXIT_FAILURE);
        printf("Error: Could not allocate memory for source string\n");
    }
    fread(source, sizeof(char), len, file);
    fclose(file);
    return source;
}

/**************** sum total *********/
float sum_total(float *arr, int ngrid)
{
    float sol = 0.0;
    for (int i=0; i<ngrid; i++)
    {
        sol += arr[i];
    }
    return sol;
}

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
    
    //printf("NPARC:%d\tNPARNUIS:%d\n",nparc,nparnuis);
    
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
        printf("Invalid device index (try '--list')\n");
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
        
        printf("Error: Failed to build program executable!\n%s\n", err_code(err));
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
        SEP;
        printf("Build Log:\n%s\n", buffer);
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
        printf("Invalid device index (try '--list')\n");
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
        
        printf("Error: Failed to build program executable!\n%s\n", err_code(err));
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
        SEP;
        printf("Build Log:\n%s\n", buffer);
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
    size_t dou_par_buff = sizeof(double) * (10);
    
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
    
    
    
    for(i =0; i<5;i++){nparc += flagcor[i];}
    for(i =0; i<3;i++){nparnuis += flagnuis[i];}
    int_par = (int*)calloc((32), sizeof(int));
    dou_par = (double*)calloc((32), sizeof(double));
    
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
        printf("Invalid device index (try '--list')\n");
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
        char buffer[2048*4];
        
        printf("Error: Failed to build program executable!\n%s\n", err_code(err));
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
        SEP;
        printf("Build Log:\n%s\n", buffer);
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
    size_t dou_par_buff = sizeof(double) * (10);
    
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


// ************* **************       ***************   OPENCL
int DevOpenCL()
{
    
    int i, j;
    char* value;
    size_t valueSize;
    cl_uint platformCount;
    cl_platform_id* platforms;
    cl_uint deviceCount;
    cl_device_id* devices;
    cl_uint maxComputeUnits;
    cl_ulong long_entries;
    size_t p_size;
    cl_device_fp_config fp;
    // get all platforms
    clGetPlatformIDs(0, NULL, &platformCount);
    platforms = (cl_platform_id*) malloc(sizeof(cl_platform_id) * platformCount);
    clGetPlatformIDs(platformCount, platforms, NULL);
    cl_device_type dt;
    
    
    for (i = 0; i < platformCount; i++) {
        
        // get all devices
        clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, 0, NULL, &deviceCount);
        devices = (cl_device_id*) malloc(sizeof(cl_device_id) * deviceCount);
        clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, deviceCount, devices, NULL);
        
        //int gpu[deviceCount];
        //int countergpu=0;
        Rprintf("-------------------------------------------------------------\n");
        // for each device print critical attributes
        for (j = 0; j < deviceCount; j++) {
            Rprintf("-------------------------------------------------------------\n");
            // print device name
            clGetDeviceInfo(devices[j], CL_DEVICE_NAME, 0, NULL, &valueSize);
            value = (char*) malloc(valueSize);
            clGetDeviceInfo(devices[j], CL_DEVICE_NAME, valueSize, value, NULL);
            Rprintf("%d.\tCL_DEVICE_NAME\tDevice: \t%s\n", j, value);
            free(value);
            
            // print hardware device version
            clGetDeviceInfo(devices[j], CL_DEVICE_VERSION, 0, NULL, &valueSize);
            value = (char*) malloc(valueSize);
            clGetDeviceInfo(devices[j], CL_DEVICE_VERSION, valueSize, value, NULL);
            Rprintf("%d.%d\tCL_DEVICE_VERSION\tHardware version: \t%s\n", j, 1, value);
            free(value);
            
            // print software driver version
            clGetDeviceInfo(devices[j], CL_DRIVER_VERSION, 0, NULL, &valueSize);
            value = (char*) malloc(valueSize);
            clGetDeviceInfo(devices[j], CL_DRIVER_VERSION, valueSize, value, NULL);
            Rprintf("%d.%d\tCL_DRIVER_VERSION\tSoftware version: \t%s\n", j, 2, value);
            free(value);
            
            // print c version supported by compiler for device
            clGetDeviceInfo(devices[j], CL_DEVICE_OPENCL_C_VERSION, 0, NULL, &valueSize);
            value = (char*) malloc(valueSize);
            clGetDeviceInfo(devices[j], CL_DEVICE_OPENCL_C_VERSION, valueSize, value, NULL);
            Rprintf("%d.%d\tCL_DEVICE_OPENCL_C_VERSION\tOpenCL C version: \t%s\n", j, 3, value);
            free(value);
            
            // print parallel compute units
            clGetDeviceInfo(devices[j], CL_DEVICE_MAX_COMPUTE_UNITS,
                            sizeof(maxComputeUnits), &maxComputeUnits, NULL);
            Rprintf("%d.%d\tCL_DEVICE_MAX_COMPUTE_UNITS\tParallel compute units: \t%d\n", j, 4, maxComputeUnits);
            
            
            clGetDeviceInfo(devices[j],CL_DEVICE_GLOBAL_MEM_SIZE,sizeof(cl_ulong),&long_entries,NULL);
            Rprintf("%d.%d\tCL_DEVICE_GLOBAL_MEM_SIZE\tGlobal Memory (MB):\t%llu\n",j, 5,long_entries/1024/1024);
            clGetDeviceInfo(devices[j],CL_DEVICE_GLOBAL_MEM_CACHE_SIZE,sizeof(cl_ulong),&long_entries,NULL);
            Rprintf("%d.%d\tCL_DEVICE_GLOBAL_MEM_CACHE_SIZE\tGlobal Memory Cache (MB):\t%llu\n",j, 6,long_entries/1024/1024);
            clGetDeviceInfo(devices[j],CL_DEVICE_LOCAL_MEM_SIZE,sizeof(cl_ulong),&long_entries,NULL);
            Rprintf("%d.%d\tCL_DEVICE_LOCAL_MEM_SIZE\tLocal Memory (KB):\t%llu\n",j, 7,long_entries/1024);
            clGetDeviceInfo(devices[j],CL_DEVICE_MAX_CLOCK_FREQUENCY,sizeof(cl_ulong),&long_entries,NULL);
            Rprintf("%d.%d\tCL_DEVICE_MAX_CLOCK_FREQUENCY\tMax clock (MHz) :\t%llu\n",j, 8,long_entries);
            clGetDeviceInfo(devices[j],CL_DEVICE_MAX_WORK_GROUP_SIZE,sizeof(size_t),&p_size,NULL);
            Rprintf("%d.%d\tCL_DEVICE_MAX_WORK_GROUP_SIZE\tMax Work Group Size:\t%zu\n",j, 9,p_size);
            clGetDeviceInfo(devices[j],CL_DEVICE_MAX_WORK_ITEM_SIZES,sizeof(size_t),&p_size,NULL);
            Rprintf("%d.%d\tCL_DEVICE_MAX_WORK_ITEM_SIZES\tMax Work Item Size:\t%zu\n",j,10,p_size);
            clGetDeviceInfo(devices[j],CL_KERNEL_WORK_GROUP_SIZE,sizeof(size_t),&p_size,NULL);
            Rprintf("%d.%d\tCL_KERNEL_WORK_GROUP_SIZE\tMax kernel Work group Size:\t%zu\n",j, 11,p_size);
            clGetDeviceInfo(devices[j],CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS,sizeof(size_t),&p_size,NULL);
            Rprintf("%d.%d\tCL_DEVICE_MAX_WORK_ITEM_DIMENSIONS\tMax Dev dim:\t%zu\n",j, 12,p_size);
            clGetDeviceInfo(devices[j],CL_DEVICE_MAX_MEM_ALLOC_SIZE,sizeof(size_t),&p_size,NULL);
            Rprintf("%d.%d\tCL_DEVICE_MAX_MEM_ALLOC_SIZE\tMax Buffer size (Mb):\t%zu\n",j, 13,p_size/1000000);
            clGetDeviceInfo(devices[j],CL_DEVICE_DOUBLE_FP_CONFIG,sizeof(cl_device_fp_config),&fp,NULL);
            Rprintf("%d.%d\tSupports double precision floating-point? %s\n",j, 14,fp != 0 ? "yes" : "no");
            clGetDeviceInfo(devices[j],CL_DEVICE_TYPE,sizeof(cl_device_type),&dt,NULL);
            Rprintf("%d.%d\tCL_DEVICE_TYPE: %s\n",j, 15,dt & CL_DEVICE_TYPE_GPU ? "GPU" : "CPU");
            //clGetDeviceInfo(devices[j],CL_DEVICE_MAX_CONSTANT_ARGS,sizeof(size_t),&p_size,NULL);
            //Rprintf("%d.%d\tCL_DEVICE_MAX_CONSTANT_ARGS\tMax kernel args:\t%zu\n",j, 15,p_size);
            Rprintf("-------------------------------------------------------------\n");
            
        }
        Rprintf("-------------------------------------------------------------\n");
        free(devices);
        
    }
    
    free(platforms);
    return 0;
    
}


// ================================ Start  Create Binary Kernel

void create_binary_kernel(int *dev, char **fname)
{
    //printf("ARCHIVO fname  %s\n", *fname);
    // Context, program, build:
    cl_int err;
    cl_device_id        device;     // compute device id
    cl_context       context;       // compute context
    cl_command_queue commands;      // compute command queue
    cl_program       program;       // compute program
    //cl_kernel kernel;
    
    
    cl_uint deviceIndex = dev[0];
    
    // Get list of devices
    cl_device_id devices[MAX_DEVICES];
    unsigned numDevices = getDeviceList(devices);
    
    // Check device index in range
    if (deviceIndex >= numDevices)
    {
        Rprintf("Invalid device index (try '--list') Compilation!\n");
        //return EXIT_FAILURE;
    }
    
    device = devices[deviceIndex];
    // Create a compute context
    /*cl_platform_id platform=NULL;
     cl_context_properties ctx_properties[] = {
     CL_CONTEXT_PLATFORM, (cl_context_properties)platform, 0
     };*/
    context = clCreateContext(0, 1, &device, NULL, NULL, &err);
    checkError(err, "Creating context");
    
    // Create a command queue
    commands = clCreateCommandQueue(context, device, 0, &err);
    checkError(err, "Creating command queue");
    // Create the compute program from the source buffer
    
    
    //fclose(fp);*/
    //char CL[5];
    char f_nameCL[100];
    //strcpy(CL, ".cl");
    strcpy(f_nameCL, *fname);
    //strcat(f_nameCL, CL);
    //printf("ARCHIVO CL  %s\n", f_nameCL);
    
    
    char *kernelsource = getKernelSource(f_nameCL);
    
    program = clCreateProgramWithSource(context,1,(const char **) &kernelsource, NULL, &err);
    if(err!=CL_SUCCESS){Rprintf("Failed clCreateProgramWithSource\n");}
    
    err = clCompileProgram(program, 1, &device, "-I ./", 0, NULL, NULL, NULL, NULL);
    //err = clBuildProgram(program, 1, &device, NULL, NULL, NULL);
    if (err != CL_SUCCESS)
    {
        size_t len;
        char buffer[2048*200];
        
        Rprintf("Error: Failed to build program executable!\n%s\n", err_code(err));
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
        SEP;
        Rprintf("Build Log:\n%s\n", buffer);
        SEP;
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_STATUS, 2048*sizeof(char), buffer, &len);
        Rprintf("Build Status:\n%s\n", buffer);
        SEP;
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_OPTIONS, 2048*sizeof(char), buffer, &len);
        Rprintf("Build Options:\n%s\n", buffer);
        SEP;
        //return EXIT_FAILURE;
    }
    //kernel = clCreateKernel(program, *fname, &err);
    
    
    //size_t work_group_size;
    //err = clGetKernelWorkGroupInfo (kernel, device, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &work_group_size, NULL);
    //checkError(err, "Getting kernel work group info");
    //Rprintf("Recommended Local: %lu\n", work_group_size);
    
    
    FILE *f;
    char *binary;
    size_t binary_size;
    
    clGetProgramInfo(program, CL_PROGRAM_BINARY_SIZES, sizeof(size_t), &binary_size, NULL);
    //if(err!=CL_SUCCESS){Rprintf("Failed to get CL_PROGRAM_BINARY_SIZES\n");}
    binary = malloc(binary_size);
    clGetProgramInfo(program, CL_PROGRAM_BINARIES, binary_size, &binary, NULL);
    //if(err!=CL_SUCCESS){Rprintf("Failed to get CL_PROGRAM_BINARIES\n");}
    
    strcat(*fname, ".bin");
   // printf("ARCHIVO CL  B %s\n", *fname);
    f = fopen(*fname, "w");
    fwrite(binary, binary_size, 1, f);
    fclose(f);
    
    
    /* Finalization */
    clFlush(commands);
    clFinish(commands);
    clReleaseProgram(program);
    clReleaseCommandQueue(commands);
    clReleaseContext(context);
    //clReleaseKernel(kernel);
    Rprintf("Successful Binary Compilation!!!\n");
}

// ================================ End  Create Binary Kernel


