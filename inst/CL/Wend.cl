void Grad_Pair_Gauss(double rho, int flag0, int flag1, int flag2, double *gradcor,double *grad,int npar, double par0,double par1,double par2, double u, double v);
double int_gen(double x,double mu, double alpha,double lag,double supp);
void integr_gen(double *x, int n, void *ex);
double wendintegral(double x, double par0, double par1, double par2);
double CorFunW_gen(double lag,double R_power1,double smooth,double scale);
double wen_time(double par0,double par1,double par2,double par3,double par4,double par5,double par6, double h,double u);
double deri_scale_s_wen_time(double par0,double par1,double par2,double par3,double par4,double par5,double par6, double h,double u);
double deri_scale_t_wen_time(double par0,double par1,double par2,double par3,double par4,double par5,double par6, double h,double u);
double deri_smooth_wen_time(double par0,double par1,double par2,double par3,double par4,double par5,double par6, double h,double u);
double deri_sill_wen_time(double *par, double h,double u);
double deri_sep_wen_time(double par0,double par1,double par2,double par3,double par4,double par5,double par6, double h,double u);
double deri_R_power_wen_time(double par0,double par1,double par2,double par3,double par4,double par5,double par6, double h,double u);
double deri_R_power_t_wen_time(double par0,double par1,double par2,double par3,double par4,double par5,double par6, double h,double u);
void GradCorrFct(double rho, int flag0, int flag1, int flag2, int flag3, int flag4,int flag5,int flag6,double *grad, double h, double u,double par0,double par1,double par2,double par3,double par4,double par5,double par6);
double CorFunBohman(double lag,double scale);
double CorFunStable(double lag, double power, double scale);
double CorFct(double h, double u, double par0,double par1,double par2,double par3,double par4,double par5,double par6);
double beta (double x, double y);



// ===================================== START Integrate  =====================================//
// https://github.com/wch/r-source/blob/trunk/src/appl/integrate.c
typedef void integr_fn(double *x, int n, void *ex);
void Rdqags(integr_fn f, void *ex, double *a, double *b,
            double *epsabs, double *epsrel,
            double *result, double *abserr, int *neval, int *ier,
            int *limit, int *lenw, int *last,  int *iwork,  double *work);


static void rdqagse(integr_fn f, void *ex, double *, double *,
                    double *, double *, int *, double *, double *,
                    int *, int *, double *, double *, double *,
                    double *, int *, int *);

static void rdqk21(integr_fn f, void *ex,
                   double *, double *, double *, double *, double *, double *);
static void rdqelg(int *, double *, double *, double *, double *, int *);


static void rdqpsrt(int *, int *, int *, double *, double *, int *, int *);



void Rdqags(integr_fn f, void *ex, double *a, double *b,
            double *epsabs, double *epsrel,
            double *result, double *abserr, int *neval, int *ier,
            int *limit, int *lenw, int *last,  int *iwork,  double *work)
{
    int l1, l2, l3;
    
    //         check validity of limit and lenw.
    
    *ier = 6;
    *neval = 0;
    *last = 0;
    *result = 0.;
    *abserr = 0.;
    if (*limit < 1 || *lenw < *limit *4) {return;}
    
    //         prepare call for dqagse.
    
    l1 = *limit;
    l2 = *limit + l1;
    l3 = *limit + l2;
    
    rdqagse(f, ex, a, b, epsabs, epsrel, limit, result, abserr, neval, ier,
            work, &work[l1], &work[l2], &work[l3], iwork, last);
    
    return;
}

static void rdqagse(integr_fn f, void *ex, double *a, double *b, double *
                    epsabs, double *epsrel, int *limit, double *result,
                    double *abserr, int *neval, int *ier, double *alist,
                    double *blist, double *rlist, double *elist, int *
                    iord, int *last)
{
    // Local variables
    bool noext, extrap;
    int k,ksgn, nres;
    int ierro;
    int ktmin, nrmax;
    int iroff1, iroff2, iroff3;
    int id;
    int numrl2;
    int jupbnd;
    int maxerr;
    double res3la[3];
    double rlist2[52];
    double abseps, area, area1, area2, area12, dres, epmach;
    double a1, a2, b1, b2, defabs, defab1, defab2, oflow, uflow, resabs, reseps;
    double error1, error2, erro12, errbnd, erlast, errmax, errsum;
    
    double correc = 0.0, erlarg = 0.0, ertest = 0.0, small = 0.0;
    
    
    // ===first executable statement  dqagse
    // Parameter adjustments
    --iord;
    --elist;
    --rlist;
    --blist;
    --alist;
    
    // Function Body
    epmach = DBL_EPSILON;
    
    //            test on validity of parameters
    //            ------------------------------
    *ier = 0;
    *neval = 0;
    *last = 0;
    *result = 0.;
    *abserr = 0.;
    alist[1] = *a;
    blist[1] = *b;
    rlist[1] = 0.;
    elist[1] = 0.;
    if (*epsabs <= 0. && *epsrel < max(epmach * 50., 5e-29)) {
        *ier = 6;
        return;
    }
    
    //           first approximation to the integral
    //           -----------------------------------
    
    uflow = DBL_MIN;
    oflow = DBL_MAX;
    ierro = 0;
    rdqk21(f, ex, a, b, result, abserr, &defabs, &resabs);
    
    //           test on accuracy.
    
    dres = fabs(*result);
    errbnd = max(*epsabs, *epsrel * dres);
    *last = 1;
    rlist[1] = *result;
    elist[1] = *abserr;
    iord[1] = 1;
    if (*abserr <= epmach * 100. * defabs && *abserr > errbnd)
    {*ier = 2;}
    if (*limit == 1)
    {*ier = 1;}
    if (*ier != 0 || (*abserr <= errbnd && *abserr != resabs)
        || *abserr == 0.) {goto L140;}
    
    //           initialization
    //           --------------
    
    rlist2[0] = *result;
    errmax = *abserr;
    maxerr = 1;
    area = *result;
    errsum = *abserr;
    *abserr = oflow;
    nrmax = 1;
    nres = 0;
    numrl2 = 2;
    ktmin = 0;
    extrap = false;
    noext = false;
    iroff1 = 0;
    iroff2 = 0;
    iroff3 = 0;
    ksgn = -1;
    if (dres >= (1. - epmach * 50.) * defabs) {
        ksgn = 1;
    }
    
    //           main do-loop
    //           ------------
    
    for (*last = 2; *last <= *limit; ++(*last)) {
        
        //           bisect the subinterval with the nrmax-th largest error estimate.
        
        a1 = alist[maxerr];
        b1 = (alist[maxerr] + blist[maxerr]) * .5;
        a2 = b1;
        b2 = blist[maxerr];
        erlast = errmax;
        rdqk21(f, ex, &a1, &b1, &area1, &error1, &resabs, &defab1);
        rdqk21(f, ex, &a2, &b2, &area2, &error2, &resabs, &defab2);
        
        //           improve previous approximations to integral
        //and error and test for accuracy.
        
        area12 = area1 + area2;
        erro12 = error1 + error2;
        errsum = errsum + erro12 - errmax;
        area = area + area12 - rlist[maxerr];
        if (!(defab1 == error1 || defab2 == error2)) {
            
            if (fabs(rlist[maxerr] - area12) <= fabs(area12) * 1e-5 &&
                erro12 >= errmax * .99) {
                if (extrap){
                    ++iroff2;}
                else //if(! extrap)
                { ++iroff1;}
            }
            if (*last > 10 && erro12 > errmax)
            {++iroff3;}
        }
        rlist[maxerr] = area1;
        rlist[*last] = area2;
        errbnd = max(*epsabs, *epsrel * fabs(area));
        
        //          test for roundoff error and eventually set error flag.
        
        if (iroff1 + iroff2 >= 10 || iroff3 >= 20)
        {*ier = 2;}
        if (iroff2 >= 5)
        {ierro = 3;}
        
        //set error flag in the case that the number of subintervals equals limit.
        if (*last == *limit)
        {*ier = 1;}
        
        //          set error flag in the case of bad integrand behaviour
        //at a point of the integration range.
        
        if (max(fabs(a1), fabs(b2)) <=
            (epmach * 100. + 1.) * (fabs(a2) + uflow * 1e3)) {
            *ier = 4;
        }
        
        //           append the newly-created intervals to the list.
        
        if (error2 > error1) {
            alist[maxerr] = a2;
            alist[*last] = a1;
            blist[*last] = b1;
            rlist[maxerr] = area2;
            rlist[*last] = area1;
            elist[maxerr] = error2;
            elist[*last] = error1;
        } else {
            alist[*last] = a2;
            blist[maxerr] = b1;
            blist[*last] = b2;
            elist[maxerr] = error1;
            elist[*last] = error2;
        }
        
        //           call subroutine dqpsrt to maintain the descending ordering
        // in the list of error estimates and select the subinterval
        //with nrmax-th largest error estimate (to be bisected next).
        
        //L30:
        rdqpsrt(limit, last, &maxerr, &errmax, &elist[1], &iord[1], &nrmax);
        
        if (errsum <= errbnd)   {goto L115;}// ===jump out of do-loop
        if (*ier != 0)        {break;}
        if (*last == 2)    { // L80:
            small = fabs(*b - *a) * .375;
            erlarg = errsum;
            ertest = errbnd;
            rlist2[1] = area;    continue;
        }
        if (noext)        {continue;}
        
        erlarg -= erlast;
        if (fabs(b1 - a1) > small) {
            erlarg += erro12;
        }
        if (!extrap) {
            
            //         test whether the interval to be bisected next is the
            // smallest interval.
            
            if (fabs(blist[maxerr] - alist[maxerr]) > small) {
                continue;
            }
            extrap = true;
            nrmax = 2;
        }
        
        if (ierro != 3 && erlarg > ertest) {
            
            //           the smallest interval has the largest error.
            // before bisecting decrease the sum of the errors over the
            // larger intervals (erlarg) and perform extrapolation.
            
            id = nrmax;
            jupbnd = *last;
            if (*last > *limit / 2 + 2) {
                jupbnd = *limit + 3 - *last;
            }
            for (k = id; k <= jupbnd; ++k) {
                maxerr = iord[nrmax];
                errmax = elist[maxerr];
                if (fabs(blist[maxerr] - alist[maxerr]) > small) {
                    goto L90;
                }
                ++nrmax;
                // L50:
            }
        }
        //           perform extrapolation.  L60:
        
        ++numrl2;
        rlist2[numrl2 - 1] = area;
        rdqelg(&numrl2, rlist2, &reseps, &abseps, res3la, &nres);
        ++ktmin;
        if (ktmin > 5 && *abserr < errsum * .001) {
            *ier = 5;
        }
        if (abseps < *abserr) {
            ktmin = 0;
            *abserr = abseps;
            *result = reseps;
            correc = erlarg;
            ertest = max(*epsabs, *epsrel * fabs(reseps));
            if (*abserr <= ertest) {
                break;
            }
        }
        
        //           prepare bisection of the smallest interval.  L70:
        
        if (numrl2 == 1) {
            noext = true;
        }
        if (*ier == 5) {
            break;
        }
        maxerr = iord[1];
        errmax = elist[maxerr];
        nrmax = 1;
        extrap = false;
        small *= .5;
        erlarg = errsum;
    L90:
        ;
    }
    
    
    // L100:    set final result and error estimate.
    //        ------------------------------------
    
    if (*abserr == oflow)     {goto L115;}
    if (*ier + ierro == 0)     {goto L110;}
    if (ierro == 3)
        *abserr += correc;
    if (*ier == 0)
        *ier = 3;
    if (*result == 0. || area == 0.) {
        if (*abserr > errsum)     {goto L115;}
        if (area == 0.)     {goto L130;}
    }
    else { // L105:
        if (*abserr / fabs(*result) > errsum / fabs(area))
        {goto L115;}
    }
    
L110://        test on divergence.
    if (ksgn == -1 && max(fabs(*result), fabs(area)) <= defabs * .01) {
        goto L130;
    }
    if (.01 > *result / area || *result / area > 100. || errsum > fabs(area)) {
        *ier = 5;
    }
    goto L130;
    
L115://        compute global integral sum.
    *result = 0.;
    for (k = 1; k <= *last; ++k)
        *result += rlist[k];
    *abserr = errsum;
L130:
    if (*ier > 2)
    {L140:
        *neval = *last * 42 - 21;}
    return;
}


static void rdqelg(int *n, double *epstab, double *
                   result, double *abserr, double *res3la, int *nres)
{
    // Local variables
    int i__, indx, ib, ib2, ie, k1, k2, k3, num, newelm, limexp;
    double delta1, delta2, delta3, e0, e1, e1abs, e2, e3, epmach, epsinf;
    double oflow, ss, res;
    double errA, err1, err2, err3, tol1, tol2, tol3;
    
    
    
    // ===first executable statement  dqelg
    // Parameter adjustments
    --res3la;
    --epstab;
    
    // Function Body
    epmach = DBL_EPSILON;
    oflow = DBL_MAX;
    ++(*nres);
    *abserr = oflow;
    *result = epstab[*n];
    if (*n < 3) {
        goto L100;
    }
    limexp = 50;
    epstab[*n + 2] = epstab[*n];
    newelm = (*n - 1) / 2;
    epstab[*n] = oflow;
    num = *n;
    k1 = *n;
    for (i__ = 1; i__ <= newelm; ++i__) {
        k2 = k1 - 1;
        k3 = k1 - 2;
        res = epstab[k1 + 2];
        e0 = epstab[k3];
        e1 = epstab[k2];
        e2 = res;
        e1abs = fabs(e1);
        delta2 = e2 - e1;
        err2 = fabs(delta2);
        tol2 = max(fabs(e2), e1abs) * epmach;
        delta3 = e1 - e0;
        err3 = fabs(delta3);
        tol3 = max(e1abs, fabs(e0)) * epmach;
        if (err2 <= tol2 && err3 <= tol3) {
            //           if e0, e1 and e2 are equal to within machine
            // accuracy, convergence is assumed.
            *result = res;//        result = e2
            *abserr = err2 + err3;//    abserr = fabs(e1-e0)+fabs(e2-e1)
            
            goto L100;    // ===jump out of do-loop
        }
        
        e3 = epstab[k1];
        epstab[k1] = e1;
        delta1 = e1 - e3;
        err1 = fabs(delta1);
        tol1 = max(e1abs, fabs(e3)) * epmach;
        
        //          if two elements are very close to each other, omit
        // a part of the table by adjusting the value of n
        
        if (err1 > tol1 && err2 > tol2 && err3 > tol3) {
            ss = 1. / delta1 + 1. / delta2 - 1. / delta3;
            epsinf = fabs(ss * e1);
            
            //           test to detect irregular behaviour in the table, and
            // eventually omit a part of the table adjusting the value of n.
            
            if (epsinf > 1e-4) {
                goto L30;
            }
        }
        
        *n = i__ + i__ - 1;
        goto L50;// ===jump out of do-loop
        
        
    L30:// compute a new element and eventually adjust the value of result.
        
        res = e1 + 1. / ss;
        epstab[k1] = res;
        k1 += -2;
        errA = err2 + fabs(res - e2) + err3;
        if (errA <= *abserr) {
            *abserr = errA;
            *result = res;
        }
    }
    
    //           shift the table.
    
L50:
    if (*n == limexp) {
        *n = (limexp / 2 << 1) - 1;
    }
    
    if (num / 2 << 1 == num) {ib = 2;} else {ib = 1;}
    ie = newelm + 1;
    for (i__ = 1; i__ <= ie; ++i__) {
        ib2 = ib + 2;
        epstab[ib] = epstab[ib2];
        ib = ib2;
    }
    if (num != *n) {
        indx = num - *n + 1;
        for (i__ = 1; i__ <= *n; ++i__) {
            epstab[i__] = epstab[indx];
            ++indx;
        }
    }
    //L80:
    if (*nres >= 4) {
        // L90:
        *abserr = fabs(*result - res3la[3]) +
        fabs(*result - res3la[2]) +
        fabs(*result - res3la[1]);
        res3la[1] = res3la[2];
        res3la[2] = res3la[3];
        res3la[3] = *result;
    } else {
        res3la[*nres] = *result;
        *abserr = oflow;
    }
    
L100:// compute error estimate
    *abserr = max(*abserr, epmach * 5. * fabs(*result));
    return;
}



static void  rdqk21(integr_fn f, void *ex, double *a, double *b, double *result,
                    double *abserr, double *resabs, double *resasc)
{
    // Initialized data
    
    double wg[5] = { .066671344308688137593568809893332,
        .149451349150580593145776339657697,
        .219086362515982043995534934228163,
        .269266719309996355091226921569469,
        .295524224714752870173892994651338 };
    double xgk[11] = { .995657163025808080735527280689003,
        .973906528517171720077964012084452,
        .930157491355708226001207180059508,
        .865063366688984510732096688423493,
        .780817726586416897063717578345042,
        .679409568299024406234327365114874,
        .562757134668604683339000099272694,
        .433395394129247190799265943165784,
        .294392862701460198131126603103866,
        .14887433898163121088482600112972,0. };
    double wgk[11] = { .011694638867371874278064396062192,
        .03255816230796472747881897245939,
        .05475589657435199603138130024458,
        .07503967481091995276704314091619,
        .093125454583697605535065465083366,
        .109387158802297641899210590325805,
        .123491976262065851077958109831074,
        .134709217311473325928054001771707,
        .142775938577060080797094273138717,
        .147739104901338491374841515972068,
        .149445554002916905664936468389821 };
    
    
    // Local variables
    double fv1[10], fv2[10], vec[21];
    double absc, resg, resk, fsum, fval1, fval2;
    double hlgth, centr, reskh, uflow;
    double fc, epmach, dhlgth;
    int j, jtw, jtwm1;
    
    
    
    // ===first executable statement  dqk21
    epmach = DBL_EPSILON;
    uflow = DBL_MIN;
    
    centr = (*a + *b) * .5;
    hlgth = (*b - *a) * .5;
    dhlgth = fabs(hlgth);
    
    //           compute the 21-point kronrod approximation to
    // the integral, and estimate the absolute error.
    
    resg = 0.;
    vec[0] = centr;
    for (j = 1; j <= 5; ++j) {
        jtw = j << 1;
        absc = hlgth * xgk[jtw - 1];
        vec[(j << 1) - 1] = centr - absc;
        // L5:
        vec[j * 2] = centr + absc;
    }
    for (j = 1; j <= 5; ++j) {
        jtwm1 = (j << 1) - 1;
        absc = hlgth * xgk[jtwm1 - 1];
        vec[(j << 1) + 9] = centr - absc;
        vec[(j << 1) + 10] = centr + absc;
    }
    //f(vec, 21, ex);
    integr_gen(vec, 21, ex);
    fc = vec[0];
    resk = wgk[10] * fc;
    *resabs = fabs(resk);
    for (j = 1; j <= 5; ++j) {
        jtw = j << 1;
        absc = hlgth * xgk[jtw - 1];
        fval1 = vec[(j << 1) - 1];
        fval2 = vec[j * 2];
        fv1[jtw - 1] = fval1;
        fv2[jtw - 1] = fval2;
        fsum = fval1 + fval2;
        resg += wg[j - 1] * fsum;
        resk += wgk[jtw - 1] * fsum;
        *resabs += wgk[jtw - 1] * (fabs(fval1) + fabs(fval2));
        // L10:
    }
    for (j = 1; j <= 5; ++j) {
        jtwm1 = (j << 1) - 1;
        absc = hlgth * xgk[jtwm1 - 1];
        fval1 = vec[(j << 1) + 9];
        fval2 = vec[(j << 1) + 10];
        fv1[jtwm1 - 1] = fval1;
        fv2[jtwm1 - 1] = fval2;
        fsum = fval1 + fval2;
        resk += wgk[jtwm1 - 1] * fsum;
        *resabs += wgk[jtwm1 - 1] * (fabs(fval1) + fabs(fval2));
        // L15:
    }
    reskh = resk * .5;
    *resasc = wgk[10] * fabs(fc - reskh);
    for (j = 1; j <= 10; ++j) {
        *resasc += wgk[j - 1] * (fabs(fv1[j - 1] - reskh) +
                                 fabs(fv2[j - 1] - reskh));
        // L20:
    }
    *result = resk * hlgth;
    *resabs *= dhlgth;
    *resasc *= dhlgth;
    *abserr = fabs((resk - resg) * hlgth);
    if (*resasc != 0. && *abserr != 0.) {
        *abserr = *resasc * min(1., pow(*abserr * 200. / *resasc, 1.5));
    }
    if (*resabs > uflow / (epmach * 50.)) {
        *abserr = max(epmach * 50. * *resabs, *abserr);
    }
    return;
}



static void rdqpsrt(int *limit, int *last, int *maxerr,
                    double *ermax, double *elist, int *iord, int *nrmax)
{
    // Local variables
    int i, j, k, ido, jbnd, isucc, jupbn;
    double errmin, errmax;
    
    
    // Parameter adjustments
    --iord;
    --elist;
    
    // Function Body
    
    //           check whether the list contains more than
    // two error estimates.
    if (*last <= 2) {
        iord[1] = 1;
        iord[2] = 2;
        goto Last;
    }
    //           this part of the routine is only executed if, due to a
    // difficult integrand, subdivision increased the error
    // estimate. in the normal case the insert procedure should
    // start after the nrmax-th largest error estimate.
    
    errmax = elist[*maxerr];
    if (*nrmax > 1) {
        ido = *nrmax - 1;
        for (i = 1; i <= ido; ++i) {
            isucc = iord[*nrmax - 1];
            if (errmax <= elist[isucc])
            {break;} // out of for-loop
            iord[*nrmax] = isucc;
            --(*nrmax);
            // L20:
        }
    }
    
    //L30:       compute the number of elements in the list to be maintained
    // in descending order. this number depends on the number of
    // subdivisions still allowed.
    if (*last > *limit / 2 + 2)
    {jupbn = *limit + 3 - *last;}
    else
        jupbn = *last;
    
    errmin = elist[*last];
    
    //          insert errmax by traversing the list top-down,
    // starting comparison from the element elist(iord(nrmax+1)).
    
    jbnd = jupbn - 1;
    for (i = *nrmax + 1; i <= jbnd; ++i) {
        isucc = iord[i];
        if (errmax >= elist[isucc]) {// ===jump out of do-loop
            // L60: insert errmin by traversing the list bottom-up.
            iord[i - 1] = *maxerr;
            for (j = i, k = jbnd; j <= jbnd; j++, k--) {
                isucc = iord[k];
                if (errmin < elist[isucc]) {
                    // goto L80; ===jump out of do-loop
                    iord[k + 1] = *last;
                    goto Last;
                }
                iord[k + 1] = isucc;
            }
            iord[i] = *last;
            goto Last;
        }
        iord[i - 1] = isucc;
    }
    
    iord[jbnd] = *maxerr;
    iord[jupbn] = *last;
    
Last:// set maxerr and ermax.
    
    *maxerr = iord[*nrmax];
    *ermax = elist[*maxerr];
    return;
}


double beta (double x, double y)
{
    return( (  tgamma(x)*tgamma(y)  )/ tgamma(x+y)  );
}


// integrand  in  generalized wendland function*/
double int_gen(double x,double mu, double alpha,double lag,double supp)
{
    double res=0.0,y;
    y=lag/supp;
    res=pow(1-x,mu-1)*pow(x*x-y*y,alpha)/beta(2*alpha+1,mu);
    return (res);
}

// function generalized wendland  to integrate

void integr_gen(double *x, int n, void *ex)
{
    int i;double mu,alpha,beta,y;
    mu =    ((double*)ex)[0];  //mu
    alpha = ((double*)ex)[1];  //alpha
    beta =     ((double*)ex)[2];  //csupp
    y =     ((double*)ex)[3];  //h
    for (i=0;i<n;i++) {x[i]=int_gen(x[i],mu,alpha,y,beta);}
    return;
}

// function computing generalized wendland
double wendintegral(double x, double par0,double par1,double par2)
{
    double ex[4], lower, upper, epsabs, epsrel, result, abserr;
    int neval, ier, subdiv, lenw, last;
    subdiv = 100;
    int iwork[subdiv];
    double work[4 * subdiv];
    epsabs = pow(DBL_EPSILON, 0.25);
    epsrel = epsabs;
    lenw = 4 * subdiv;             // as instructed in WRE
    ex[0] = par0; ex[1] = par1; ex[2] = par2;ex[3]=x;
    lower=x/(ex[2]);
    upper=1;
    // Compute the integral
    if(x<=par2) {
        Rdqags(integr_gen, (void *) &ex,
               &lower, &upper, &epsabs, &epsrel, &result,
               &abserr, &neval, &ier, &subdiv, &lenw, &last, iwork, work);
        
    }else   {result=0;}
    return(result);
}



// ===================================== END Integrate  =====================================//



///********************************************************************
//************ gradient of the pairwise CL ****************
//********************************************************************



void Grad_Pair_Gauss(double rho, int flag0, int flag1, int flag2,double *gradcor, double *grad,
                     int npar, double par0,double par1,double par2, double u, double v)
{
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
    for(j=i;j<npar;j++){grad[j]=C*gradcor[h];h++;}
    return;
}

#define EPS1 1.0e-40
//#define SQE 3.162278e-30
#define SQE 1.0e-15


//---------START WENDLAND FUNCTIONS-----------





// END Wendland covariance


// START DERIVATIVES Wendland covariance

// SCALE_S:
double deri_scale_s_wen_time(double par0,double par1,double par2,double par3,double par4,double par5,double par6, double h,double u)
{
    double paro11[7];
    paro11[0]=par0;
    paro11[1]=par1;
    paro11[2]=par2;
    paro11[3]=par3;
    paro11[4]=par4;
    paro11[5]=par5;
    paro11[6]=par6;

    double delta=SQE*par3;
    
    double paro1[7];
    paro1[0]=par0;
    paro1[1]=2;
    paro1[2]=par2;
    paro1[3]=par3 + delta;
    paro1[4]=par4;
    paro1[5]=par5;
    paro1[6]=par6;
    double grad=(wen_time(paro1[0],paro1[1],paro1[2],paro1[3],paro1[4],paro1[5],paro1[6],h,u)-wen_time(paro11[0],paro11[1],paro11[2],paro11[3],paro11[4],paro11[5],paro11[6],h,u))/delta;
    return(grad);
}



// SCALE_T:
double deri_scale_t_wen_time(double par0,double par1,double par2,double par3,double par4,double par5,double par6, double h,double u)
{
    double paro11[7];
    paro11[0]=par0;
    paro11[1]=par1;
    paro11[2]=par2;
    paro11[3]=par3;
    paro11[4]=par4;
    paro11[5]=par5;
    paro11[6]=par6;
    
    //double EPS = 1.0e-10;
    
    //double delta=sqrt(EPS1)*par4;
    double delta=SQE*par4;
    
    double paro1[7];
    paro1[0]=par0;
    paro1[1]=2.0;
    paro1[2]=par2;
    paro1[3]=par3 ;
    paro1[4]=par4+ delta;
    paro1[5]=par5;
    paro1[6]=par6;
    double grad=(wen_time(paro1[0],paro1[1],paro1[2],paro1[3],paro1[4],paro1[5],paro1[6],h,u)-wen_time(paro11[0],paro11[1],paro11[2],paro11[3],paro11[4],paro11[5],paro11[6],h,u))/delta;
    //printf("CL: scale_s:%0.20f scale_s1:%0.20f\n",scale_s,scale_s+ delta);
    //printf("%f %f %f %f\n",wen_time(par1,h,u),wen_time(par,h,u),delta,EPS);
    //free(par1);
    //printf("SCALE_S CL: grad:%f\n",grad);
    return(grad);
    
}



// SMOOTH:
double deri_smooth_wen_time(double par0,double par1,double par2,double par3,double par4,double par5,double par6, double h,double u)
{
    double paro11[7];
    paro11[0]=par0;
    paro11[1]=par1;
    paro11[2]=par2;
    paro11[3]=par3;
    paro11[4]=par4;
    paro11[5]=par5;
    paro11[6]=par6;
    
    //double EPS = 1.0e-10;
    
    //double delta=sqrt(EPS1)*par6;
    double delta=SQE*par6;
    
    double paro1[7];
    paro1[0]=par0;
    paro1[1]=2.0;
    paro1[2]=par2;
    paro1[3]=par3 ;
    paro1[4]=par4;
    paro1[5]=par5;
    paro1[6]=par6+ delta;
    double grad=(wen_time(paro1[0],paro1[1],paro1[2],paro1[3],paro1[4],paro1[5],paro1[6],h,u)-wen_time(paro11[0],paro11[1],paro11[2],paro11[3],paro11[4],paro11[5],paro11[6],h,u))/delta;
    //printf("CL: scale_s:%0.20f scale_s1:%0.20f\n",scale_s,scale_s+ delta);
    //printf("%f %f %f %f\n",wen_time(par1,h,u),wen_time(par,h,u),delta,EPS);
    //free(par1);
    //printf("SCALE_S CL: grad:%f\n",grad);
    return(grad);
}


// SEP:
double deri_sep_wen_time(double par0,double par1,double par2,double par3,double par4,double par5,double par6, double h,double u)
{
    double paro11[7];
    paro11[0]=par0;
    paro11[1]=par1;
    paro11[2]=par2;
    paro11[3]=par3;
    paro11[4]=par4;
    paro11[5]=par5;
    paro11[6]=par6;
    
    //double EPS = 1.0e-10;
    
    //double delta=sqrt(EPS1)*par5;
    double delta=SQE*par5;
    
    double paro1[7];
    paro1[0]=par0;
    paro1[1]=2.0;
    paro1[2]=par2;
    paro1[3]=par3 ;
    paro1[4]=par4;
    paro1[5]=par5+ delta;
    paro1[6]=par6;
    double grad=(wen_time(paro1[0],paro1[1],paro1[2],paro1[3],paro1[4],paro1[5],paro1[6],h,u)-wen_time(paro11[0],paro11[1],paro11[2],paro11[3],paro11[4],paro11[5],paro11[6],h,u))/delta;
    //printf("CL: scale_s:%0.20f scale_s1:%0.20f\n",scale_s,scale_s+ delta);
    //printf("%f %f %f %f\n",wen_time(par1,h,u),wen_time(par,h,u),delta,EPS);
    //free(par1);
    //printf("SCALE_S CL: grad:%f\n",grad);
    return(grad);
}



// R_power:
double deri_R_power_wen_time(double par0,double par1,double par2,double par3,double par4,double par5,double par6, double h,double u)
{
    double paro11[7];
    paro11[0]=par0;
    paro11[1]=par1;
    paro11[2]=par2;
    paro11[3]=par3;
    paro11[4]=par4;
    paro11[5]=par5;
    paro11[6]=par6;
    
    //double EPS = 1.0e-10;
    
    //double delta=sqrt(EPS1)*par0;
    double delta=SQE*par0;
    
    double paro1[7];
    paro1[0]=par0+ delta;
    paro1[1]=2.0;
    paro1[2]=par2;
    paro1[3]=par3 ;
    paro1[4]=par4;
    paro1[5]=par5;
    paro1[6]=par6;
    double grad=(wen_time(paro1[0],paro1[1],paro1[2],paro1[3],paro1[4],paro1[5],paro1[6],h,u)-wen_time(paro11[0],paro11[1],paro11[2],paro11[3],paro11[4],paro11[5],paro11[6],h,u))/delta;
    //printf("CL: scale_s:%0.20f scale_s1:%0.20f\n",scale_s,scale_s+ delta);
    //printf("%f %f %f %f\n",wen_time(par1,h,u),wen_time(par,h,u),delta,EPS);
    //free(par1);
    //printf("SCALE_S CL: grad:%f\n",grad);
    return(grad);
}



// R_power_t:
double deri_R_power_t_wen_time(double par0,double par1,double par2,double par3,double par4,double par5,double par6, double h,double u)
{
    double paro11[7];
    paro11[0]=par0;
    paro11[1]=par1;
    paro11[2]=par2;
    paro11[3]=par3;
    paro11[4]=par4;
    paro11[5]=par5;
    paro11[6]=par6;
    
    //double EPS = 1.0e-10;
    
   // double delta=sqrt(EPS1)*par2;
    double delta=SQE*par2;
    
    double paro1[7];
    paro1[0]=par0;
    paro1[1]=2.0;
    paro1[2]=par2+ delta;
    paro1[3]=par3 ;
    paro1[4]=par4;
    paro1[5]=par5;
    paro1[6]=par6;
    double grad=(wen_time(paro1[0],paro1[1],paro1[2],paro1[3],paro1[4],paro1[5],paro1[6],h,u)-wen_time(paro11[0],paro11[1],paro11[2],paro11[3],paro11[4],paro11[5],paro11[6],h,u))/delta;
    //printf("CL: scale_s:%0.20f scale_s1:%0.20f\n",scale_s,scale_s+ delta);
    //printf("%f %f %f %f\n",wen_time(par1,h,u),wen_time(par,h,u),delta,EPS);
    //free(par1);
    //printf("SCALE_S CL: grad:%f\n",grad);
    return(grad);
}


// END DERIVATIVES Wendland covariance


//---------END WENDLAND FUNCTIONS-----------

//  Derivatives with respect ot the correlations parameters:


void GradCorrFct(double rho, int flag0, int flag1, int flag2, int flag3, int flag4,int flag5,int flag6,  double *grad, double h, double u,double par0,double par1,double par2,double par3,double par4,double par5,double par6)
{
    
    
    int i=0;
 
    
    if(flag0==1){//power2_s parameter NO
        grad[i]=deri_R_power_wen_time(par0,par1,par2,par3,par4,par5,par6,h,u);i++;}
    if(flag1==1){//power_s parameter NO
        grad[i]=0.0;i++;}
    if(flag2==1){//power2_t parameter NO
        grad[i]=deri_R_power_t_wen_time(par0,par1,par2,par3,par4,par5,par6, h,u);i++;}
    if(flag3==1){//scale_s parameter SI
        grad[i]=deri_scale_s_wen_time(par0,par1,par2,par3,par4,par5,par6, h,u);i++;}
    if(flag4==1){//scale_t parameter SI
        grad[i]=deri_scale_t_wen_time(par0,par1,par2,par3,par4,par5,par6, h,u);i++;}
    if(flag5==1){//sep parameter NO
        grad[i]=deri_sep_wen_time(par0,par1,par2,par3,par4,par5,par6, h,u);i++;}
    //smooth parameter SI
    if(flag6==1)
    {grad[i]=deri_smooth_wen_time(par0,par1,par2,par3,par4,par5,par6, h,u);
        //printf("CL param: %f  %f  %f  %f  %f  %f  %f\n",parodia[0],parodia[1],parodia[2],parodia[3],parodia[4],parodia[5],parodia[6]);
    }
    
    return;
}


//
//correlation  models
//
double CorFunBohman(double lag,double scale)
{
    double rho=0.0,x=0.0;
    x=lag/scale;
    if(x<=1)
    {
        if (x>0.0) rho=(1-x)*(sin(2*M_PI*x)/(2*M_PI*x))+(1-cos(2*M_PI*x))/(2*M_PI*M_PI*x);
        else   rho=1;
    }
    else rho=0.0;
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



// generalized wendland function
double CorFunW_gen(double lag,double R_power1,double smooth,double scale)  // mu alpha beta
{
    
    //Rprintf("PARAM A: %f %f %f %f\n",lag,R_power1,smooth,scale);
    double rho=0.0,x=0.0;
    if(smooth==0) {
        
        x=lag/scale;
        if(x<=1) rho=pow(1-x,R_power1);
        else rho=0;
        return(rho);
    }
    
    if(smooth==1) {
        
        x=lag/scale;
        if(x<=1) rho=pow(1-x,R_power1+1)*(1+x*(R_power1+1));
        else rho=0;
        return(rho);
    }
    if(smooth==2) {
        
        x=lag/scale;
        if(x<=1) rho=pow(1-x,R_power1+2)*(1+x*(R_power1+2)+x*x*(R_power1*R_power1 +4*R_power1 +3 )/3  );
        else rho=0;
        return(rho);
    }
     
     if(smooth>0) {
         x=lag;
         //double param0=R_power1;double param1=smooth;double param2=scale;  //mu,alpha //beta
         //rho=wendintegral(x,param0,param1,param2);
         //printf("%f %f %f %f \n",x,R_power1,smooth,scale);
         rho=wendintegral(x,R_power1,smooth,scale);
     }
    return(rho);
}


double wen_time(double par0,double par1,double par2,double par3,double par4,double par5,double par6, double h,double u)
{
    double R_power_s=2.0;
    double R_power=par0;
    double R_power_t=par2;
    double scale_s=par3;
    double scale_t=par4;
    double sep=par5;
    double smooth=par6;
    double arg=pow(1+pow(h/scale_s,R_power_s/2),-1/(R_power_s/2));
    double
    rho=pow(arg,R_power)*CorFunW_gen(u,R_power_t,smooth,scale_t*pow(arg,sep));
    return(rho);
}

double CorFct(double h, double u, double par0,double par1,double par2,double par3,double par4,double par5,double par6)
{
    double rho=0.0;
    
    rho=wen_time(par0,par1,par2,par3,par4,par5,par6,h,u);
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
    
    int flagcor5 = int_par[16];
    int flagcor6 = int_par[17];
    
    double grad[npar];
    double gradcor[nparc];
    
    
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
    
    double parcor5 = dou_par[10];//parcor
    double parcor6 = dou_par[11];//parcor
    
    int l = get_global_id(0);
    int t = get_global_id(1);
    
    int ls = get_global_size(0);
    int ms = get_global_size(1);
    
    int m=0,v =0;
    double lagt=0.0f,lags=0.0, rho=0.0,weights=0.0;
    double ww0 =0.0,ww1 =0.0,ww2 =0.0,ww3 =0.0;
    
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
    //int j = (ntime*gidx+gidy);
    bool isValid = true;
    
    if(l >= npts) isValid = false;
    
    if(t >= ntime) isValid = false;
    
    if(isValid)
        
    {
        mom_cond0[i] = 0.0;mom_cond1[i] = 0.0;mom_cond2[i] = 0.0;mom_cond3[i] = 0.0;
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
                        
                    rho=CorFct(0.0,lagt,parcor0,parcor1,parcor2,parcor3,parcor4,parcor5,parcor6);
                        //Computing the gradient of the corr parameters
                        GradCorrFct(rho,flagcor0,flagcor1,flagcor2,flagcor3,flagcor4,flagcor5,flagcor6,gradcor,0,lagt,parcor0,parcor1,parcor2,parcor3,parcor4,parcor5,parcor6);
                        //Compute the gradient of the composite likelihood:
                        Grad_Pair_Gauss(rho,flagnuis0,flagnuis1,flagnuis2,gradcor,grad,npar,nuis0,nuis1,nuis2,sdata[(t+ntime*l)],sdata[(v+ntime*l)]);
        
                        ww0 =1.0,ww1 =1.0,ww2 =1.0,ww3 =1.0;
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
                        rho=CorFct(lags,lagt,parcor0,parcor1,parcor2,parcor3,parcor4,parcor5,parcor6);
                        //Computing the gradient of the corr parameters
                       //printf("B: CL gradcor: %f\t%f\t%f\n\n",gradcor[0],gradcor[1],gradcor[2]);
                       
                        if(isnan(rho))
                        {
                            printf("++++\n rho: %f\n",rho);
                            printf("PARAMS \nlags %f,lagt %f,parcor0 %f,parcor1 %f,parcor2 %f,parcor3 %f,parcor4 %f,parcor5 %f,parcor6 %f\n",lags,lagt,parcor0,parcor1,parcor2,parcor3,parcor4,parcor5,parcor6);
                        }
                        GradCorrFct(rho,flagcor0,flagcor1,flagcor2,flagcor3,flagcor4,flagcor5,flagcor6,gradcor,0,lagt,parcor0,parcor1,parcor2,parcor3,parcor4,parcor5,parcor6);
                        
                        //Compute the gradient of the composite likelihood:
                    Grad_Pair_Gauss(rho,flagnuis0,flagnuis1,flagnuis2,gradcor,grad,npar,nuis0,nuis1,nuis2,sdata[(t+ntime*l)],sdata[(v+ntime*m)]);
                    
                        
                        ww0 =1.0,ww1 =1.0,ww2 =1.0,ww3 =1.0;
                        mom_cond0[i]+=ww0*grad[0];
                        mom_cond1[i]+=ww1*grad[1];
                        mom_cond2[i]+=ww2*grad[2];
                        mom_cond3[i]+=ww3*grad[3];
                        //printf("B: grad[0] %f,grad[1] %f,grad[2] %f,grad[3] %f\n",grad[0],grad[1],grad[2],grad[3]);
                    }
                }
            
            }
            
        }
    }
}
