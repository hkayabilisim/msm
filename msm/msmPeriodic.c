#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <float.h>
#include "msmLibrary.h"
#define MYPI 3.141592653589793238462643

double Bspline(int k,double u) {
    if (k==1) {
        if (u >= 0 && u < 1)
            return 1.0;
        else
            return 0.0;
    } else {
        return u/(double)(k-1) * Bspline(k-1,u) + (k-u)/(double)(k-1) * Bspline(k-1,u-1.0);
    }
}

/* A non-optimized implementation of gamma(rho) */
double gama(double rho,int p) {
    double out = 0.0;
    if (rho >= 1.0)
        out = 1.0 / rho ;
    else {
        double  rho2 = rho * rho ;
        for (int k = 0 ; k < p ; k++) {
            double binom = tgamma(0.5)/(tgamma(k+1)*tgamma(0.5-k)) ;
            out += binom * pow(rho2-1.0 , k);
        }
    }
    return out;
}

double g1(double r,double a,int p) {
    return (1.0 / a) * gama(r / a , p) ;
}
double g0(double r,double a,int p) {
    return 1.0 / r - g1(r,a,p) ;
}
double gLR(double r, double beta) {
    if (r < DBL_EPSILON)
        return 2.0 * beta / sqrt(MYPI) ;
    else
        return erf(beta * r) / r ;
}
double g1star(double r,double a,int p,double beta) {
    return g1(r,a,p)  - gLR(r,beta) ;
}

msm_stats msmperiodic(int N,double *r,double *q,int p,double abar,int mu,double *uappr,double *fappr,
                      double Ax, double Ay, double Az)
{
    msm_stats stats;
    clock_t begin;
    double rmin[MAXDIM],rmax[MAXDIM];
    double Ox,Oy,Oz,hx,hy,hz;
    double wmax,wmin,a,aL;
    double **qgrd,**edir,**elng,**edir_part1,**edir_part2,**elng_part1,**elng_part2;
    double *flong,*fshort;
    double potential_short_real;
    int m,n,v,i,j,k,l,L,kx,ky,kz,nx,ny,nz,mx,my,mz,px,py,pz,kk,jj,ii;
    int slsizex[LMAX],slsizey[LMAX],slsizez[LMAX];
    int slminx[LMAX] ,slminy[LMAX], slminz[LMAX];
    int slmaxx[LMAX] ,slmaxy[LMAX], slmaxz[LMAX];
    int Mx,My,Mz, MLx,MLy,MLz, ML2x,ML2y,ML2z,Mlx[LMAX], Mly[LMAX], Mlz[LMAX],DIM;
    double htilde,volume;
    const struct MasterBaseFunction *base,*based;
    double PHI[MAXDIM][MAX_POLY_DEGREE],PHID[MAXDIM][MAX_POLY_DEGREE];
    double PHIATKNOTS[2][3] = {{4.0/6.0 , 1./6. , 0} , /* p = 4 */
        {66./120., 26./120. , 1./120.}}; /* p = 6 */

    /* TODO: remove it before release the code */
     for (i = 0 ; i < N * 3 ; i++)
        r[i] += 0.5 ; 

    double beta = 1.0; /* S(r) = erf(beta * r) */

    /* TODO: section: Initialization */
    clock_t total_start = clock();

    begin = clock();
    /* Base function */
    base  = & CubicBSpline;
    based = & CubicBSplineD;
    if (p == 4)
    {
        base  = & CubicBSpline;
        based = & CubicBSplineD;
    } else if ( p == 6)
    {
        base  = & QuinticBSpline ;
        based = & QuinticBSplineD;
    }

    stats.N = N;
    stats.p = p;

    if (p != 4)
    {
        fprintf(stderr,"p must be 4\n");
        return stats;
    }

    int p2 = p/2;

    /* The minimum box containing all the particles */
    for ( DIM = DIMX; DIM < MAXDIM ; DIM++)
    {
        rmin[DIM] =  VERYLARGENUMBER;
        rmax[DIM] = -VERYLARGENUMBER;
        for ( i = 0 ; i < N ; i++)
        {
            if (r[MAXDIM*i + DIM ] < rmin[DIM]) rmin[DIM] = r[MAXDIM*i+DIM];
            if (r[MAXDIM*i + DIM ] > rmax[DIM]) rmax[DIM] = r[MAXDIM*i+DIM];
        }
    }
    double rminx = rmin[DIMX]; double rmaxx = rmax[DIMX] ;
    double rminy = rmin[DIMY]; double rmaxy = rmax[DIMY] ;
    double rminz = rmin[DIMZ]; double rmaxz = rmax[DIMZ] ;
    
    /* TODO: remove it before releasing the code */
    rminx = 0.0 ; rmaxx = 1.0;
    rminy = 0.0 ; rmaxy = 1.0;
    rminz = 0.0 ; rmaxz = 1.0;
    
    stats.rminx = rminx;
    stats.rminy = rminy;
    stats.rminz = rminz;
    stats.rmaxx = rmaxx;
    stats.rmaxy = rmaxy;
    stats.rmaxz = rmaxz;

    stats.Ax = Ax ;
    stats.Ay = Ay ;
    stats.Az = Az ;

    wmin =  VERYLARGENUMBER;
    wmax = -VERYLARGENUMBER;
    for (DIM = DIMX; DIM < MAXDIM ; DIM++)
    {
        if (rmax[DIM] > wmax) wmax = rmax[DIM];
        if (rmin[DIM] < wmin) wmin = rmin[DIM];
    }

    volume = Ax * Ay * Az ;
    stats.volume = volume;
    htilde = pow(volume/N,1.0/3.0);
    //printf("...htilde:%8.4f\n",htilde);
    stats.h = htilde;
    Ox = rminx ; /* - (Ax - (rmaxx - rminx))  ;  TODO: think about margin (0.5) */
    Oy = rminy ; /* - (Ay - (rmaxy - rminy))  ;*/
    Oz = rminz ; /* - (Az - (rmaxz - rminz))  ;*/

    Mx = chooseM(htilde,Ax);
    My = chooseM(htilde,Ay);
    Mz = chooseM(htilde,Az);

    /* TODO: remove it */
    /* Mx = 2; My = 2 ; Mz = 2; */

    hx = Ax / Mx ;
    hy = Ay / My ;
    hz = Az / Mz ;

    //printf("...Mx:%2d My:%2d Mz:%2d\n",Mx,My,Mz);
    //printf("...hx:%8.4f hy:%8.4f hz:%8.4f\n",hx,hy,hz);

    stats.Mx = Mx;
    stats.My = My;
    stats.Mz = Mz;
    stats.hx = hx ;
    stats.hy = hy ;
    stats.hz = hz ;

    L = chooseL(Mx,My,Mz,N);
    stats.L = L;
    //printf("...L:%d\n",L);

    Mlx[0] = Mx;
    Mly[0] = My;
    Mlz[0] = Mz;
    for (l = 1 ; l <= L ; l++) {
        Mlx[l] = Mx / pow(2,l-1) ;
        Mly[l] = My / pow(2,l-1) ;
        Mlz[l] = Mz / pow(2,l-1) ;
    }
    for (l = 0 ; l <= L ; l++) {
        stats.Mlx[l] = Mlx[l] ;
        stats.Mly[l] = Mly[l] ;
        stats.Mlz[l] = Mlz[l] ;
    }
    MLx = Mlx[L] ;
    MLy = Mly[L] ;
    MLz = Mlz[L] ;
    ML2x = MLx/2;
    ML2y = MLy/2;
    ML2z = MLz/2;
    stats.MLx = MLx;
    stats.MLy = MLy;
    stats.MLz = MLz;
    stats.ML2x = ML2x;
    stats.ML2y = ML2y;
    stats.ML2z = ML2z;

    flong  = (double *)calloc(MAXDIM*N, sizeof(double) );
    fshort = (double *)calloc(MAXDIM*N, sizeof(double));

    qgrd = (double **)calloc(L,sizeof(double *));
    elng = (double **)calloc(L,sizeof(double *));
    elng_part1 = (double **)calloc(L,sizeof(double *));
    elng_part2 = (double **)calloc(L,sizeof(double *));

    edir = (double **)calloc(L,sizeof(double *));
    edir_part1 = (double **)calloc(L,sizeof(double *));
    edir_part2 = (double **)calloc(L,sizeof(double *));
    for (l = 0; l < L ; l++)
    {
        qgrd[l] = (double *)calloc(Mlx[l]*Mly[l]*Mlz[l],sizeof(double));
        elng[l] = (double *)calloc(Mlx[l]*Mly[l]*Mlz[l],sizeof(double));
        elng_part1[l] = (double *)calloc(Mlx[l]*Mly[l]*Mlz[l],sizeof(double));
        elng_part2[l] = (double *)calloc(Mlx[l]*Mly[l]*Mlz[l],sizeof(double));
        edir[l] = (double *)calloc(Mlx[l]*Mly[l]*Mlz[l],sizeof(double));
        edir_part1[l] = (double *)calloc(Mlx[l]*Mly[l]*Mlz[l],sizeof(double));
        edir_part2[l] = (double *)calloc(Mlx[l]*Mly[l]*Mlz[l],sizeof(double));
    }

    /* Cutoff */
    a = abar*htilde;
    /* TODO: remove the following line */
    /* a = abar * hx ; */

    aL = pow(2,L)*a;
    //printf("...a:%8.3f\n",a);
    //printf("...abar:%8.3f\n",abar);

    stats.alpha = abar;
    stats.a = a;
    stats.s = p-1;
    stats.mu = mu;
    stats.aL = aL;
    stats.cputime_initialization = (double)(clock()-begin)/CLOCKS_PER_SEC;

    /* Choose beta s.t erfc(beta r)/r < epsilon when r = aL
    beta = choosebeta(aL);
    Bob's new idea
    beta = sqrt(p*log(abar))/aL ; */
    double epsilon = 1e-12;
    double hL = pow(2,L-1)*hx;

    beta = choosebetaNew(epsilon, aL,hL);
    
    /*  B-spline interpolation of sinusoidal (Page 12)
     *  c_alpha(k_alpha) coefficients where alpha = x,y,z
     *  and k_alpha */
    begin = clock();
    double *cvecx;
    double *cvecy;
    double *cvecz;
    double kmax_double = choosekmax(epsilon, beta,hL);
    int kmax = ceil(fabs(kmax_double));
    stats.kmax = kmax;
    //printf("...... beta: %20.12e kmax: %20.12e --> %d\n",beta,kmax_double,kmax);
    {
        /* int kmax = maxOf3Integers(ML2x,ML2y,ML2z) + 1; */
        cvecx = (double *)calloc(kmax+1,sizeof(double));
        cvecy = (double *)calloc(kmax+1,sizeof(double));
        cvecz = (double *)calloc(kmax+1,sizeof(double));

        for ( k = 0 ; k < kmax + 1 ; k++) {
            double cx = PHIATKNOTS[p2-2][0];
            double cy = PHIATKNOTS[p2-2][0];
            double cz = PHIATKNOTS[p2-2][0];
            for ( m = 1 ; m <= p2 - 1 ; m++) {
                cx += 2.0*cos(2*PI*k*m/(double)MLx) * PHIATKNOTS[p2-2][m];
                cy += 2.0*cos(2*PI*k*m/(double)MLy) * PHIATKNOTS[p2-2][m];
                cz += 2.0*cos(2*PI*k*m/(double)MLz) * PHIATKNOTS[p2-2][m];
            }
            cvecx[k] = 1./cx;
            cvecy[k] = 1./cy;
            cvecz[k] = 1./cz;
        }
        stats.cvecx = cvecx;
        stats.cvecy = cvecy;
        stats.cvecz = cvecz;
        stats.clenx = ML2x + 1;
        stats.cleny = ML2y + 1;
        stats.clenz = ML2z + 1;
    }
    
    stats.beta = beta;
    /* Correction for nonzero net charge systems (Eq. 3) */
    double csr = PI / (beta*beta*volume);
    stats.csr = csr;
    /********************
     * Fourier Coefficients X(k;A)
     *
     * Equation 4.
     ********************/
    int chinx = 2*kmax+1;
    int chiny = 2*kmax+1;
    int chinz = 2*kmax+1;
    int chisize = chinx * chiny * chinz ;
    stats.chinx = chinx;
    stats.chiny = chiny;
    stats.chinz = chinz;
    stats.chisize = chisize;

    /* Accessing a 3-d matrix with a 1-d array */
    double *chi = (double *)calloc(chisize,sizeof(double));
    double chisum  = 0.0;
    for ( kz = 0 ; kz < chinz ; kz++) {
        int kzoff = chinx * chiny * kz;
        int kzz = kz - kmax ;
        double k2x = kzz * kzz / (Az * Az) ; /*  k2 will be |inv(A)*k|^2 */
        for ( ky = 0 ; ky < chiny ; ky++) {
            int kyoff = kzoff + chinx * ky;
            int kyy = ky - kmax;
            double k2y = k2x + kyy * kyy / (Ay * Ay);
            for ( kx = 0 ; kx < chinx ; kx++) {
                int kxoff = kyoff + kx;
                int kxx = kx - kmax;
                if (kxx == 0 && kyy == 0 && kzz == 0)
                    chi[kxoff] = 0;
                else {
                    double k2 = k2y + kxx * kxx / (Ax * Ax);
                    chi[kxoff] = exp(-PI * PI * k2 / (beta * beta)) / (PI * k2 * volume) ;
                }
                chisum += chi[kxoff];
            }
        }
    }
    stats.chi    = chi;
    stats.chisum = chisum;

    /********************
     * Stencil for L+1
     * (Fourier stencil)
     * Equation 19
     ********************/
    //printf("... Top-level stencil\n");

    int sLnx = 2*(MLx-1) + 1 ;
    int sLny = 2*(MLy-1) + 1 ;
    int sLnz = 2*(MLz-1) + 1 ;
    int sLsize = sLnx * sLny * sLnz ;
    stats.sLnx = sLnx;
    stats.sLny = sLny;
    stats.sLnz = sLnz;
    stats.sLsize = sLsize;
    double *sL = (double *)calloc(sLsize,sizeof(double));
    /* TODO: speedup: use the symmetry of stencil L */
    for ( nz = - (MLz-1) ; nz <= (MLz-1) ; nz++) {
        int nzoff = sLnx * sLny * (nz + MLz - 1);
        for ( ny = - (MLy-1) ; ny <= (MLy-1) ; ny++) {

            int nyoff = nzoff + sLnx * (ny + MLy - 1);
            for ( nx = - (MLx-1) ; nx <= (MLx-1) ; nx++) {
                /* double tic = clock(); */
                int nxoff = nyoff + (nx + MLx - 1);
                /* Summation in Equation 19 */
                double chisum = 0.0;
                for ( kz = - kmax ; kz < kmax + 1 ; kz++) {
                    int kzoff = chinx * chiny * ( kz + kmax );
                    double cprodz = cvecz[abs(kz)]; /* c(k) part */
                    double eprodz = kz*nz/(double)MLz;   /* exponential part */
                    for ( ky = - kmax ; ky < kmax + 1 ; ky++) {
                        int kyoff = kzoff + chinx * ( ky + kmax);
                        double cprody = cprodz * cvecy[abs(ky)];
                        double eprody = eprodz + ky*ny/(double)MLy;
                        for ( kx = - kmax ; kx < kmax + 1 ; kx++) {
                            int kxoff = kyoff + (kx + kmax);
                            double cprodx = cprody * cvecx[abs(kx)];
                            double eprodx = eprody + kx*nx/(double)MLx;
                           
                                                                
                            chisum += chi[kxoff] * cprodx * cprodx * cos(2*PI*eprodx);
                            /*if (nz == 0 && ny == 0 && nx == 0)
                                printf("%2d %2d %2d chi:%25.16e c:%25.16e exp:%25.16e sum:%25.16e\n",
                                       kz,ky,kx,chi[kxoff],cprodx * cprodx,cos(2*PI*eprodx),chisum); */
                        }
                    }
                }
                /* double toc = (double)(clock()-tic)/CLOCKS_PER_SEC;
                printf("... ... nx [%d <-- %d --> %d] (%f)\n",- (MLx-1) ,nx,(MLx-1),toc); */
                //printf("%2d %2d %2d --> %25.16e\n",nx,ny,nz,chisum);
                sL[nxoff] = chisum;
            }
        }
    }
    
    /*
    for ( k = 0 ; k < sLnz ; k++) {
        printf("Slice z=%d\n",k);
        for ( j = 0 ; j < sLny ; j++) {
            for ( i = 0 ; i < sLnx ; i++) {
                printf("%9.4f ",sL[k*sLny*sLnx + j*sLnx + i]) ;
            }
            printf("\n");
        }
    } */

    stats.cputime_stencilL = (double)(clock()-begin)/CLOCKS_PER_SEC;

    /********************
     * Short-range stencils for < L + 1
     *
     * Equation 26
     TODO: section: Stencil short-range
     ********************/
    
    begin = clock();
    double **sl;
    sl = (double **)calloc(L,sizeof(double *));
    /* TODO: remove before production */
    double **sl_part1, **sl_part2;
    sl_part1 = (double **)calloc(L,sizeof(double *));
    sl_part2 = (double **)calloc(L,sizeof(double *));


    //double ERFBETA0 = 2*beta/sqrt(PI);
    for (l = 0 ; l < L  ; l++)
    {
        //printf("... Stencil short-range : %d < %d\n",l,L);

        double al = pow(2.0,l+1)*a;
        double alminus1 = pow(2.0,l)*a;
        double al_1 = 1./al;

        double khatx = (int)ceil(2.0*a/hx);
        double khaty = (int)ceil(2.0*a/hy);
        double khatz = (int)ceil(2.0*a/hz);

        int mlx    = Mlx[l+1] ;
        int mly    = Mly[l+1] ;
        int mlz    = Mlx[l+1] ;

        double mlp1x = mlx/2.0;
        double mlp1y = mly/2.0;
        double mlp1z = mlz/2.0;
        
        int minx  ;
        int miny  ;
        int minz  ;

        int maxx  ;
        int maxy  ;
        int maxz  ;

        if (khatx > mlp1x) { minx = floor(-mlp1x) + 1; maxx = floor(mlp1x) ; } else { minx = -khatx + 1; maxx = + khatx - 1 ; }
        if (khaty > mlp1y) { miny = floor(-mlp1y) + 1; maxy = floor(mlp1y) ; } else { miny = -khaty + 1; maxy = + khaty - 1 ; }
        if (khatz > mlp1z) { minz = floor(-mlp1z) + 1; maxz = floor(mlp1z) ; } else { minz = -khatz + 1; maxz = + khatz - 1 ; }
        
        int sizex = maxx - minx + 1 ;
        int sizey = maxy - miny + 1 ;
        int sizez = maxz - minz + 1 ;

        slminx[l] = minx ;
        slminy[l] = miny ;
        slminz[l] = minz ;
        slmaxx[l] = maxx ;
        slmaxy[l] = maxy ;
        slmaxz[l] = maxz ;
        slsizex[l] = sizex ;
        slsizey[l] = sizey ;
        slsizez[l] = sizez ;

        stats.slsizex[l] = slsizex[l] ;
        stats.slsizey[l] = slsizey[l] ;
        stats.slsizez[l] = slsizez[l] ;

        int size  = sizex * sizey * sizez ;
        sl[l] = (double *)calloc(size,sizeof(double));
        /* TODO: remove */
        sl_part1[l] = (double *)calloc(size,sizeof(double));
        sl_part2[l] = (double *)calloc(size,sizeof(double));


        /*  This is used to understant */
        LinkedList *sllist = msm4g_linkedlist_new();
        /* int redundant = 0;
         int nonredundant = 0; */
        /*
         N LOOP, Stencil
         */
        for ( nz = minz ; nz <= maxz ; nz++ ) {
            int nzoff = sizex * sizey * (nz - minz) ;
            //int nzt = nz + minz;
            for ( ny = miny ; ny <= maxy  ; ny++ ) {
                int nyoff = nzoff + sizex * (ny - miny) ;
                //int nyt = ny + miny ;
                for ( nx = minx ; nx <= maxx ; nx++ ) {
                    int nxoff = nyoff + (nx - minx);
                    //int nxt = nx + minx;
                    //int hash = nxt*nxt + nyt*nyt + nzt*nzt ;
        
                    //int recordnumber = msm4g_linkedlist_search(sllist, hash);
                    /* if exists */
                    //if ( recordnumber > -1) {
                    //    double value = msm4g_linkedlist_getvalue(sllist, recordnumber);
                    //    sl[l][nxoff] = value;
                    //    continue;
                    // }
                   
                    /*
                     K LOOP, Interpolation operators Equation 28
                     Range of p: infinite norm is used i.e.
                     Take p s.t. |A(Hk-p)|_inf < a_l
                     So there will be p s.t. g_l(r) > a_l
                     That's why I should check if r = 0
                     */
                    for ( kz = - mu - p2; kz <= mu + p2 ; kz++) {
                    //for ( kz = -2; kz <= 2 ; kz++) {
                        
                        for ( ky = - mu - p2 ; ky <= mu + p2 ; ky++) {
                        //for ( ky = -2; ky <= 2 ; ky++) {
                            
                            for ( kx = - mu - p2 ; kx <= mu + p2 ; kx++) {
                            //for ( kx = -2; kx <= 2 ; kx++) {

                                /*
                                 P LOOP
                                 */
                                double g = 0.0;
                                double part1 = 0.0;
                                double part2 = 0.0;
                                /*
                                int pzmin = -al/Az + (nz+kz)/(double)mlz ;
                                int pzmax =  al/Az + (nz+kz)/(double)mlz ;
                                int pymin = -al/Ay + (ny+ky)/(double)mly ;
                                int pymax =  al/Ay + (ny+ky)/(double)mly ;
                                int pxmin = -al/Ax + (nx+kx)/(double)mlx ;
                                int pxmax =  al/Ax + (nx+kx)/(double)mlx ; */
                                
                                int pzmin = -5; int pzmax = 5;
                                int pymin = -5; int pymax = 5;
                                int pxmin = -5; int pxmax = 5;
                                /* TODO: remove above before release */
                                for ( pz = pzmin ; pz <= pzmax ; pz++) {
                                    for ( py = pymin ; py <= pymax ; py++) {
                                        for ( px = pxmin ; px <= pxmax ; px++) {
                                            /* double r2 =
                                                    pow(Ax*(px-(nx+kx)/(double)mlx),2) +
                                                    pow(Ay*(py-(ny+ky)/(double)mly),2) +
                                                    pow(Az*(pz-(nz+kz)/(double)mlz),2); */
                                            double r2 = pow((nx + kx ) * hx - Ax * px,2) +
                                                        pow((ny + ky ) * hy - Ay * py,2) +
                                                        pow((nz + kz ) * hz - Az * pz,2) ;
                                            double r = sqrt(r2);

                                            /*
                                             * real-sum of gL* Equation 28 part 1
                                             */
                                            if (l == L - 1) {
                                                //double rho = r * al_1  ;
                                                double gL_value = g1(r,a,p);
                                                double gLR_value = gLR(r,beta);
                                                /*
                                                if ( r > alminus1 )
                                                    gL = 1.0/r ;
                                                else
                                                    gL = al_1 * 2.0 * TAUP4(2*rho) ;
                                                
                                                
                                                if (r < DBL_EPSILON)
                                                    gLR = ERFBETA0 ;
                                                else
                                                    gLR = erf(beta*r)/r ; */
                                                
                                                g     += gL_value - gLR_value ;
                                                part1 += gL_value  ;
                                                part2 += gLR_value ;
                                               
                                            } else {
                                                /* Equation 28 part 2 */
                                                if (r < alminus1) {
                                                    double rho = r * al_1  ;
                                                    g += al_1 * (2* TAUP4(2*rho) - TAUP4(rho));
                                                } else if ( r < al) {
                                                    double rho = r * al_1 ;
                                                    g += al_1 * (1./r - TAUP4(rho));
                                                }
                                            }
                                        }
                                    }
                                }
                                sl[l][nxoff]       += OMEGA[p2-2][abs(kx)] * OMEGA[p2-2][abs(ky)] * OMEGA[p2-2][abs(kz)] * g  ;
                                sl_part1[l][nxoff] += OMEGA[p2-2][abs(kx)] * OMEGA[p2-2][abs(ky)] * OMEGA[p2-2][abs(kz)] * part1  ;
                                sl_part2[l][nxoff] += OMEGA[p2-2][abs(kx)] * OMEGA[p2-2][abs(ky)] * OMEGA[p2-2][abs(kz)] * part2  ;
                                //sl[l][nxoff] += g;
                            }
                        }
                    }
                    //msm4g_linkedlist_addvalue(sllist, hash, sl[l][nxoff]);
                }
            }
        }
        
        printf("Stencil[%d]\n",l);
        for ( nz = 0 ; nz < sizez ; nz++ ) {
            int nzoff = sizex * sizey * nz ;
            for ( ny = 0 ; ny < sizey  ; ny++ ) {
                int nyoff = nzoff + sizex * ny ;
                for ( nx = 0 ; nx < sizex ; nx++ ) {
                    int nxoff = nyoff + nx;
                    printf("%2d x %2d x %2d %25.16f \n",nz-1,ny-1,nx-1,sl[l][nxoff]);
                }
            }
        }
        
        /* Print stencil */
        /*
        for ( nz = 0 ; nz < sizez ; nz++ ) {
            int nzoff = sizex * sizey * nz ;
            int nzt = nz + minz;
            for ( ny = 0 ; ny < sizey  ; ny++ ) {
                int nyoff = nzoff + sizex * ny ;
                int nyt = ny + miny ;
                for ( nx = 0 ; nx < sizex ; nx++ ) {
                    int nxoff = nyoff + nx;
                    int nxt = nx + minx;
                    printf("sl[%d %d %d] = %12.8f \n",nxt,nyt,nzt,sl[l][nxoff]);
                }
            }
        } */
        
        
        msm4g_linkedlist_destroy(sllist);
    }
    
    stats.cputime_stencill = (double)(clock()-begin)/CLOCKS_PER_SEC;

   
    /********************
     * Anterpolation
     * TODO: section: Anterpolation
     ********************/
    
    begin = clock();
    //printf("... Anterpolation\n");
    int ia = 0;
    int ib = Mx - 1;
    int ja = 0;
    int jb = My - 1 ;
    int ka = 0;
    int kb = Mz - 1;
    int ni = Mx;
    int nj = My;
    int nk = Mz;
    /*
    for (n = 0 ; n < N ; n++) {
        int i0 = floor((r[n*MAXDIM + DIMX] - Ox)/ hx) - (p2 - 1);
        int j0 = floor((r[n*MAXDIM + DIMY] - Oy)/ hy) - (p2 - 1);
        int k0 = floor((r[n*MAXDIM + DIMZ] - Oz)/ hz) - (p2 - 1);
        double rx = r[n*MAXDIM + DIMX] - hx*i0 - Ox;
        double ry = r[n*MAXDIM + DIMY] - hy*j0 - Oy;
        double rz = r[n*MAXDIM + DIMZ] - hz*k0 - Oz;


        for ( v = 0 ; v < p ; v++) {
            PHI[DIMX][v] = base->region[v](rx/hx - v);
            PHI[DIMY][v] = base->region[v](ry/hy - v);
            PHI[DIMZ][v] = base->region[v](rz/hz - v);
        }
        for ( k = 0 ; k < p ; k++) {
            int kk = k0 + k; 
            if      (kk < ka) do { kk += nk; } while (kk < ka);
            else if (kk > kb) do { kk -= nk; } while (kk > kb);
            int koff = ni*nj*kk;
            for ( j = 0 ; j < p ; j++) {
                int jj = j0 + j; 
                if      (jj < ja) do { jj += nj; } while (jj < ja);
                else if (jj > jb) do { jj -= nj; } while (jj > jb);

                int joff = koff + ni*jj;
                for ( i = 0 ; i < p ; i++) {
                    int ii =  i0 + i; 
                    if      (ii < ia) do { ii += ni; } while (ii < ia);
                    else if (ii > ib) do { ii -= ni; } while (ii > ib);
                    int ioff = joff + ii;
                    qgrd[0][ioff] +=  PHI[DIMX][i]*PHI[DIMY][j]*PHI[DIMZ][k] * q[n] ;
                }
            }
        }
    } */

    // Anterpolation corrected
    {
        for (int i = 0 ; i < N ; i++) {
            double rx = r[i * MAXDIM + DIMX] ;
            double ry = r[i * MAXDIM + DIMY] ;
            double rz = r[i * MAXDIM + DIMZ] ;
            double sx = ( rx - Ox ) / Ax ; // Position in [0 1]
            double sy = ( ry - Oy ) / Ay ;
            double sz = ( rz - Oz ) / Az ;
            for (int mz = - Mz/2 + 1 ; mz <= Mz/2 ; mz++) { 
                int mzoff = Mx * My * (mz + Mz/2 - 1) ;
                for (int my = - My/2 + 1 ; my <= My/2 ; my++) { 
                    int myoff = mzoff + Mx * (my + My/2 - 1);
                    for (int mx = - Mx/2 + 1 ; mx <= Mx/2 ; mx++) { 
                        int mxoff = myoff + (mx + Mx/2 - 1) ;
                        int pmax = p/2 ;
                        for (int px = - pmax ; px <= pmax ; px++) {
                            double phix = Bspline(p, Mx * (sx - px ) - mx + p2);
                            for (int py = - pmax ; py <= pmax ; py++) {
                                double phiy = Bspline(p, My * (sy - py ) - my + p2);
                                for (int pz = - pmax ; pz <= pmax ; pz++) {
                                    double phiz = Bspline(p, Mz * (sz - pz ) - mz + p2);
                                    qgrd[0][mxoff] += phix * phiy * phiz * q[i];
                                }
                            }
                        }
                    }
                }
            }
        }

    }
    
    stats.cputime_anterpolation = (double)(clock()-begin)/CLOCKS_PER_SEC;

    printf("grid charge at first level\n");
    double qgrdsum = 0.0;
    for ( k = 0 ; k < Mlz[0] ; k++ ) {
        for ( j = 0 ; j < Mly[0] ; j++ ) {
            for ( i = 0 ; i < Mlz[0] ; i++) {
                double charge = qgrd[0][ k * Mlx[0] * Mly[0] + j * Mlx[0] + i];
                qgrdsum += charge;
                printf("%25.16e\n",charge);
            }
        }
    }
    /* qgrd[0][0] =  1.0;
    qgrd[0][1] = -1.0;
    qgrd[0][2] = -1.0;
    qgrd[0][3] =  1.0;
    qgrd[0][4] = -1.0;
    qgrd[0][5] =  1.0;
    qgrd[0][6] =  1.0;
    qgrd[0][7] = -1.0;
    for (int i = 0 ; i < Mx * My * Mz ; i++)
        printf("%2d %25.16e\n",i,qgrd[0][i]); */
    
    
/*    qgrd[l] = (double *)calloc(Mlx[l]*Mly[l]*Mlz[l],sizeof(double)); */

    /********************
     * Restriction
     *
     *
     ********************/
    /* TODO: section: Restriction */
    begin = clock();
    for ( l = 0 ; l < L - 1  ; l++) {
        //printf("... Restriction : %d < %d\n",l,L-1);

        int mlx = Mlx[l] ;
        int mly = Mly[l] ;
        int mlz = Mlz[l] ;
        int mlxp1 = Mlx[l+1] ;
        int mlyp1 = Mly[l+1] ;
        int mlzp1 = Mlz[l+1] ;
        /* coarse level */
        int mxmin = - (p2 - 1);
        int mxmax = mlxp1 + (p2 - 1);
        int mymin = - (p2 - 1);;
        int mymax = mlyp1 + (p2 - 1);
        int mzmin = - (p2 - 1);
        int mzmax = mlzp1 + (p2 - 1);
        /* Lower and upper limits for grid charges */
        int mia = 0;
        int mib = mlxp1 - 1;
        int mja = 0;
        int mjb = mlyp1 - 1 ;
        int mka = 0;
        int mkb = mlzp1 - 1;
        /* Number of grid charges along each direction */
        int mi = mlxp1;
        int mj = mlyp1;
        int mk = mlzp1;

        /* fine level */
        /* Lower and upper limits for grid charges */
        int nia = 0;
        int nib = mlx - 1;
        int nja = 0;
        int njb = mly - 1 ;
        int nka = 0;
        int nkb = mlz - 1;
        /* Number of grid charges along each direction */
        int ni = mlx;
        int nj = mly;
        int nk = mlz;

        /* coarse level */
        for ( mz = mzmin ; mz <= mzmax ; mz++) {
            int mzz = mz;
            if      (mzz < mka) do { mzz += mk; } while (mzz < mka);
            else if (mzz > mkb) do { mzz -= mk; } while (mzz > mkb);
            int mzoff = mi * mj * mzz;
            for ( my = mymin; my <= mymax; my++) {
                int myy = my;
                if      (myy < mja) do { myy += mj; } while (myy < mja);
                else if (myy > mjb) do { myy -= mj; } while (myy > mjb);
                int myoff = mzoff + mi * myy;
                for ( mx = mxmin ; mx <= mxmax; mx++) {
                    int mxx = mx;
                    if      (mxx < mia) do { mxx += mi; } while (mxx < mia);
                    else if (mxx > mib) do { mxx -= mi; } while (mxx > mib);
                    int mxoff = myoff + mxx;

                    /* fine level */
                    for ( nz = 2*mz - p2 ; nz <= 2*mz + p2; nz++) {
                        int nzz = nz;
                        if      (nzz < nka) do { nzz += nk; } while (nzz < nka);
                        else if (nzz > nkb) do { nzz -= nk; } while (nzz > nkb);
                        int nzoff = ni * nj * nzz;

                        for ( ny = 2*my - p2 ; ny <= 2*my + p2 ; ny++) {
                            int nyy = ny;
                            if      (nyy < nja) do { nyy += nj; } while (nyy < nja);
                            else if (nyy > njb) do { nyy -= nj; } while (nyy > njb);
                            int nyoff = nzoff + ni * nyy;
                            for ( nx = 2*mx -p2 ; nx <= 2*mx + p2 ; nx++) {
                                int nxx = nx;
                                if      (nxx < nia) do { nxx += ni; } while (nxx < nia);
                                else if (nxx > nib) do { nxx -= ni; } while (nxx > nib);
                                int nxoff = nyoff + nxx;

                                qgrd[l+1][mxoff] += qgrd[l][nxoff] *
                                JN[p2-2][abs(nx-2*mx)] * JN[p2-2][abs(ny-2*my)] * JN[p2-2][abs(nz-2*mz)];
                            }
                        }
                    }
                }
            }
        }
    }
    stats.cputime_restriction = (double)(clock()-begin)/CLOCKS_PER_SEC;

    /*****************************
     * TODO: section: Grid-to-grid mapping
     *****************************/
    begin = clock();
    for (l = 0 ; l < L ; l++)
    {
        //printf("... Grid to grid mapping : %d < %d\n",l,L);
        /* loop: m
         Lower and upper limits for grid charges */
        int mia = 0;
        int mib = Mlx[l] - 1;
        int mja = 0;
        int mjb = Mly[l] - 1 ;
        int mka = 0;
        int mkb = Mlz[l] - 1;
        /* Number of grid charges along each direction */
        int mi = Mlx[l];
        int mj = Mly[l];
        //int mk = Mlz[l];

        /*  Size of the stencil */
        int si = slsizex[l];
        int sj = slsizey[l];
        int sk = slsizez[l];

        int sia = slminx[l];
        int sib = slmaxx[l];
        int sja = slminy[l];
        int sjb = slmaxy[l];
        int ska = slminz[l];
        int skb = slmaxz[l];
        int mx,my,mz,nx,ny,nz;
        for ( mz = mka ; mz <= mkb ; mz++) {
            int mzoff = mi * mj * mz;
            for ( my = mja ; my <= mjb ; my++) {
                int myoff = mzoff + mi * my;
                for ( mx = mia ; mx <= mib ; mx++) {
                    int mxoff = myoff + mx;

                    double esum = 0.0;
                    double esum_part1 = 0.0 ;
                    double esum_part2 = 0.0 ;
                    /* loop: n */
                    for ( nz = mka ; nz <= mkb ; nz++) {
                        int nzoff = mi * mj * nz;
                        for ( ny = mja ; ny <= mjb ; ny++) {
                            int nyoff = nzoff + mi * ny;
                            for ( nx = mia ; nx <= mib ; nx++) {
                                int nxoff = nyoff + nx;

                                
                                int sx = mx - nz ;
                                int sy = my - ny ;
                                int sz = mz - nz ;
                                if      (sx < sia) do {sx += si; } while (sx < sia);
                                else if (sx > sib) do {sx -= si; } while (sx > sib); 
                                if      (sy < sja) do {sy += sj; } while (sy < sja);
                                else if (sy > sjb) do {sy -= sj; } while (sy > sjb); 
                                if      (sz < ska) do {sz += sk; } while (sz < ska);
                                else if (sz > skb) do {sz -= sk; } while (sz > skb); 

                                int szoff = si * sj * (sz - ska) ;
                                int syoff = szoff + si * (sy - sja) ;
                                int sxoff = syoff + (sx - sia);

                                esum +=  sl[l][sxoff] * qgrd[l][nxoff] ;
                                esum_part1 +=  sl_part1[l][sxoff] * qgrd[l][nxoff] ;
                                esum_part2 +=  sl_part2[l][sxoff] * qgrd[l][nxoff] ;
                            }
                        }
                    }
                    edir[l][mxoff] = esum/pow(2,l);
                    edir_part1[l][mxoff] = esum_part1/pow(2,l);
                    edir_part2[l][mxoff] = esum_part2/pow(2,l);
                }
            }
        }
    }
    
    stats.cputime_grid_to_grid_mapping = (double)(clock()-begin)/CLOCKS_PER_SEC;

    
    /* for ( k = 0 ; k < Mlz[0] ; k++ ) {
        printf("edir slice: %d\n",k);
        for ( j = 0 ; j < Mly[0] ; j++ ) {
            for ( i = 0 ; i < Mlz[0] ; i++) {
                double pot = edir[0][ k * Mlx[0] * Mly[0] + j * Mlx[0] + i];
                printf("%9.4e ",pot);
            }
            printf("\n");
        }
    } */
    
    /* TODO: section: Top-level fourier calculation */
    begin=clock();
    double potential_long_fourier = 0.0;
    {
        //printf("... Top level fourier calculation\n");

        int mia = 0;
        int mib = MLx - 1;
        int mja = 0;
        int mjb = MLy - 1 ;
        int mka = 0;
        int mkb = MLz - 1;
        /* Number of grid charges along each direction */
        int mi = MLx;
        int mj = MLy;
        int mk = MLz;
        int si = sLnx;
        int sj = sLny;
        int sk = sLnz;
        int sia = - (MLx - 1);
        int sib =   (MLx - 1);
        int sja = - (MLy - 1);
        int sjb =   (MLy - 1);
        int ska = - (MLz - 1);
        int skb =   (MLz + 1);
 
        int mz,my,mx,nz,ny,nx;
        for ( mz = mka ; mz <= mkb ; mz++) {
            int mzoff = mi * mj * mz;
            for ( my = mja ; my <= mjb ; my++) {
                int myoff = mzoff + mi * my;
                for ( mx = mia ; mx <= mib ; mx++) {
                    int mxoff = myoff + mx;

                    for ( nz = mka ; nz <= mkb ; nz++) {
                        int nzoff = mi * mj * nz;
                        for ( ny = mja ; ny <= mjb ; ny++) {
                            int nyoff = nzoff + mi * ny;
                            for ( nx = mia ; nx <= mib ; nx++) {
                                int nxoff = nyoff + nx;

                                int sx = mx - nz ;
                                int sy = my - ny ;
                                int sz = mz - nz ;
                                if      (sx < sia) do {sx += si; } while (sx < sia);
                                else if (sx > sib) do {sx -= si; } while (sx > sib); 
                                if      (sy < sja) do {sy += sj; } while (sy < sja);
                                else if (sy > sjb) do {sy -= sj; } while (sy > sjb); 
                                if      (sz < ska) do {sz += sk; } while (sz < ska);
                                else if (sz > skb) do {sz -= sk; } while (sz > skb); 
 
                                int szoff = si * sj * (sz - ska) ;
                                int syoff = szoff + si * (sy - sja) ;
                                int sxoff = syoff + (sx - sia);

                                elng[L-1][mxoff] += qgrd[L-1][nxoff] * sL[sxoff] ;

                            }
                        }
                    }

                }
            }
        }

        for ( i = 0 ; i < mi * mj * mk ; i++)
        {
            potential_long_fourier += - 0.5 * elng[L-1][i] * qgrd[L-1][i];
        }
        stats.potential_long_fourier = potential_long_fourier;
        
        for ( i = 0 ; i < mi * mj * mk ; i++)
        {
            elng[L-1][i] = edir[L-1][i];
            elng_part1[L-1][i] = edir_part1[L-1][i];
            elng_part2[L-1][i] = edir_part2[L-1][i];
        }
    }
    stats.cputime_top_level = (double)(clock()-begin)/CLOCKS_PER_SEC;



    /* TODO: section: Prolongation */
    begin = clock();
    for ( l = L-2 ; l >= 0  ; l--)
    {
        //printf("... Prolongation : %d >= %d\n",l,0);

        /* coarse level */
        int mymin = - (p2 - 1);;
        int mymax = Mly[l+1] + (p2 - 1);
        int mxmin = - (p2 - 1);
        int mxmax = Mlx[l+1] + (p2 - 1);
        
        int mzmin = - (p2 - 1);
        int mzmax = Mlz[l+1] + (p2 - 1);
        /* Lower and upper limits for grid charges */
        int mia = 0;
        int mib = Mlx[l+1] - 1;
        int mja = 0;
        int mjb = Mly[l+1] - 1 ;
        int mka = 0;
        int mkb = Mlz[l+1] - 1;
        /* Number of grid charges along each direction */
        int mi = Mlx[l+1];
        int mj = Mly[l+1];
        int mk = Mlz[l+1];

        /* fine level
         Lower and upper limits for grid charges */
        int nia = 0;
        int nib = Mlx[l] - 1;
        int nja = 0;
        int njb = Mly[l] - 1 ;
        int nka = 0;
        int nkb = Mlz[l] - 1;
        /* Number of grid charges along each direction */
        int ni = Mlx[l];
        int nj = Mly[l];
        int nk = Mlz[l];

        /* coarse level */
        for ( mz = mzmin ; mz <= mzmax ; mz++)
        {
            int mzz = mz;
            if      (mzz < mka) do { mzz += mk; } while (mzz < mka);
            else if (mzz > mkb) do { mzz -= mk; } while (mzz > mkb);
            int mzoff = mi * mj * mzz;
            for ( my = mymin; my <= mymax; my++)
            {
                int myy = my;
                if      (myy < mja) do { myy += mj; } while (myy < mja);
                else if (myy > mjb) do { myy -= mj; } while (myy > mjb);
                int myoff = mzoff + mi * myy;
                for ( mx = mxmin ; mx <= mxmax; mx++)
                {
                    int mxx = mx;
                    if      (mxx < mia) do { mxx += mi; } while (mxx < mia);
                    else if (mxx > mib) do { mxx -= mi; } while (mxx > mib);
                    int mxoff = myoff + mxx;

                    /* fine level */
                    for ( nz = 2*mz - p2 ; nz <= 2*mz + p2; nz++)
                    {
                        int nzz = nz;
                        if      (nzz < nka) do { nzz += nk; } while (nzz < nka);
                        else if (nzz > nkb) do { nzz -= nk; } while (nzz > nkb);
                        int nzoff = ni * nj * nzz;

                        for ( ny = 2*my - p2 ; ny <= 2*my + p2 ; ny++)
                        {
                            int nyy = ny;
                            if      (nyy < nja) do { nyy += nj; } while (nyy < nja);
                            else if (nyy > njb) do { nyy -= nj; } while (nyy > njb);
                            int nyoff = nzoff + ni * nyy;
                            for ( nx = 2*mx -p2 ; nx <= 2*mx + p2 ; nx++)
                            {
                                int nxx = nx;
                                if      (nxx < nia) do { nxx += ni; } while (nxx < nia);
                                else if (nxx > nib) do { nxx -= ni; } while (nxx > nib);
                                int nxoff = nyoff + nxx;

                                elng[l][nxoff] += elng[l+1][mxoff] *
                                JN[p2-2][abs(nx-2*mx)] * JN[p2-2][abs(ny-2*my)] * JN[p2-2][abs(nz-2*mz)];

                            }
                        }
                    }
                }
            }
        }
        for ( i = 0 ; i < ni * nj * nk; i++)
            elng[l][i] += edir[l][i];
    }

    stats.cputime_prolongation = (double)(clock()-begin)/CLOCKS_PER_SEC;

    /* TODO: section: Interpolation */
    begin = clock();
    //printf("... Interpolation\n");

    for ( n = 0 ; n < N ; n++) /* for each particle */
    {
        int i0 = floor((r[n*MAXDIM + DIMX] - Ox)/ hx) - (p2 - 1);
        int j0 = floor((r[n*MAXDIM + DIMY] - Oy)/ hy) - (p2 - 1);
        int k0 = floor((r[n*MAXDIM + DIMZ] - Oz)/ hz) - (p2 - 1);
        double rx = r[n*MAXDIM + DIMX] - hx*i0 - Ox;
        double ry = r[n*MAXDIM + DIMY] - hy*j0 - Oy;
        double rz = r[n*MAXDIM + DIMZ] - hz*k0 - Oz;

        for ( v = 0 ; v < p ; v++)
        {
            PHI[DIMX][v] = base->region[v](rx/hx - v);
            PHI[DIMY][v] = base->region[v](ry/hy - v);
            PHI[DIMZ][v] = base->region[v](rz/hz - v);
            PHID[DIMX][v] = based->region[v](rx/hx - v)/hx;
            PHID[DIMY][v] = based->region[v](ry/hy - v)/hy;
            PHID[DIMZ][v] = based->region[v](rz/hz - v)/hz;
        }
        for ( k = 0 ; k < p ; k++)
        {
            int kk = k0 + k + p2 - 1 ;
            if      (kk < ka) do { kk += nk; } while (kk < ka);
            else if (kk > kb) do { kk -= nk; } while (kk > kb);
            int koff = ni*nj*kk;
            for ( j = 0 ; j < p ; j++)
            {
                int jj = j0 + j + p2 - 1 ;
                if      (jj < ja) do { jj += nj; } while (jj < ja);
                else if (jj > jb) do { jj -= nj; } while (jj > jb);
                int joff = koff + ni*jj;
                for ( i = 0 ; i < p ; i++)
                {
                    int ii =  i0 + i + p2 - 1 ;
                    if      (ii < ia) do { ii += ni; } while (ii < ia);
                    else if (ii > ib) do { ii -= ni; } while (ii > ib);
                    int ioff = joff + ii;
                    flong[n*MAXDIM + DIMX] += PHID[DIMX][i]*PHI [DIMY][j] * PHI [DIMZ][k]* elng[0][ioff];
                    flong[n*MAXDIM + DIMY] += PHI [DIMX][i]*PHID[DIMY][j] * PHI [DIMZ][k]* elng[0][ioff];
                    flong[n*MAXDIM + DIMZ] += PHI [DIMX][i]*PHI [DIMY][j] * PHID[DIMZ][k]* elng[0][ioff];
                }
            }
        }
        flong[n*MAXDIM + DIMX] *= q[n];
        flong[n*MAXDIM + DIMY] *= q[n];
        flong[n*MAXDIM + DIMZ] *= q[n];
    }

    stats.cputime_interpolation = (double)(clock()-begin)/CLOCKS_PER_SEC;

    begin=clock();
    //printf("... Short-range calculation \n");

    potential_short_real = 0.0;
    int nshortrangecal = 0;

    int useCellsInShortRange = 0;

    if (useCellsInShortRange)
    {
        /* Put the particles into cells */
        int nx= ceil((rmax[DIMX]-rmin[DIMX])/a);
        int ny= ceil((rmax[DIMY]-rmin[DIMY])/a);
        int nz= ceil((rmax[DIMZ]-rmin[DIMZ])/a);
        LinkedList **cells;
        cells = (LinkedList **)calloc(nx*ny*nz,sizeof(LinkedList *));
        for ( i = 0 ; i < nx*ny*nz ; i++)
            cells[i] = msm4g_linkedlist_new();

        for ( n = 0 ; n < N ; n++)
        {
           
            double rx = r[MAXDIM*n + DIMX];
            double ry = r[MAXDIM*n + DIMY];
            double rz = r[MAXDIM*n + DIMZ];

            int ix = floor((rx - rmin[DIMX])/a);
            int iy = floor((ry - rmin[DIMY])/a);
            int iz = floor((rz - rmin[DIMZ])/a);

            if (ix < 0 ) ix = 0;
            if (iy < 0 ) iy = 0;
            if (iz < 0 ) iz = 0;
            if (ix > nx - 1) ix = nx - 1;
            if (iy > ny - 1) iy = ny - 1;
            if (iz > nz - 1) iz = nz - 1;


            int icell = iz*ny*nx+iy*nx+ix;
            if (icell < 0 || icell > nx*ny*nz - 1)
            {
                fprintf(stderr,"error");
            }
            msm4g_linkedlist_add(cells[icell], n);
        }
        /* Computing cell interactions
         Cell neighbors */
        for ( k = 0 ; k < nz ; k++)
        {
            for ( j = 0 ; j < ny ; j++)
            {
                for ( i = 0 ; i < nx ; i++)
                {
                    int currcellindex = k * ny * nx + j * nx + i;
                    /* Within cell */
                    LinkedList *cell = cells[currcellindex];
                    LinkedListElement *curr = cell->head;
                    while (curr != NULL)
                    {
                        int i = curr->data;
                        LinkedListElement *next = curr->next;
                        while (next != NULL)
                        {
                            int j = next->data;
                            potential_short_real += q[i] * q[j] * short_range_periodic(i,j,r,q,p-1,a,fshort,Ax,Ay,Az);
                            nshortrangecal++;
                            next = next->next;
                        }
                        curr = curr->next;
                    }

                    /* Between cells */
                    for ( kk = k - 1 ; kk <= k + 1  ; kk++)
                    {
                        for ( jj = j - 1 ; jj <= j + 1 ; jj++)
                        {
                            for ( ii = i - 1 ; ii <= i + 1 ; ii++)
                            {
                                if ( kk <  0   || jj <  0  || ii < 0 ||
                                    kk >= nz  || jj >= ny || ii >= nx)
                                    continue;
                                int targetcellindex = kk * ny * nx + jj * nx + ii;
                                if (targetcellindex <= currcellindex) continue;


                                LinkedList *targetcell = cells[targetcellindex];

                                LinkedListElement *currentParticle = cell->head;
                                while (currentParticle!= NULL)
                                {
                                    int currentParticleIndex = currentParticle->data;

                                    LinkedListElement *targetParticle = targetcell->head;
                                    while (targetParticle != NULL)
                                    {
                                        int targetParticleIndex = targetParticle->data;

                                        potential_short_real +=  q[currentParticleIndex] * q[targetParticleIndex] *
                                           short_range_periodic(currentParticleIndex,targetParticleIndex,r,q,p-1,a,fshort,Ax,Ay,Az);
                                        nshortrangecal++;


                                        targetParticle = targetParticle->next;
                                    }

                                    currentParticle = currentParticle->next;
                                }



                            }
                        }
                    }



                }
            }
        }

        for ( i = 0 ; i < nx*ny*nz ; i++) {
            msm4g_linkedlist_destroy(cells[i]);
        }
        free(cells);

    } else /* Short-range without cells (for debugging) */
    {
        for ( i = 0 ; i < N ; i++)
        {
            for ( j = 0 ; j < N ; j++)
            {
                if ( i == j)
                    continue;
                potential_short_real += 0.5 * q[i] * q[j] * short_range_periodic(i,j,r,q,p-1,a,fshort,Ax,Ay,Az);
            }
        }


    }

    double potential_short_self = 0.0 ;
    /* Self-interaction term */
    for ( i = 0 ; i < N ; i++)
    {
        potential_short_self += 0.5 * q[i] * q[i] * short_range_periodic(i,i,r,q,p-1,a,fshort,Ax,Ay,Az);
    }

    /* Background correction */
    double qsum = 0.0, qsum2 = 0.0 ;
    for ( i = 0 ; i < N ; i++)
    {
        qsum  += q[i];
        qsum2 += q[i] * q[i] ;
    }
    double potential_short_csr =  - 0.5 * qsum * qsum * csr  ;

    double potential_short_total = potential_short_real + potential_short_self + potential_short_csr ;

    stats.potential_short_total = potential_short_total;
    stats.potential_short_real = potential_short_real;
    stats.potential_short_self = potential_short_self ;
    stats.potential_short_csr = potential_short_csr ;

    /* potential_long */
    begin=clock();
    double potential_long_real = 0.0;
    double potential_long_real_part1 = 0.0;
    double potential_long_real_part2 = 0.0;

    for ( i = 0 ; i < ni * nj * nk ; i++) {
        potential_long_real += 0.5 * elng[0][i] * qgrd[0][i] ;
        potential_long_real_part1 += 0.5 * elng_part1[0][i] * qgrd[0][i] ;
        potential_long_real_part2 += 0.5 * elng_part2[0][i] * qgrd[0][i] ;
    }
    
    printf("%-20s : %25.16e\n","    ulong_real_g1" ,potential_long_real_part1);
    printf("%-20s : %25.16e\n","    ulong_real_gLR" ,potential_long_real_part2);

    double potential_long_self = -0.5 * qsum2 * (1.0/a) * TAUP4(0.0);
    /* double potential_long_self = -0.5 * qsum2 * 2*beta/sqrt(PI); */

    /* Fourier is already incorporated into potential_long_real */
    /* TODO: decision to be made: should I include Fourier to the energy? */

    double potential_long_total = potential_long_fourier + potential_long_real + potential_long_self ;
    //double potential_long_total = potential_long_real + potential_long_self ;

    stats.cputime_potential_long = (double)(clock()-begin)/CLOCKS_PER_SEC;
    stats.potential_long_total = potential_long_total;
    stats.potential_long_real = potential_long_real;
    stats.potential_long_self = potential_long_self ;
    double potential_total = potential_long_total + potential_short_total;

    stats.potential_total = potential_total ;
    *uappr = potential_total ;

    //printf("... total potential %20.14e\n",potential_total);
    for ( i = 0 ; i < N ; i++)
    {
        fappr[i*MAXDIM + DIMX] = fshort[i*MAXDIM + DIMX] + flong[i*MAXDIM + DIMX];
        fappr[i*MAXDIM + DIMY] = fshort[i*MAXDIM + DIMY] + flong[i*MAXDIM + DIMY];
        fappr[i*MAXDIM + DIMZ] = fshort[i*MAXDIM + DIMZ] + flong[i*MAXDIM + DIMZ];
    }

    FILE *fp = fopen("msm.acclong","w");
    fprintf(fp,"%d\n",N);
    for (i=0;i<N;i++) fprintf(fp,"%20.14e\n",flong[3*i]);
    for (i=0;i<N;i++) fprintf(fp,"%20.14e\n",flong[3*i+1]);
    for (i=0;i<N;i++) fprintf(fp,"%20.14e\n",flong[3*i+2]);
    fclose(fp);
    
    fp = fopen("msm.accshort","w");
    fprintf(fp,"%d\n",N);
    for (i=0;i<N;i++) fprintf(fp,"%20.14e\n",fshort[3*i]);
    for (i=0;i<N;i++) fprintf(fp,"%20.14e\n",fshort[3*i+1]);
    for (i=0;i<N;i++) fprintf(fp,"%20.14e\n",fshort[3*i+2]);
    fclose(fp);
    
    stats.cputime_short_range = (double)(clock()-begin)/CLOCKS_PER_SEC;
    stats.short_range_interactions =  nshortrangecal;

    /* Sanity checks
     *
     */
    /* Sum of charges should stay same each level */
    for ( l = 0; l < L ; l++) {
        double sum = 0.0 ;
        for ( i = 0 ; i < Mlx[l]*Mly[l]*Mlz[l] ; i++)
          sum += qgrd[l][i];
        stats.netgridcharge[l] = sum ;
    }
    stats.netcharge = qsum ;

    stats.cputime_total = (double)(clock()-total_start)/CLOCKS_PER_SEC;

    /* Deallocations */
    free(sL);
    for ( l = 0; l < L ; l++) {
        free(qgrd[l]);
    }
    for ( l = 0; l < L ; l++) {
        free(edir[l]);
        free(sl[l]);
    }
    free(qgrd);
    free(edir);
    free(sl);
    free(flong);
    free(fshort);
    /* TODO: I should free it after print_stats
     free(cvecx);
     free(cvecy);
    free(cvecz);
    free(chi); */

    return stats;
}

int main(int argc,char *argv[])
{
    double *r, *q, *f, abar, Lx, Ly, Lz, potential;
    int n,p,mu,i;
    
    if (argc != 5 ) {
        printf("Usage: %s datafile p mu abar\n",
               argv[0]);
        return 1;
    }
    p = atoi(argv[2]);
    mu = atoi(argv[3]);
    abar=atof(argv[4]);
    
    data_read(argv[1], &q, &r, &n, &Lx, &Ly, &Lz);
    f = calloc(n*3, sizeof(double));
    potential = msmPeriodicRunner(r, q, f, n, p,mu,abar,Lx,Ly,Lz);
    printf("%20.15e\n", potential);

    FILE *fp = fopen("msm.pot","w");
    fprintf(fp,"%20.15e\n", potential);
    fclose(fp);
    fp = fopen("msm.acc","w");
    fprintf(fp,"%d\n",n);
    for (i=0;i<n;i++) fprintf(fp,"%20.14e\n",f[3*i]);
    for (i=0;i<n;i++) fprintf(fp,"%20.14e\n",f[3*i+1]);
    for (i=0;i<n;i++) fprintf(fp,"%20.14e\n",f[3*i+2]);
    fclose(fp);

    free(r);
    free(q);
    free(f);
    return 0;
}

double msmPeriodicRunner(double *r,double *q,double *force,int n,int p,int mu,double abar,double Lx,double Ly,double Lz)
{
    double potential;
    msm_stats statistics;
    statistics = msmperiodic(n, r, q, p,abar,mu,&potential,force, Lx, Ly, Lz);
    print_stats(statistics);
    return statistics.potential_total;
}
