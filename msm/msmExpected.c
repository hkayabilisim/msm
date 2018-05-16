#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <float.h>
#include "msmLibrary.h"
#define MYPI 3.141592653589793238462643

/*
static const double wprime[13] = {
    3.4641016151377539e+00, 
    -1.7320508075688767e+00,
    6.7949192431122685e-01,
    -2.3978297178318450e-01,
    7.9713982076653617e-02,
    -2.5502951436845282e-02,
    7.9437840692459776e-03,
    -2.4260315207973015e-03,
    7.2976833805943851e-04,
    -2.1633348486544432e-04,
    6.1410409127506203e-05, 
    -1.4537930745802521e-05,
    1.9569521248038610e-06}; */

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
double chi(double k,double beta,double volume) {
    return exp(-MYPI*MYPI*k*k/(beta*beta)) / (MYPI*k*k*volume) ;
}



double msmExpected(double *r, double *q, double *f, int N, int p, int mu,double abar, double Ax, double Ay, double Az,int pmax,int kmax) {
    int    L = 1;
    double utotal = 0.0;
    double volume = Ax * Ay * Az ;
    double htilde = pow(volume/N,1.0/3.0);
    double Mx = chooseM(htilde,Ax);
    double My = chooseM(htilde,Ay);
    double Mz = chooseM(htilde,Az);
    double hx = Ax / Mx ;
    double hy = Ay / My ;
    double hz = Az / Mz ;
    double a  = abar * htilde;
    double aL = pow(2,L)*a;
    double epsilon = 1e-12;
    double hL = pow(2,L-1)*hx;
    double beta = choosebetaNew(epsilon, aL,hL);
    
    /* (A) Real-Space */
    double ushort_real = 0.0;
    for (int i = 0 ; i < N ; i++) {
        for (int j = 0 ; j < N ; j++) {
            if (i == j) continue;
            for (int px = -pmax ; px <= pmax ; px++ ) {
                for (int py = -pmax; py <= pmax ; py++ ) {
                    for (int pz = -pmax ; pz < pmax ; pz++ ) {
                        double rx = r[i*3  ] - r[j*3  ] - Ax * px ;
                        double ry = r[i*3+1] - r[j*3+1] - Ay * py ;
                        double rz = r[i*3+2] - r[j*3+2] - Az * pz ;
                        double rlen2 = rx * rx + ry * ry + rz * rz ;
                        double rlen  = sqrt(rlen2);
                        ushort_real += 0.5 * q[i] * q[j] * g0(rlen,a,p) ;
                    }
                }
            }
        }
    }
    
    /* (B) Self-Interaction */
    double ushort_self = 0.0;
    for (int i = 0 ; i < N ; i++) {
        for (int px = -pmax ; px <= pmax ; px++ ) {
            for (int py = -pmax ; py <= pmax ; py++ ) {
                for (int pz = -pmax ; pz <= pmax ; pz++ ) {
                    if (px == 0 && py == 0 && pz == 0) continue;
                    double rx = Ax * px ;
                    double ry = Ay * py ;
                    double rz = Az * pz ;
                    double rlen2 = rx * rx + ry * ry + rz * rz ;
                    double rlen  = sqrt(rlen2);
                    ushort_self += 0.5 * q[i]*q[i] *g0(rlen,a,p);
                }
            }
        }
    }
    
    /* (C) CSR */
    double csr = MYPI / (beta * beta * volume );
    double qsum = 0.0;
    for (int i = 0 ; i < N ; i++) {
        qsum += q[i];
    }
    double ushort_csr = -0.5 * qsum * qsum * csr ;
    
    /* (D) Long-range Real-Space */
    double ulong_real = 0.0 ;
    double ulong_real_g1  = 0.0 ;
    double ulong_real_gLR = 0.0 ;
    double ulong_real_interpolated = 0.0;
    double ulong_real_interpolated_g1 = 0.0;
    double ulong_real_interpolated_gLR = 0.0;
    for (int i = 0 ; i < N ; i++ ) {
        double rix = r[i*3  ];
        double riy = r[i*3+1];
        double riz = r[i*3+2];
        for (int j = 0 ; j < N ; j++ ) {
            double rjx = r[j*3  ];
            double rjy = r[j*3+1];
            double rjz = r[j*3+2];

            double psum_g1 = 0.0;
            double psum_gLR = 0.0;
            double psum_interpolated_g1 = 0.0;
            double psum_interpolated_gLR = 0.0;
            
            for (int px = -pmax ; px <= pmax ; px++) {
                for (int py = -pmax ; py <= pmax ; py++) {
                    for (int pz = -pmax ; pz <= pmax ; pz++ ) {  
                        double rx = rix - rjx - Ax * px ;
                        double ry = riy - rjy - Ay * py ;
                        double rz = riz - rjz - Az * pz ;
                        double rlen2 = rx * rx + ry * ry + rz * rz ;
                        double rlen = sqrt(rlen2);

                        double g1_true = g1(rlen,a,p);
                        double gLR_true = gLR(rlen,beta);

                        psum_g1  += g1_true;
                        psum_gLR += gLR_true;

                        double g1_interpolated = 0.0;
                        double gLR_interpolated = 0.0;
/*
                        for (int mx = -Mx/2 + 1 ; mx <= Mx/2 ; mx++) {
                            for (int my = -My/2 + 1 ; my <= My/2 ; my++) {
                                for (int mz = -Mz/2 + 1 ; mz <= Mz/2 ; mz++) {
                                    double phim = 0;
                                    for (int ppx = -1 ; ppx <= 1 ; ppx++) {
                                        for (int ppy = -1 ; ppy <= 1 ; ppy++) {
                                            for (int ppz = -1 ; ppz <= 1 ; ppz++) {
                                               phim += Bspline(p,rix/hx - Mx*ppx - mx + p/2) *
                                                       Bspline(p,riy/hy - My*ppy - my + p/2) *
                                                       Bspline(p,riz/hz - Mz*ppz - mz + p/2) ;
                                            }
                                        }
                                    }
                                    if (fabs(phim) < DBL_EPSILON) {
                                        continue;
                                    }
                                    for (int nx = -Mx/2 + 1 ; nx <= Mx/2 ; nx++) {
                                        for (int ny = -My/2 + 1 ; ny < My/2 ; ny++) {
                                            for (int nz = -Mz/2 + 1; nz <= Mz/2 ; nz++) {
                                                double phin = 0;
                                                for (int ppx = -1 ; ppx <= 1 ; ppx++) {
                                                    for (int ppy = -1 ; ppy <= 1 ; ppy++) {
                                                        for (int ppz = -1 ; ppz <= 1 ; ppz++) {
                                                          phin += Bspline(p,rjx/hx - Mx*ppx - nx + p/2) *
                                                                  Bspline(p,rjy/hy - My*ppy - ny + p/2) *
                                                                  Bspline(p,rjz/hz - Mz*ppz - nz + p/2) ;
                                                        }
                                                    }
                                                }
                                                if (fabs(phin) < DBL_EPSILON) {
                                                    continue;
                                                }
                                                for (int kx = -mu-p/2 ; kx <= mu+p/2 ; kx++) {
                                                    for (int ky = -mu-p/2 ; ky <= mu+p/2; ky++) {
                                                        for (int kz = -mu-p/2 ; kz < mu+p/2 ; kz++) {
                                                            double w = wprime[abs(kx)]*wprime[abs(ky)]*wprime[abs(kz)];
                                                            
                                                            double rx = (mx - nx + kx ) * hx - Ax * px ;
                                                            double ry = (my - ny + ky ) * hy - Ay * py ;
                                                            double rz = (mz - nz + kz ) * hz - Az * pz ;
                                                            double rlen2 = rx * rx + ry * ry + rz * rz ;
                                                            double rlen = sqrt(rlen2);
                                                            
                                                            g1_interpolated  += phim * phin * w * g1(rlen,a,p);
                                                            gLR_interpolated += phim * phin * w * gLR(rlen,beta);
                                                            
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        } */
 
                        psum_interpolated_g1  += g1_interpolated  ;
                        psum_interpolated_gLR += gLR_interpolated ;
                        //printf("%2d %2d %10.5f %25.16e %25.16e %25.16e %25.16e\n",i,j,rlen,g1_true,g1_interpolated,gLR_true,gLR_interpolated);
                    }
                }
            }
            ulong_real += 0.5 * q[i] * q[j] * (psum_g1 - psum_gLR) ;
            ulong_real_g1  += 0.5 * q[i] * q[j] * (psum_g1) ;
            ulong_real_gLR += 0.5 * q[i] * q[j] * (psum_gLR) ;

            ulong_real_interpolated     += 0.5 * q[i] * q[j] * (psum_interpolated_g1 - psum_interpolated_gLR) ;
            ulong_real_interpolated_g1  += 0.5 * q[i] * q[j] * (psum_interpolated_g1) ;
            ulong_real_interpolated_gLR += 0.5 * q[i] * q[j] * (psum_interpolated_gLR) ;
        }
    }
    
    /* (D) Interpolated Long-range Real-Space */

    /*
    double ulong_real_interpolated = 0.0;
    double ulong_real_interpolated_g1 = 0.0;
    double ulong_real_interpolated_gLR = 0.0;
    for (int i = 0 ; i < N ; i++ ) {
        double rix = r[i*3  ];
        double riy = r[i*3+1];
        double riz = r[i*3+2];
        for (int j = 0 ; j < N ; j++ ) {
            double rjx = r[j*3  ];
            double rjy = r[j*3+1];
            double rjz = r[j*3+2];
   
            double psum = 0.0;
            double psum_g1 = 0.0;
            double psum_gLR = 0.0;
            double psum_interpolated = 0.0;
            double psum_interpolated_g1 = 0.0;
            double psum_interpolated_gLR = 0.0;

            double rijnorm = sqrt(pow(rix-rjx,2)+pow(riy-rjy,2)+pow(riz-rjz,2));
            
            for (int px = -pmax ; px <= pmax ; px++) {
                for (int py = -pmax ; py <= pmax ; py++) {
                    for (int pz = -pmax ; pz <= pmax ; pz++ ) {  
                  
            
                        double rx = rix - rjx - Ax * px ;
                        double ry = riy - rjy - Ay * py ;
                        double rz = riz - rjz - Az * pz ;
                        double rlen2 = rx * rx + ry * ry + rz * rz ;
                        double rlen = sqrt(rlen2);
                        psum_g1  += g1(rlen,a,p);
                        psum_gLR += gLR(rlen,beta);
                        
                        
                        
                        for (int mx = -Mx/2 + 1 ; mx <= Mx/2 ; mx++) {
                            for (int my = -My/2 + 1 ; my <= My/2 ; my++) {
                                for (int mz = -Mz/2 + 1 ; mz <= Mz/2 ; mz++) {
                                    double phim = 0;
                                    for (int ppx = -1 ; ppx <= 1 ; ppx++) {
                                        for (int ppy = -1 ; ppy <= 1 ; ppy++) {
                                            for (int ppz = -1 ; ppz <= 1 ; ppz++) {
                                               phim += Bspline(p,rix/hx - Mx*ppx - mx + p/2) *
                                                       Bspline(p,riy/hy - My*ppy - my + p/2) *
                                                       Bspline(p,riz/hz - Mz*ppz - mz + p/2) ;
                                            }
                                        }
                                    }
                                    if (fabs(phim) < DBL_EPSILON) {
                                        continue;
                                    }
                                    for (int nx = -Mx/2 + 1 ; nx <= Mx/2 ; nx++) {
                                        for (int ny = -My/2 + 1 ; ny < My/2 ; ny++) {
                                            for (int nz = -Mz/2 + 1; nz <= Mz/2 ; nz++) {
                                                double phin = 0;
                                                for (int ppx = -1 ; ppx <= 1 ; ppx++) {
                                                    for (int ppy = -1 ; ppy <= 1 ; ppy++) {
                                                        for (int ppz = -1 ; ppz <= 1 ; ppz++) {
                                                          phin += Bspline(p,rjx/hx - Mx*ppx - nx + p/2) *
                                                                  Bspline(p,rjy/hy - My*ppy - ny + p/2) *
                                                                  Bspline(p,rjz/hz - Mz*ppz - nz + p/2) ;
                                                        }
                                                    }
                                                }
                                                if (fabs(phin) < DBL_EPSILON) {
                                                    continue;
                                                }
                                                for (int kx = -mu-p/2 ; kx <= mu+p/2 ; kx++) {
                                                    for (int ky = -mu-p/2 ; ky <= mu+p/2; ky++) {
                                                        for (int kz = -mu-p/2 ; kz < mu+p/2 ; kz++) {
                                                            double w = wprime[abs(kx)]*wprime[abs(ky)]*wprime[abs(kz)];
                                                            
                                                            double rx = (mx - nx + kx ) * hx - Ax * px ;
                                                            double ry = (my - ny + ky ) * hy - Ay * py ;
                                                            double rz = (mz - nz + kz ) * hz - Az * pz ;
                                                            double rlen2 = rx * rx + ry * ry + rz * rz ;
                                                            double rlen = sqrt(rlen2);
                                                            psum_interpolated += phim * phin * w * g1star(rlen,a,p,beta)  ;
                                                            
                                                            psum_interpolated_g1  += phim * phin * w * g1(rlen,a,p);
                                                            psum_interpolated_gLR += phim * phin * w * gLR(rlen,beta);
                                                            
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        
                    }
                }
            }
            printf("%2d %2d %10.5e %25.16e %25.16e %25.16e %25.16e %25.16e %25.16e\n",i,j,
                        rijnorm,psum_g1 , psum_interpolated_g1  ,fabs(psum_g1 -psum_interpolated_g1 )/fabs(psum_g1),
                        psum_gLR, psum_interpolated_gLR ,fabs(psum_gLR-psum_interpolated_gLR)/fabs(psum_gLR));
            
            ulong_real += 0.5 * q[i] * q[j] * (psum_g1 - psum_gLR) ;
            ulong_real_g1  += 0.5 * q[i] * q[j] * (psum_g1) ;
            ulong_real_gLR += 0.5 * q[i] * q[j] * (psum_gLR) ;
            ulong_real_interpolated    += 0.5 * q[i] * q[j] * (psum_interpolated_g1 - psum_interpolated_gLR);
            ulong_real_interpolated_g1  += 0.5 * q[i] * q[j] * (psum_interpolated_g1);
            ulong_real_interpolated_gLR += 0.5 * q[i] * q[j] * (psum_interpolated_gLR);
        }
    } */
    
    /* (E) Long-range Reciprocal-Space */
    double ulong_four = 0.0 ;
    double cos_sum = 0.0;
    double sin_sum = 0.0 ;
    for (int i = 0 ; i < N ; i++) {
        for (int j = 0 ; j < N ; j++ ) {
            for (int kx = -kmax ; kx <= kmax ; kx++) {
                for (int ky = -kmax ; ky <= kmax ; ky++) {
                    for (int kz = -kmax ; kz <= kmax ; kz++) {
                        if (kx == 0 && ky == 0 && kz == 0 ) continue;
                        double klen2 = (kx/Ax) * (kx/Ax) + (ky/Ay)*(ky/Ay) + (kz/Az)*(kz/Az) ;
                        double klen = sqrt(klen2);
                        double dotprodx = (kx/Ax) * (r[i*3  ]-r[j*3  ]) ;
                        double dotprody = (ky/Ay) * (r[i*3+1]-r[j*3+1]) ;
                        double dotprodz = (kz/Az) * (r[i*3+2]-r[j*3+2]) ;
                        double dotprod = dotprodx + dotprody + dotprodz ;
                        cos_sum += 0.5 * q[i] * q[j] * chi(klen,beta,volume) * cos(2*MYPI*dotprod);
                        sin_sum += 0.5 * q[i] * q[j] * chi(klen,beta,volume) * sin(2*MYPI*dotprod);
                    }
                }
            }
        }
    }
    ulong_four = cos_sum ;
    /* (F) Long-Range self-interaction */
    double ulong_self = 0;
    for (int i = 0 ; i < N ; i++) {
        ulong_self += -0.5 * q[i] * q[i] * g1(0.0,a,p) ;
    }
    
    utotal = ushort_real + ushort_self + ushort_csr + ulong_real + ulong_four + ulong_self;
    
    
    /* Printing report */
    printf("%-20s : %25.16e\n","h  (estimated)" ,htilde);
    printf("%-20s : %25.16e\n","hx (calculated)",hx);
    printf("%-20s : %25.16e\n","hy (calculated)",hy);
    printf("%-20s : %25.16e\n","hz (calculated)",hz);
    printf("%-20s : %25.16e\n","beta",beta);
    printf("%-20s : %25.16e\n","volume",volume);
    printf("%-20s : %25.16e\n","csr",csr);
    printf("%-20s : %25.16e\n","qsum",qsum);
    printf("%-20s : %25.16e\n","cos_sum",cos_sum);
    printf("%-20s : %25.16e\n","sin_sum",sin_sum);
    printf("%-20s : %25.16e\n","(A) ushort_real",ushort_real);
    printf("%-20s : %25.16e\n","(B) ushort_self",ushort_self);
    printf("%-20s : %25.16e\n","(C) ushort_csr" ,ushort_csr);
    printf("%-20s : %25.16e\n","(D) ulong_real" ,ulong_real);
    printf("%-20s : %25.16e\n","    ulong_real_g1"  ,ulong_real_g1);
    printf("%-20s : %25.16e\n","    ulong_real_gLR" ,ulong_real_gLR);
    printf("%-20s : %25.16e\n","    ulong_real*" ,ulong_real_interpolated);
    printf("%-20s : %25.16e\n","    ulong_real*_g1" ,ulong_real_interpolated_g1);
    printf("%-20s : %25.16e\n","    ulong_real*_gLR" ,ulong_real_interpolated_gLR);
    printf("%-20s : %25.16e\n","(E) ulong_four" ,ulong_four);
    printf("%-20s : %25.16e\n","(F) ulong_self" ,ulong_self);
    printf("%-20s : %25.16e\n","    utotal"     ,utotal);
    
    return utotal;
}

int main(int argc,char *argv[]) {
    double *r, *q, *f, abar, Lx, Ly, Lz;
    int n,p,mu,pmax,kmax;
    
    if (argc != 7 ) {
        printf("Usage: %s datafile p mu abar pmax kmax\n",
               argv[0]);
        return 1;
    }
    p = atoi(argv[2]);
    mu = atoi(argv[3]);
    abar=atof(argv[4]);
    pmax=atof(argv[5]);
    kmax=atof(argv[6]);
    
    data_read(argv[1], &q, &r, &n, &Lx, &Ly, &Lz);
    f = calloc(n*3, sizeof(double));
    msmExpected(r, q, f, n, p,mu,abar,Lx,Ly,Lz,pmax,kmax);
    //for (int i=0; i <= 40 ; i++)
    //    printf("%8.2f %8.4f\n",i/10.0,Bspline(4,i/10.0));
    free(r);
    free(q);
    free(f);
    return 0;
}


