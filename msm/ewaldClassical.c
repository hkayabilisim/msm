#include "ewaldClassical.h"
#include "msmLibrary.h"
#include <math.h>

double ewaldClassical(double *r,double *q,int N,
                      double Lx,double Ly,double Lz,
                      double alpha,double cutoffDirect,double  cutoffFourier)
{
    double ereal = 0.0;
    double efour = 0.0;
    double eself = 0.0;
    double etota = 0.0;
    int kmax ; 
    /* double PI = 3.14159265358979323846; */
    double V = Lx*Ly*Lz;
    
    int i,j,l,kx,ky,kz;
    double rx,ry,rz,rjx,rjy,rjz,qj,rlx,rly,rlz,ql,rjlx,rjly,rjlz,rjl2,rjl,e;
    double kxL,kyL,kzL,k2,creal,cimag,dotprod,csquare,factor,q2sum,qsum,ecorr;
    
    for (j = 0 ; j < N ; j++)
    {
        rjx = r[j*3]   ; 
        rjy = r[j*3+1] ; 
        rjz = r[j*3+2] ; 
        qj  = q[j]     ; 
        
        for (l = 0 ; l < N ; l++)
        {
            if (j == l) continue;
            rlx = r[l*3]   ; 
            rly = r[l*3+1] ; 
            rlz = r[l*3+2] ; 
            ql  = q[l]     ; 
            
            rjlx = rjx - rlx;
            rjly = rjy - rly;
            rjlz = rjz - rlz;
            
            /* Minimum image convention */
            /* rjlx = rjlx - Lx * floor(rjlx / Lx ) ;
            rjly = rjly - Ly * floor(rjly / Ly ) ;
            rjlz = rjlz - Lz * floor(rjlz / Lz ) ; */
            
            rjl2 = rjlx*rjlx + rjly*rjly + rjlz*rjlz;
            rjl  = sqrt(rjl2);
            
            if (rjl <= cutoffDirect)
            {
                e = 0.5 * qj * ql * erfc(alpha * rjl) / rjl ;
                ereal += e;
            }
        }
    }
    kmax = (int) ceil(cutoffFourier);
    for (kx = -kmax ; kx <= kmax ; kx++) {
        for (ky = -kmax ; ky <= kmax ; ky++) {
            for (kz = -kmax ; kz <= kmax ; kz++) {
                if (kx == 0 && ky == 0 && kz == 0)
                    continue;
                
                kxL = kx/Lx;
                kyL = ky/Ly;
                kzL = kz/Lz;
                
                k2 = kx * kx + ky * ky + kz * kz ;
                
                if (k2 >= cutoffFourier*cutoffFourier)
                    continue;
                
                creal = 0.0;
                cimag = 0.0;
                for (j = 0 ; j < N ; j++)
                {
                    rx = r[j*3    ];
                    ry = r[j*3 + 1];
                    rz = r[j*3 + 2];
                    
                    dotprod = kxL * rx  + kyL * ry + kzL * rz ;
                    creal += q[j] * cos(2 * PI * dotprod );
                    cimag += q[j] * sin(2 * PI * dotprod );
                }
                csquare = creal * creal + cimag * cimag ;
                
                k2 = (kxL * kxL + kyL * kyL + kzL * kzL) ;
                factor = exp(-PI * PI * k2 / (alpha * alpha) ) / k2 ;
                efour +=  0.5 * (1.0 / (PI * V)) * factor * csquare ;
            }
        }
    }
    
    q2sum = 0.0;
    for (i = 0 ; i < N ; i++)
    {
        q2sum += q[i] * q[i];;
    }
    eself = - 0.5 * q2sum * (2 * alpha  / sqrt(PI));
    
    qsum = 0.0;
    for (i = 0 ; i < N ; i++)
        qsum += q[i];
    ecorr = - 0.5 * qsum * qsum * (PI / (alpha*alpha*V));
    etota = ereal + efour + eself + ecorr;
    fprintf(stderr,"\e[1;34m Classical \e[0m\n");
    fprintf(stderr,"\e[1;34m alpha = %12.4e,L=%12.4e  \e[0m\n",alpha,Lx);
    fprintf(stderr,"\e[1;34m %12.4e %12.4e %12.4e %12.4e = %12.4e \e[0m\n",
            ereal,efour,eself,ecorr,etota); 
    return  etota;
}

void ewaldClassicalRunner(char *dataFile, double alpha, double cutoffDirect, double cutoffFourier)
{
    double *r, *q, Lx, Ly, Lz, potential;
    int n;
    
    data_read(dataFile, &q, &r, &n, &Lx, &Ly, &Lz);
    
    potential = ewaldClassical(r, q, n, Lx,Ly,Lz,alpha,cutoffDirect,cutoffFourier);
    printf("%20.15e\n", potential);
    
    free(r);
    free(q);
}

int main(int argc, char *argv[]) {
    if (argc != 5 ) {
        printf("Usage: %s datafile alpha cutoffDirect cutoffFourier\n",
               argv[0]);
        return 1;
    }
    ewaldClassicalRunner(argv[1],atof(argv[2]),atof(argv[3]),atof(argv[4]));
    return 0;
}
