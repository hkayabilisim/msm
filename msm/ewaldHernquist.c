#include "ewaldHernquist.h"
#include <math.h>

double ewaldHernquist(double *r, double *q, double *f, int n, double L,double alpha, double cutoffDirect, double cutoffFourier) {
    double x, y, z, q2sum, U;
    int i, j,nxmin,nxmax,nymin,nymax,nzmin,nzmax,nx,ny,nz,hmax,hx,hy,hz,h;
    double ereal,efour,eself,ecorr,etota;
    double rx,ry,rz,rlen,h2,dotprod,rlen2,rlen3;
    /* double PI = 3.14159265358979323846; */
    double H1 = 0,H2=0,H3=0,H4=0,H5=0,H6=0,H7=0;
    
    ereal = 0.0;
    efour = 0.0;
    eself = 0.0;
    ecorr = 0.0;
    
    U =   2.8372975 / L;
    
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i == j) continue;
            x = r[3 * i]     - r[3 * j];
            y = r[3 * i + 1] - r[3 * j + 1];
            z = r[3 * i + 2] - r[3 * j + 2];
            /* fprintf(stderr,"|r_%d - r_%d|\n",i,j); */
            
            nxmin = (int) floor(-cutoffDirect + x / L);
            nymin = (int) floor(-cutoffDirect + y / L);
            nzmin = (int) floor(-cutoffDirect + z / L);
            nxmax = (int)  ceil( cutoffDirect + x / L);
            nymax = (int)  ceil( cutoffDirect + y / L);
            nzmax = (int)  ceil( cutoffDirect + z / L);
            nxmin = -4; nxmax = 4; nymin = -4 ; nymax = 4; nzmin = -4; nzmax = 4; 
            for ( nx = nxmin; nx <= nxmax; nx++ ) {
                for ( ny = nymin; ny <= nymax; ny++) {
                    for ( nz = nzmin; nz <= nzmax; nz++) {
                        rx = x - nx * L;
                        ry = y - ny * L;
                        rz = z - nz * L;
                        rlen = sqrt(rx * rx + ry * ry + rz * rz);
                        /* if (rlen <= cutoffDirect) { */
                            ereal += 0.5 * q[i] * q[j] * erfc(alpha * rlen) / rlen;
                            H1 += 0.5 * q[i] * q[j] * erfc(alpha * rlen) / rlen;
                        /* } */
                    }
                }
            }
            
            /* Fourier sum */
            /* hmax = (int) ceil(cutoffFourier); */
            hmax = 4;
            for (hx = -hmax; hx <= hmax; hx++) {
                for (hy = -hmax; hy <= hmax; hy++) {
                    for (hz = -hmax; hz <= hmax; hz++) {
                        if (hx == 0 && hy == 0 && hz == 0)
                            continue;
                        h2 = hx * hx + hy * hy + hz * hz;
                        /* if (h2 >= cutoffFourier*cutoffFourier)
                            continue; */
                        dotprod = hx * x + hy * y + hz * z;
                        efour +=  0.5 * q[i] * q[j] * (1.0 / (L * PI * h2))
                        * exp(-PI * PI * h2 / (alpha * alpha * L * L))
                        * cos(2.0 * PI * dotprod / L);
                        H2 += 0.5 * q[i] * q[j] * (1.0 / (L * PI * h2))
                        * exp(-PI * PI * h2 / (alpha * alpha * L * L))
                        * cos(2.0 * PI * dotprod / L);
                        
                    }
                }
            }
            /* Correction term*/
            ecorr +=  0.5 * q[i] * q[j] * PI / (alpha * alpha * L * L * L);
        }
    }

    q2sum = 0.0;
    for (i = 0; i < n; i++)
        q2sum += q[i] * q[i];
    
    hmax = 10;
    for (hx = -hmax; hx <= hmax; hx++) {
        for (hy = -hmax; hy <= hmax; hy++) {
           for (hz = -hmax; hz <= hmax; hz++) {
                if (hx == 0 && hy == 0 && hz == 0)
                    continue;
                h2 = hx * hx + hy * hy + hz * hz;
                h = sqrt(h2);
                H3 +=  0.5 * q2sum *  (1.0 / (L * PI * h2))
                * exp(-PI * PI * h2 / (alpha * alpha * L * L));
                H4 += 0.5 * q2sum * erfc(alpha * L * h)/ (L*h);
            }
        }
    }
    
  
    H5 = 0.5 * q2sum * PI / (alpha*alpha*L*L*L) ;
    for (i = 0 ; i < n ; i++) {
        for (j = 0 ; j < n ; j++) {
            if (i == j) continue;
            H6 += 0.5 * q[i] * q[j] * PI /(alpha*alpha*L*L*L) ;
        }
    }

    H7 = 0.5 * q2sum * 2 * alpha / sqrt(PI) ;
    
    eself =  0.5 * q2sum * U ;
    
    etota = -ereal - efour + ecorr + eself ; 


    /* Now it is time to calculate forces 
     * Eq 2.14a and b in Hernquist 1991 */
    for (i=0;i<n;i++) {
	for (j=0;j<n;j++) { 
	    if (j == i) continue;
            x = r[3 * i]     - r[3 * j];
            y = r[3 * i + 1] - r[3 * j + 1];
            z = r[3 * i + 2] - r[3 * j + 2];
            /* fprintf(stderr,"|r_%d - r_%d|\n",i,j); */
            
            nxmin = (int) floor(-cutoffDirect + x / L);
            nymin = (int) floor(-cutoffDirect + y / L);
            nzmin = (int) floor(-cutoffDirect + z / L);
            nxmax = (int)  ceil( cutoffDirect + x / L);
            nymax = (int)  ceil( cutoffDirect + y / L);
            nzmax = (int)  ceil( cutoffDirect + z / L);
            nxmin = -4; nxmax = 4; nymin = -4 ; nymax = 4; nzmin = -4; nzmax = 4; 
            for ( nx = nxmin; nx <= nxmax; nx++ ) {
                for ( ny = nymin; ny <= nymax; ny++) {
                    for ( nz = nzmin; nz <= nzmax; nz++) {
                        rx = x - nx * L;
                        ry = y - ny * L;
                        rz = z - nz * L;
                        rlen2 = rx * rx + ry * ry + rz * rz;
                        rlen  = sqrt(rlen2);
			rlen3 = rlen * rlen2 ;

			double factor = erfc(alpha*rlen) + 
				 (2*alpha/sqrt(PI))  * rlen * exp(-alpha*alpha*rlen2);
			f[3*i]    -=  q[j]*rx * factor / rlen3 ;
			f[3*i +1] -=  q[j]*ry * factor / rlen3 ;
			f[3*i +2] -=  q[j]*rz * factor / rlen3 ;
			

		    }
		}
	    }



    for (hx = -hmax; hx <= hmax; hx++) {
        for (hy = -hmax; hy <= hmax; hy++) {
           for (hz = -hmax; hz <= hmax; hz++) {
                if (hx == 0 && hy == 0 && hz == 0)
                    continue;
                h2 = hx * hx + hy * hy + hz * hz;
		double dotprod = hx * x + hy * y + hz * z;

		f[3*i + 0 ] -= q[j]*hx * (2.0/(L*L)) * exp(-PI*PI*h2/(alpha*alpha*L*L)) 
			 * sin(2*PI*dotprod/L) / h2;
		f[3*i + 1 ] -= q[j]*hy * (2.0/(L*L)) * exp(-PI*PI*h2/(alpha*alpha*L*L)) 
			 * sin(2*PI*dotprod/L) / h2;
		f[3*i + 2 ] -= q[j]*hz * (2.0/(L*L)) * exp(-PI*PI*h2/(alpha*alpha*L*L)) 
			 * sin(2*PI*dotprod/L) / h2;
            }
        }
    }
    
 
 		
	}
    }

    /*
    fprintf(stderr,"\e[1;34m Hernquist \e[0m\n");
    fprintf(stderr,"\e[1;34m alpha = %12.4e, L=%12.4e  \e[0m\n",alpha,L);
    fprintf(stderr,"\e[1;34m eself = %12.4e\e[0m\n",eself);
    fprintf(stderr,"\e[1;34m H1 = %12.4e  \e[0m\n",H1);
    fprintf(stderr,"\e[1;34m H2 = %12.4e  \e[0m\n",H2);
    fprintf(stderr,"\e[1;34m H3 = %12.4e  \e[0m\n",H3);
    fprintf(stderr,"\e[1;34m H4 = %12.4e  \e[0m\n",H4);
    fprintf(stderr,"\e[1;34m H5 = %12.4e  \e[0m\n",H5);
    fprintf(stderr,"\e[1;34m H6 = %12.4e  \e[0m\n",H6);
    fprintf(stderr,"\e[1;34m H7 = %12.4e  \e[0m\n",H7);
    fprintf(stderr,"\e[1;34m Hsum = %12.4e  \e[0m\n",H1 + H2 + H3 + H4 - H5 - H6 - H7);
    fprintf(stderr,"\e[1;34m L*(H3 + H4 - H5 -H7) = %12.4e  \e[0m\n",(H3 + H4 - H5 - H7)*L*2/q2sum );

    fprintf(stderr,"\e[1;34m %12.4e %12.4e %12.4e %12.4e = %12.4e \e[0m\n",
            ereal,efour,eself,ecorr,etota); */
    printf("%-20s : %25.16e\n","    utotal"     ,etota);
    return  etota;
}

void ewaldHernQuistRunner(char *dataFile, double alpha, double cutoffDirect, double cutoffFourier)
{
    double *r, *q, *f,Lx, Ly, Lz, potential;
    int n;
    
    data_read(dataFile, &q, &r, &n, &Lx, &Ly, &Lz);
    f = (double *)calloc(3*n,sizeof(double));
    potential = ewaldHernquist(r, q, f, n, Lx,alpha,cutoffDirect,cutoffFourier);
    /*printf("%20.15e\n", potential);

    FILE *fp = fopen("hernquist.pot","w");
    fprintf(fp,"%20.15e\n", potential);
    fclose(fp);
    
    fp = fopen("hernquist.acc","w");
    fprintf(fp,"%d\n",n);
    for (i=0;i<n;i++) fprintf(fp,"%20.16e\n",f[3*i]);
    for (i=0;i<n;i++) fprintf(fp,"%20.16e\n",f[3*i+1]);
    for (i=0;i<n;i++) fprintf(fp,"%20.16e\n",f[3*i+2]);
    fclose(fp); */
    free(r);
    free(f);
    free(q);
}

int main(int argc, char *argv[]) {
    if (argc != 5 ) {
        printf("Usage: %s datafile alpha cutoffDirect cutoffFourier\n",
               argv[0]);
        return 1;
    }
    ewaldHernQuistRunner(argv[1],atof(argv[2]),atof(argv[3]),atof(argv[4]));
    return 0;
}
