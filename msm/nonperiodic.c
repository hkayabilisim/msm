#include "msmLibrary.h"
#include <stdio.h>
int main(int argc,char *argv[])
{
	double *r,*q,Lx,Ly,Lz,u,d,d2;
	int N,i,j;
	data_read(argv[1],&q,&r,&N,&Lx,&Ly,&Lz);

    for (i = 0 ; i < N ; i++)
        q[i] = 1.0;

	u = 0.0 ;
	for (i=0;i<N;i++) {
		for (j = 0 ; j< N;j++) {
			if (i == j) continue;
			d2 = pow(r[3*i]  - r[3*j],2) + 
			     pow(r[3*i+1]- r[3*j+1],2) +
			     pow(r[3*i+2]- r[3*j+2],2) ;
			d = sqrt(d2);
			u += 0.5 * q[i]*q[j]/d;
	}
	}

	printf("u is %15.8e\n",u);
	free(q);
	free(r);
	return 0;
}

