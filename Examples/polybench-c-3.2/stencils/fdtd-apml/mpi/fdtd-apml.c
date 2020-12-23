#include "mpi.h" 
#include <string.h> 
#define TAG 123
/*
Copyright (C) 1991-2014 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it andor
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, see
   <http:www.gnu.org/licenses/>. 
*/
/*
This header is separate from features.h so that the compiler can
   include it implicitly at the start of every compilation.  It must
   not itself include <features.h> or any other header that includes
   <features.h> because the implicit include comes before any feature
   test macros that may be defined in a source file before it first
   explicitly includes a system header.  GCC knows the name of this
   header in order to preinclude it. 
*/
/*
glibc's intent is to support the IEC 559 math functionality, real
   and complex.  If the GCC (4.9 and later) predefined macros
   specifying compiler intent are available, use them to determine
   whether the overall intent is to support these features; otherwise,
   presume an older compiler has intent to support these features and
   define these macros by default. 
*/
/*
wchar_t uses ISOIEC 10646 (2nd ed., published 2011-03-15) /
   Unicode 6.0. 
*/
/* We do not support C11 <threads.h>.  */
/*

 * fdtd-apml.c: This file is part of the PolyBench/C 3.2 test suite.
 *
 *
 * Contact: Louis-Noel Pouchet <pouchet@cse.ohio-state.edu>
 * Web address: http:polybench.sourceforge.net

*/
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
/* Include polybench common header. */
#include <polybench.h>
/* Include benchmark-specific header. */
/* Default data type is double, default size is 256x256x256. */
#include "fdtd-apml.h"
/* Array initialization. */
static void init_array(int cz, int cxm, int cym, double * mui, double * ch, double Ax[((32+1)+0)][((32+1)+0)], double Ry[((32+1)+0)][((32+1)+0)], double Ex[((32+1)+0)][((32+1)+0)][((32+1)+0)], double Ey[((32+1)+0)][((32+1)+0)][((32+1)+0)], double Hz[((32+1)+0)][((32+1)+0)][((32+1)+0)], double czm[((32+1)+0)], double czp[((32+1)+0)], double cxmh[((32+1)+0)], double cxph[((32+1)+0)], double cymh[((32+1)+0)], double cyph[((32+1)+0)])
{
	int i, j, k;
	( * mui)=2341;
	( * ch)=42;
	for (i=0; i<=cz; i ++ )
	{
		czm[i]=((((double)i)+1)/cxm);
		czp[i]=((((double)i)+2)/cxm);
	}
	for (i=0; i<=cxm; i ++ )
	{
		cxmh[i]=((((double)i)+3)/cxm);
		cxph[i]=((((double)i)+4)/cxm);
	}
	for (i=0; i<=cym; i ++ )
	{
		cymh[i]=((((double)i)+5)/cxm);
		cyph[i]=((((double)i)+6)/cxm);
	}
	for (i=0; i<=cz; i ++ )
	{
		for (j=0; j<=cym; j ++ )
		{
			Ry[i][j]=(((((double)i)*(j+1))+10)/cym);
			Ax[i][j]=(((((double)i)*(j+2))+11)/cym);
			for (k=0; k<=cxm; k ++ )
			{
				Ex[i][j][k]=((((((double)i)*(j+3))+k)+1)/cxm);
				Ey[i][j][k]=((((((double)i)*(j+4))+k)+2)/cym);
				Hz[i][j][k]=((((((double)i)*(j+5))+k)+3)/cz);
			}
		}
	}
	return ;
}

/*
DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output.
*/
static void print_array(int cz, int cxm, int cym, double Bza[((32+1)+0)][((32+1)+0)][((32+1)+0)], double Ex[((32+1)+0)][((32+1)+0)][((32+1)+0)], double Ey[((32+1)+0)][((32+1)+0)][((32+1)+0)], double Hz[((32+1)+0)][((32+1)+0)][((32+1)+0)])
{
	int i, j, k;
	for (i=0; i<=cz; i ++ )
	{
		for (j=0; j<=cym; j ++ )
		{
			for (k=0; k<=cxm; k ++ )
			{
				fprintf(stderr, "%0.2lf ", Bza[i][j][k]);
				fprintf(stderr, "%0.2lf ", Ex[i][j][k]);
				fprintf(stderr, "%0.2lf ", Ey[i][j][k]);
				fprintf(stderr, "%0.2lf ", Hz[i][j][k]);
				if ((((i*cxm)+j)%20)==0)
				{
					fprintf(stderr, "\n");
				}
			}
		}
	}
	fprintf(stderr, "\n");
	return ;
}

/*
Main computational kernel. The whole function will be timed,
   including the call and return.
*/
/*
rcv buffer for accessed data - clf
*/
double rcvBuffclf[((32+1)+0)][((32+1)+0)];
/*
rcv buffer for accessed data - tmp
*/
double rcvBufftmp[((32+1)+0)][((32+1)+0)];
/*
rcv buffer for accessed data - Hz
*/
double rcvBuffHz[((32+1)+0)][((32+1)+0)][((32+1)+0)];
/*
rcv buffer for accessed data - Bza
*/
double rcvBuffBza[((32+1)+0)][((32+1)+0)][((32+1)+0)];
/*
rcv buffer for accessed data - Hz
*/
double rcvBuffHz[((32+1)+0)][((32+1)+0)][((32+1)+0)];
/*
rcv buffer for accessed data - Bza
*/
double rcvBuffBza[((32+1)+0)][((32+1)+0)][((32+1)+0)];
/*
rcv buffer for accessed data - Hz
*/
double rcvBuffHz[((32+1)+0)][((32+1)+0)][((32+1)+0)];
/*
rcv buffer for accessed data - Bza
*/
double rcvBuffBza[((32+1)+0)][((32+1)+0)][((32+1)+0)];
/*
rcv buffer for accessed data - Hz
*/
double rcvBuffHz[((32+1)+0)][((32+1)+0)][((32+1)+0)];
/*
rcv buffer for accessed data - Bza
*/
double rcvBuffBza[((32+1)+0)][((32+1)+0)][((32+1)+0)];
int main(int argc, char * * argv)
{
	/* Retrieve problem size. */
	int cz = 32;
	int cym = 32;
	int cxm = 32;
	/* Variable declarationallocation. */
	double mui;
	double ch;
	double Ax[((32+1)+0)][((32+1)+0)];
	double Ry[((32+1)+0)][((32+1)+0)];
	double clf[((32+1)+0)][((32+1)+0)];
	double tmp[((32+1)+0)][((32+1)+0)];
	double Bza[((32+1)+0)][((32+1)+0)][((32+1)+0)];
	double Ex[((32+1)+0)][((32+1)+0)][((32+1)+0)];
	double Ey[((32+1)+0)][((32+1)+0)][((32+1)+0)];
	double Hz[((32+1)+0)][((32+1)+0)][((32+1)+0)];
	double czm[((32+1)+0)];
	double czp[((32+1)+0)];
	double cxmh[((32+1)+0)];
	double cxph[((32+1)+0)];
	double cymh[((32+1)+0)];
	double cyph[((32+1)+0)];
	/* Initialize array(s). */
	int iz, iy, ix;
	int _ret_val_0;
	/*
	MPI variables declaration  begins
	*/
	MPI_Status status;
	int rank;
	int size;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	/*
	MPI variables declaration ends
	*/
	init_array(cz, cxm, cym,  & mui,  & ch, Ax, Ry, Ex, Ey, Hz, czm, czp, cxmh, cxph, cymh, cyph);
	/* Start timer. */
	;
	/* Run kernel. */
	#pragma scop 
	#pragma MPI  parallel 0
	#pragma MPI  accessed_data(clf[iz][iy],tmp[iz][iy],Hz[iz][iy][ix],Bza[iz][iy][ix],Hz[iz][iy][32],Bza[iz][iy][32],Hz[iz][32][ix],Bza[iz][32][ix],Hz[iz][32][32],Bza[iz][32][32])
	for (iz=0; iz<32; iz ++ )
	/*
	Computation partitioning and mapping code begins 
	*/
	{
		int Ic;
		int Rc;
		int Lc;
		int Uc;
		int Lcp;
		int Ucp;
		Lc=0;
		Uc=32;
		Ic=(Uc-Lc);
		Rc=(Ic%size);
		if (rank==0)
		{
			Lcp=Lc;
			Ucp=((Lcp+(Ic/size))+Rc);
		}
		else
		{
			Lcp=((Lc+(rank*(Ic/size)))+Rc);
			Ucp=(Lcp+(Ic/size));
		}
		if ((iz>=Lcp)&&(iz<Ucp))
		{
			for (iy=0; iy<32; iy ++ )
			{
				for (ix=0; ix<32; ix ++ )
				{
					clf[iz][iy]=(((Ex[iz][iy][ix]-Ex[iz][iy+1][ix])+Ey[iz][iy][ix+1])-Ey[iz][iy][ix]);
					tmp[iz][iy]=(((cymh[iy]/cyph[iy])*Bza[iz][iy][ix])-((ch/cyph[iy])*clf[iz][iy]));
					Hz[iz][iy][ix]=((((cxmh[ix]/cxph[ix])*Hz[iz][iy][ix])+(((mui*czp[iz])/cxph[ix])*tmp[iz][iy]))-(((mui*czm[iz])/cxph[ix])*Bza[iz][iy][ix]));
					Bza[iz][iy][ix]=tmp[iz][iy];
				}
				clf[iz][iy]=(((Ex[iz][iy][32]-Ex[iz][iy+1][32])+Ry[iz][iy])-Ey[iz][iy][32]);
				tmp[iz][iy]=(((cymh[iy]/cyph[iy])*Bza[iz][iy][32])-((ch/cyph[iy])*clf[iz][iy]));
				Hz[iz][iy][32]=((((cxmh[32]/cxph[32])*Hz[iz][iy][32])+(((mui*czp[iz])/cxph[32])*tmp[iz][iy]))-(((mui*czm[iz])/cxph[32])*Bza[iz][iy][32]));
				Bza[iz][iy][32]=tmp[iz][iy];
				for (ix=0; ix<32; ix ++ )
				{
					clf[iz][iy]=(((Ex[iz][32][ix]-Ax[iz][ix])+Ey[iz][32][ix+1])-Ey[iz][32][ix]);
					tmp[iz][iy]=(((cymh[32]/cyph[iy])*Bza[iz][iy][ix])-((ch/cyph[iy])*clf[iz][iy]));
					Hz[iz][32][ix]=((((cxmh[ix]/cxph[ix])*Hz[iz][32][ix])+(((mui*czp[iz])/cxph[ix])*tmp[iz][iy]))-(((mui*czm[iz])/cxph[ix])*Bza[iz][32][ix]));
					Bza[iz][32][ix]=tmp[iz][iy];
				}
				clf[iz][iy]=(((Ex[iz][32][32]-Ax[iz][32])+Ry[iz][32])-Ey[iz][32][32]);
				tmp[iz][iy]=(((cymh[32]/cyph[32])*Bza[iz][32][32])-((ch/cyph[32])*clf[iz][iy]));
				Hz[iz][32][32]=((((cxmh[32]/cxph[32])*Hz[iz][32][32])+(((mui*czp[iz])/cxph[32])*tmp[iz][iy]))-(((mui*czm[iz])/cxph[32])*Bza[iz][32][32]));
				Bza[iz][32][32]=tmp[iz][iy];
			}
		}
	}
	/*
	COMMUNICATION CODE STARTS
	*/
	if (size>1)
	{
		int Ld;
		int Ud;
		int Id;
		int Rd;
		int CD_L_p;
		int CD_U_p;
		int CD_S_p;
		int PSM;
		/*
		For accessed data - clf- gathering
		*/
		{
			Ld=0;
			Ud=32;
			Id=(Ud-Ld);
			Rd=(Id%size);
			PSM=33;
			if (rank==0)
			{
				CD_L_p=Ld;
				CD_U_p=((CD_L_p+(Id/size))+Rd);
			}
			else
			{
				CD_L_p=((Ld+(rank*(Id/size)))+Rd);
				CD_U_p=(CD_L_p+(Id/size));
			}
			if ((rank==0)&&(Rd!=0))
			{
				CD_L_p=(CD_L_p+Rd);
				memcpy(rcvBuffclf, clf, sizeof(clf));
			}
			CD_S_p=((CD_U_p-CD_L_p)*PSM);
			MPI_Gather(&clf[CD_L_p][0], CD_S_p, MPI_DOUBLE, &rcvBuffclf[CD_L_p][0], CD_S_p, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			if (rank==0)
			{
				memcpy(clf, rcvBuffclf, sizeof(rcvBuffclf));
			}
		}
		/*
		For accessed data - tmp- gathering
		*/
		{
			Ld=0;
			Ud=32;
			Id=(Ud-Ld);
			Rd=(Id%size);
			PSM=33;
			if (rank==0)
			{
				CD_L_p=Ld;
				CD_U_p=((CD_L_p+(Id/size))+Rd);
			}
			else
			{
				CD_L_p=((Ld+(rank*(Id/size)))+Rd);
				CD_U_p=(CD_L_p+(Id/size));
			}
			if ((rank==0)&&(Rd!=0))
			{
				CD_L_p=(CD_L_p+Rd);
				memcpy(rcvBufftmp, tmp, sizeof(tmp));
			}
			CD_S_p=((CD_U_p-CD_L_p)*PSM);
			MPI_Gather(&tmp[CD_L_p][0], CD_S_p, MPI_DOUBLE, &rcvBufftmp[CD_L_p][0], CD_S_p, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			if (rank==0)
			{
				memcpy(tmp, rcvBufftmp, sizeof(rcvBufftmp));
			}
		}
		/*
		For accessed data - Hz- gathering
		*/
		{
			Ld=0;
			Ud=32;
			Id=(Ud-Ld);
			Rd=(Id%size);
			PSM=1089;
			if (rank==0)
			{
				CD_L_p=Ld;
				CD_U_p=((CD_L_p+(Id/size))+Rd);
			}
			else
			{
				CD_L_p=((Ld+(rank*(Id/size)))+Rd);
				CD_U_p=(CD_L_p+(Id/size));
			}
			if ((rank==0)&&(Rd!=0))
			{
				CD_L_p=(CD_L_p+Rd);
				memcpy(rcvBuffHz, Hz, sizeof(Hz));
			}
			CD_S_p=((CD_U_p-CD_L_p)*PSM);
			MPI_Gather(&Hz[CD_L_p][0][0], CD_S_p, MPI_DOUBLE, &rcvBuffHz[CD_L_p][0][0], CD_S_p, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			if (rank==0)
			{
				memcpy(Hz, rcvBuffHz, sizeof(rcvBuffHz));
			}
		}
		/*
		For accessed data - Bza- gathering
		*/
		{
			Ld=0;
			Ud=32;
			Id=(Ud-Ld);
			Rd=(Id%size);
			PSM=1089;
			if (rank==0)
			{
				CD_L_p=Ld;
				CD_U_p=((CD_L_p+(Id/size))+Rd);
			}
			else
			{
				CD_L_p=((Ld+(rank*(Id/size)))+Rd);
				CD_U_p=(CD_L_p+(Id/size));
			}
			if ((rank==0)&&(Rd!=0))
			{
				CD_L_p=(CD_L_p+Rd);
				memcpy(rcvBuffBza, Bza, sizeof(Bza));
			}
			CD_S_p=((CD_U_p-CD_L_p)*PSM);
			MPI_Gather(&Bza[CD_L_p][0][0], CD_S_p, MPI_DOUBLE, &rcvBuffBza[CD_L_p][0][0], CD_S_p, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			if (rank==0)
			{
				memcpy(Bza, rcvBuffBza, sizeof(rcvBuffBza));
			}
		}
		/*
		For accessed data - Hz- gathering
		*/
		{
			Ld=0;
			Ud=32;
			Id=(Ud-Ld);
			Rd=(Id%size);
			PSM=1089;
			if (rank==0)
			{
				CD_L_p=Ld;
				CD_U_p=((CD_L_p+(Id/size))+Rd);
			}
			else
			{
				CD_L_p=((Ld+(rank*(Id/size)))+Rd);
				CD_U_p=(CD_L_p+(Id/size));
			}
			if ((rank==0)&&(Rd!=0))
			{
				CD_L_p=(CD_L_p+Rd);
				memcpy(rcvBuffHz, Hz, sizeof(Hz));
			}
			CD_S_p=((CD_U_p-CD_L_p)*PSM);
			MPI_Gather(&Hz[CD_L_p][0][32], CD_S_p, MPI_DOUBLE, &rcvBuffHz[CD_L_p][0][32], CD_S_p, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			if (rank==0)
			{
				memcpy(Hz, rcvBuffHz, sizeof(rcvBuffHz));
			}
		}
		/*
		For accessed data - Bza- gathering
		*/
		{
			Ld=0;
			Ud=32;
			Id=(Ud-Ld);
			Rd=(Id%size);
			PSM=1089;
			if (rank==0)
			{
				CD_L_p=Ld;
				CD_U_p=((CD_L_p+(Id/size))+Rd);
			}
			else
			{
				CD_L_p=((Ld+(rank*(Id/size)))+Rd);
				CD_U_p=(CD_L_p+(Id/size));
			}
			if ((rank==0)&&(Rd!=0))
			{
				CD_L_p=(CD_L_p+Rd);
				memcpy(rcvBuffBza, Bza, sizeof(Bza));
			}
			CD_S_p=((CD_U_p-CD_L_p)*PSM);
			MPI_Gather(&Bza[CD_L_p][0][32], CD_S_p, MPI_DOUBLE, &rcvBuffBza[CD_L_p][0][32], CD_S_p, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			if (rank==0)
			{
				memcpy(Bza, rcvBuffBza, sizeof(rcvBuffBza));
			}
		}
		/*
		For accessed data - Hz- gathering
		*/
		{
			Ld=0;
			Ud=32;
			Id=(Ud-Ld);
			Rd=(Id%size);
			PSM=1089;
			if (rank==0)
			{
				CD_L_p=Ld;
				CD_U_p=((CD_L_p+(Id/size))+Rd);
			}
			else
			{
				CD_L_p=((Ld+(rank*(Id/size)))+Rd);
				CD_U_p=(CD_L_p+(Id/size));
			}
			if ((rank==0)&&(Rd!=0))
			{
				CD_L_p=(CD_L_p+Rd);
				memcpy(rcvBuffHz, Hz, sizeof(Hz));
			}
			CD_S_p=((CD_U_p-CD_L_p)*PSM);
			MPI_Gather(&Hz[CD_L_p][32][0], CD_S_p, MPI_DOUBLE, &rcvBuffHz[CD_L_p][32][0], CD_S_p, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			if (rank==0)
			{
				memcpy(Hz, rcvBuffHz, sizeof(rcvBuffHz));
			}
		}
		/*
		For accessed data - Bza- gathering
		*/
		{
			Ld=0;
			Ud=32;
			Id=(Ud-Ld);
			Rd=(Id%size);
			PSM=1089;
			if (rank==0)
			{
				CD_L_p=Ld;
				CD_U_p=((CD_L_p+(Id/size))+Rd);
			}
			else
			{
				CD_L_p=((Ld+(rank*(Id/size)))+Rd);
				CD_U_p=(CD_L_p+(Id/size));
			}
			if ((rank==0)&&(Rd!=0))
			{
				CD_L_p=(CD_L_p+Rd);
				memcpy(rcvBuffBza, Bza, sizeof(Bza));
			}
			CD_S_p=((CD_U_p-CD_L_p)*PSM);
			MPI_Gather(&Bza[CD_L_p][32][0], CD_S_p, MPI_DOUBLE, &rcvBuffBza[CD_L_p][32][0], CD_S_p, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			if (rank==0)
			{
				memcpy(Bza, rcvBuffBza, sizeof(rcvBuffBza));
			}
		}
		/*
		For accessed data - Hz- gathering
		*/
		{
			Ld=0;
			Ud=32;
			Id=(Ud-Ld);
			Rd=(Id%size);
			PSM=1089;
			if (rank==0)
			{
				CD_L_p=Ld;
				CD_U_p=((CD_L_p+(Id/size))+Rd);
			}
			else
			{
				CD_L_p=((Ld+(rank*(Id/size)))+Rd);
				CD_U_p=(CD_L_p+(Id/size));
			}
			if ((rank==0)&&(Rd!=0))
			{
				CD_L_p=(CD_L_p+Rd);
				memcpy(rcvBuffHz, Hz, sizeof(Hz));
			}
			CD_S_p=((CD_U_p-CD_L_p)*PSM);
			MPI_Gather(&Hz[CD_L_p][32][32], CD_S_p, MPI_DOUBLE, &rcvBuffHz[CD_L_p][32][32], CD_S_p, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			if (rank==0)
			{
				memcpy(Hz, rcvBuffHz, sizeof(rcvBuffHz));
			}
		}
		/*
		For accessed data - Bza- gathering
		*/
		{
			Ld=0;
			Ud=32;
			Id=(Ud-Ld);
			Rd=(Id%size);
			PSM=1089;
			if (rank==0)
			{
				CD_L_p=Ld;
				CD_U_p=((CD_L_p+(Id/size))+Rd);
			}
			else
			{
				CD_L_p=((Ld+(rank*(Id/size)))+Rd);
				CD_U_p=(CD_L_p+(Id/size));
			}
			if ((rank==0)&&(Rd!=0))
			{
				CD_L_p=(CD_L_p+Rd);
				memcpy(rcvBuffBza, Bza, sizeof(Bza));
			}
			CD_S_p=((CD_U_p-CD_L_p)*PSM);
			MPI_Gather(&Bza[CD_L_p][32][32], CD_S_p, MPI_DOUBLE, &rcvBuffBza[CD_L_p][32][32], CD_S_p, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			if (rank==0)
			{
				memcpy(Bza, rcvBuffBza, sizeof(rcvBuffBza));
			}
		}
	}
	#pragma endscop 
	/* Stop and print timer. */
	;
	;
	/*
	Prevent dead-code elimination. All live-out data must be printed
	     by the function call in argument.
	*/
	if (rank==0)
	{
		print_array(cz, cxm, cym, Bza, Ex, Ey, Hz);
	}
	/* Be clean. */
	;
	;
	;
	;
	;
	;
	;
	;
	;
	;
	;
	;
	;
	;
	_ret_val_0=0;
	MPI_Finalize();
	return _ret_val_0;
}

/*
=>========================MPI Analysis===================
=> Benchmark =
=> TOTAL NO. OF PARALLEL LOOPS FOUND = 1
=> INNER PARALLEL LOOPS FOUND =0
=> OUTER PARALLEL LOOPS FOUND = 1
=>==========================================================
*/
