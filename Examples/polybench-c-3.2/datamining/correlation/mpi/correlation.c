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

 * correlation.c: This file is part of the PolyBench/C 3.2 test suite.
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
/* Default data type is double, default size is 4000. */
#include "correlation.h"
/* Array initialization. */
static void init_array(int m, int n, double * float_n, double data[(32+0)][(32+0)])
{
	int i, j;
	( * float_n)=1.2;
	for (i=0; i<m; i ++ )
	{
		for (j=0; j<n; j ++ )
		{
			data[i][j]=((((double)i)*j)/32);
		}
	}
	return ;
}

/*
DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output.
*/
static void print_array(int m, double symmat[(32+0)][(32+0)])
{
	int i, j;
	for (i=0; i<m; i ++ )
	{
		for (j=0; j<m; j ++ )
		{
			fprintf(stderr, "%0.2lf ", symmat[i][j]);
			if ((((i*m)+j)%20)==0)
			{
				fprintf(stderr, "\n");
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
rcv buffer for accessed data - symmat
*/
double rcvBuffsymmat[(32+0)][(32+0)];
/*
rcv buffer for accessed data - symmat
*/
double rcvBuffsymmat[(32+0)][(32+0)];
/*
rcv buffer for accessed data - symmat
*/
double rcvBuffsymmat[(32+0)][(32+0)];
int main(int argc, char * * argv)
{
	/* Retrieve problem size. */
	int n = 32;
	int m = 32;
	/* Variable declarationallocation. */
	double float_n;
	double data[(32+0)][(32+0)];
	double symmat[(32+0)][(32+0)];
	double mean[(32+0)];
	double stddev[(32+0)];
	/* Initialize array(s). */
	int i, j, j1, j2;
	double eps = 0.1F;
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
	init_array(m, n,  & float_n, data);
	/* Start timer. */
	;
	/* Run kernel. */
	#pragma scop 
	/* Determine mean of column vectors of input data matrix */
	#pragma MPI  parallel 0
	#pragma MPI  accessed_data(mean[j])
	for (j=0; j<32; j ++ )
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
		if ((j>=Lcp)&&(j<Ucp))
		{
			mean[j]=0.0;
			for (i=0; i<32; i ++ )
			{
				mean[j]+=data[i][j];
			}
			mean[j]/=float_n;
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
		int p;
		for (p=0; p<size; p ++ )
		{
			/*
			For accessed data - mean - all-to-all broadcasting
			*/
			{
				Ld=0;
				Ud=32;
				Id=(Ud-Ld);
				Rd=(Id%size);
				if (p==0)
				{
					CD_L_p=Ld;
					CD_U_p=((CD_L_p+(Id/size))+Rd);
				}
				else
				{
					CD_L_p=((Ld+(p*(Id/size)))+Rd);
					CD_U_p=(CD_L_p+(Id/size));
				}
				PSM=1;
				CD_S_p=((CD_U_p-CD_L_p)*PSM);
				MPI_Bcast(&mean[CD_L_p], CD_S_p, MPI_DOUBLE, p, MPI_COMM_WORLD);
			}
		}
	}
	/* Determine standard deviations of column vectors of data matrix. */
	#pragma MPI  parallel 0
	#pragma MPI  accessed_data(stddev[j])
	for (j=0; j<32; j ++ )
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
		if ((j>=Lcp)&&(j<Ucp))
		{
			stddev[j]=0.0;
			for (i=0; i<32; i ++ )
			{
				stddev[j]+=((data[i][j]-mean[j])*(data[i][j]-mean[j]));
			}
			stddev[j]/=float_n;
			stddev[j]=sqrt(stddev[j]);
			/*
			The following in an inelegant but usual way to handle
				 near-zero std. dev. values, which below would cause a zero-
				 divide.
			*/
			stddev[j]=((stddev[j]<=eps) ? 1.0 : stddev[j]);
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
		int p;
		for (p=0; p<size; p ++ )
		{
			/*
			For accessed data - stddev - all-to-all broadcasting
			*/
			{
				Ld=0;
				Ud=32;
				Id=(Ud-Ld);
				Rd=(Id%size);
				if (p==0)
				{
					CD_L_p=Ld;
					CD_U_p=((CD_L_p+(Id/size))+Rd);
				}
				else
				{
					CD_L_p=((Ld+(p*(Id/size)))+Rd);
					CD_U_p=(CD_L_p+(Id/size));
				}
				PSM=1;
				CD_S_p=((CD_U_p-CD_L_p)*PSM);
				MPI_Bcast(&stddev[CD_L_p], CD_S_p, MPI_DOUBLE, p, MPI_COMM_WORLD);
			}
		}
	}
	/* Center and reduce the column vectors. */
	#pragma MPI  parallel 0
	#pragma MPI  accessed_data(data[i][j])
	for (i=0; i<32; i ++ )
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
		if ((i>=Lcp)&&(i<Ucp))
		{
			for (j=0; j<32; j ++ )
			{
				data[i][j]-=mean[j];
				data[i][j]/=(sqrt(float_n)*stddev[j]);
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
		int p;
		for (p=0; p<size; p ++ )
		{
			/*
			For accessed data - data - all-to-all broadcasting
			*/
			{
				Ld=0;
				Ud=32;
				Id=(Ud-Ld);
				Rd=(Id%size);
				if (p==0)
				{
					CD_L_p=Ld;
					CD_U_p=((CD_L_p+(Id/size))+Rd);
				}
				else
				{
					CD_L_p=((Ld+(p*(Id/size)))+Rd);
					CD_U_p=(CD_L_p+(Id/size));
				}
				PSM=32;
				CD_S_p=((CD_U_p-CD_L_p)*PSM);
				MPI_Bcast(&data[CD_L_p][0], CD_S_p, MPI_DOUBLE, p, MPI_COMM_WORLD);
			}
		}
	}
	/* Calculate the m m correlation matrix. */
	#pragma MPI  parallel 0
	#pragma MPI  accessed_data(symmat[j1][j1],symmat[j1][j2],symmat[j2][j1])
	for (j1=0; j1<31; j1 ++ )
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
		Uc=31;
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
		if ((j1>=Lcp)&&(j1<Ucp))
		{
			symmat[j1][j1]=1.0;
			for (j2=(j1+1); j2<32; j2 ++ )
			{
				symmat[j1][j2]=0.0;
				for (i=0; i<32; i ++ )
				{
					symmat[j1][j2]+=(data[i][j1]*data[i][j2]);
				}
				symmat[j2][j1]=symmat[j1][j2];
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
		For accessed data - symmat- gathering
		*/
		{
			Ld=0;
			Ud=31;
			Id=(Ud-Ld);
			Rd=(Id%size);
			PSM=32;
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
				memcpy(rcvBuffsymmat, symmat, sizeof(symmat));
			}
			CD_S_p=((CD_U_p-CD_L_p)*PSM);
			MPI_Gather(&symmat[CD_L_p][0], CD_S_p, MPI_DOUBLE, &rcvBuffsymmat[CD_L_p][0], CD_S_p, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			if (rank==0)
			{
				memcpy(symmat, rcvBuffsymmat, sizeof(rcvBuffsymmat));
			}
		}
		/*
		For accessed data - symmat- gathering
		*/
		{
			Ld=0;
			Ud=31;
			Id=(Ud-Ld);
			Rd=(Id%size);
			PSM=32;
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
				memcpy(rcvBuffsymmat, symmat, sizeof(symmat));
			}
			CD_S_p=((CD_U_p-CD_L_p)*PSM);
			MPI_Gather(&symmat[CD_L_p][1], CD_S_p, MPI_DOUBLE, &rcvBuffsymmat[CD_L_p][1], CD_S_p, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			if (rank==0)
			{
				memcpy(symmat, rcvBuffsymmat, sizeof(rcvBuffsymmat));
			}
		}
		/*
		For accessed data - symmat- gathering
		*/
		{
			Ld=0;
			Ud=31;
			Id=(Ud-Ld);
			Rd=(Id%size);
			PSM=1;
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
				memcpy(rcvBuffsymmat, symmat, sizeof(symmat));
			}
			CD_S_p=((CD_U_p-CD_L_p)*PSM);
			MPI_Gather(&symmat[j2][CD_L_p], CD_S_p, MPI_DOUBLE, &rcvBuffsymmat[j2][CD_L_p], CD_S_p, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			if (rank==0)
			{
				memcpy(symmat, rcvBuffsymmat, sizeof(rcvBuffsymmat));
			}
		}
	}
	symmat[32-1][32-1]=1.0;
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
		print_array(m, symmat);
	}
	/* Be clean. */
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
=> TOTAL NO. OF PARALLEL LOOPS FOUND = 4
=> INNER PARALLEL LOOPS FOUND =0
=> OUTER PARALLEL LOOPS FOUND = 4
=>==========================================================
*/
