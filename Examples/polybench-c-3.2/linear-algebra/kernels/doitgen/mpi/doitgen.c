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

 * doitgen.c: This file is part of the PolyBench/C 3.2 test suite.
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
#include "doitgen.h"
/* Array initialization. */
static void init_array(int nr, int nq, int np, double A[(256+0)][(256+0)][(256+0)], double C4[(256+0)][(256+0)])
{
	int i, j, k;
	for (i=0; i<nr; i ++ )
	{
		for (j=0; j<nq; j ++ )
		{
			for (k=0; k<np; k ++ )
			{
				A[i][j][k]=(((((double)i)*j)+k)/np);
			}
		}
	}
	for (i=0; i<np; i ++ )
	{
		for (j=0; j<np; j ++ )
		{
			C4[i][j]=((((double)i)*j)/np);
		}
	}
	return ;
}

/*
DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output.
*/
static void print_array(int nr, int nq, int np, double A[(256+0)][(256+0)][(256+0)])
{
	int i, j, k;
	for (i=0; i<nr; i ++ )
	{
		for (j=0; j<nq; j ++ )
		{
			for (k=0; k<np; k ++ )
			{
				fprintf(stderr, "%0.2lf ", A[i][j][k]);
				if ((i%20)==0)
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
rcv buffer for accessed data - sum
*/
double rcvBuffsum[(256+0)][(256+0)][(256+0)];
/*
rcv buffer for accessed data - A
*/
double rcvBuffA[(256+0)][(256+0)][(256+0)];
int main(int argc, char * * argv)
{
	/* Retrieve problem size. */
	int nr = 256;
	int nq = 256;
	int np = 256;
	/* Variable declarationallocation. */
	double A[(256+0)][(256+0)][(256+0)];
	double sum[(256+0)][(256+0)][(256+0)];
	double C4[(256+0)][(256+0)];
	/* Initialize array(s). */
	int r, q, p, s;
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
	init_array(nr, nq, np, A, C4);
	/* Start timer. */
	;
	;
	/* Run kernel. */
	polybench_timer_start();
	#pragma scop 
	#pragma MPI  parallel 0
	#pragma MPI  accessed_data(sum[r][q][p],A[r][q][p])
	for (r=0; r<256; r ++ )
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
		Uc=256;
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
		if ((r>=Lcp)&&(r<Ucp))
		{
			for (q=0; q<256; q ++ )
			{
				for (p=0; p<256; p ++ )
				{
					sum[r][q][p]=0;
					for (s=0; s<256; s ++ )
					{
						sum[r][q][p]=(sum[r][q][p]+(A[r][q][s]*C4[s][p]));
					}
				}
				for (p=0; p<256; p ++ )
				{
					A[r][q][p]=sum[r][q][p];
				}
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
		mpipar_comm_timer_start();
		/*
		For accessed data - sum- gathering
		*/
		{
			Ld=0;
			Ud=256;
			Id=(Ud-Ld);
			Rd=(Id%size);
			PSM=65536;
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
				memcpy(rcvBuffsum, sum, sizeof(sum));
			}
			CD_S_p=((CD_U_p-CD_L_p)*PSM);
			MPI_Gather(&sum[CD_L_p][0][0], CD_S_p, MPI_DOUBLE, &rcvBuffsum[CD_L_p][0][0], CD_S_p, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			if (rank==0)
			{
				memcpy(sum, rcvBuffsum, sizeof(rcvBuffsum));
			}
		}
		/*
		For accessed data - A- gathering
		*/
		{
			Ld=0;
			Ud=256;
			Id=(Ud-Ld);
			Rd=(Id%size);
			PSM=65536;
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
				memcpy(rcvBuffA, A, sizeof(A));
			}
			CD_S_p=((CD_U_p-CD_L_p)*PSM);
			MPI_Gather(&A[CD_L_p][0][0], CD_S_p, MPI_DOUBLE, &rcvBuffA[CD_L_p][0][0], CD_S_p, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			if (rank==0)
			{
				memcpy(A, rcvBuffA, sizeof(rcvBuffA));
			}
		}
		mpipar_comm_timer_stop();
	}
	polybench_timer_stop();
	#pragma endscop 
	if (rank==0)
	{
		mpipar_comm_overhead_print();
	}
	/* Stop and print timer. */
	;
	;
	;
	;
	/*
	Prevent dead-code elimination. All live-out data must be printed
	     by the function call in argument.
	*/
	if ((argc>42)&&( ! strcmp(argv[0], "")))
	{
		print_array(nr, nq, np, A);
	}
	/* Be clean. */
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
