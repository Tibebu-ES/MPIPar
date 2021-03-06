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

 * mvt.c: This file is part of the PolyBench/C 3.2 test suite.
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
#include "mvt.h"
/* Array initialization. */
static void init_array(int n, double x1[(8000+0)], double x2[(8000+0)], double y_1[(8000+0)], double y_2[(8000+0)], double A[(8000+0)][(8000+0)])
{
	int i, j;
	for (i=0; i<n; i ++ )
	{
		x1[i]=(((double)i)/n);
		x2[i]=((((double)i)+1)/n);
		y_1[i]=((((double)i)+3)/n);
		y_2[i]=((((double)i)+4)/n);
		for (j=0; j<n; j ++ )
		{
			A[i][j]=((((double)i)*j)/8000);
		}
	}
	return ;
}

/*
DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output.
*/
static void print_array(int n, double x1[(8000+0)], double x2[(8000+0)])
{
	int i;
	for (i=0; i<n; i ++ )
	{
		fprintf(stderr, "%0.2lf ", x1[i]);
		fprintf(stderr, "%0.2lf ", x2[i]);
		if ((i%20)==0)
		{
			fprintf(stderr, "\n");
		}
	}
	return ;
}

/*
Main computational kernel. The whole function will be timed,
   including the call and return.
*/
/*
rcv buffer for accessed data - x1
*/
double rcvBuffx1[(8000+0)];
/*
rcv buffer for accessed data - x2
*/
double rcvBuffx2[(8000+0)];
int main(int argc, char * * argv)
{
	/* Retrieve problem size. */
	int n = 8000;
	/* Variable declarationallocation. */
	double A[(8000+0)][(8000+0)];
	double x1[(8000+0)];
	double x2[(8000+0)];
	double y_1[(8000+0)];
	double y_2[(8000+0)];
	/* Initialize array(s). */
	int i, j;
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
	init_array(n, x1, x2, y_1, y_2, A);
	/* Start timer. */
	;
	;
	/* Run kernel. */
	polybench_timer_start();
	#pragma scop 
	#pragma MPI  parallel 0
	#pragma MPI  accessed_data(x1[i])
	for (i=0; i<8000; i ++ )
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
		Uc=8000;
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
			for (j=0; j<8000; j ++ )
			{
				x1[i]=(x1[i]+(A[i][j]*y_1[j]));
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
		For accessed data - x1- gathering
		*/
		{
			Ld=0;
			Ud=8000;
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
				memcpy(rcvBuffx1, x1, sizeof(x1));
			}
			CD_S_p=((CD_U_p-CD_L_p)*PSM);
			MPI_Gather(&x1[CD_L_p], CD_S_p, MPI_DOUBLE, &rcvBuffx1[CD_L_p], CD_S_p, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			if (rank==0)
			{
				memcpy(x1, rcvBuffx1, sizeof(rcvBuffx1));
			}
		}
		mpipar_comm_timer_stop();
	}
	#pragma MPI  parallel 0
	#pragma MPI  accessed_data(x2[i])
	for (i=0; i<8000; i ++ )
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
		Uc=8000;
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
			for (j=0; j<8000; j ++ )
			{
				x2[i]=(x2[i]+(A[j][i]*y_2[j]));
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
		For accessed data - x2- gathering
		*/
		{
			Ld=0;
			Ud=8000;
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
				memcpy(rcvBuffx2, x2, sizeof(x2));
			}
			CD_S_p=((CD_U_p-CD_L_p)*PSM);
			MPI_Gather(&x2[CD_L_p], CD_S_p, MPI_DOUBLE, &rcvBuffx2[CD_L_p], CD_S_p, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			if (rank==0)
			{
				memcpy(x2, rcvBuffx2, sizeof(rcvBuffx2));
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
		print_array(n, x1, x2);
	}
	/* Be clean. */
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
=> TOTAL NO. OF PARALLEL LOOPS FOUND = 2
=> INNER PARALLEL LOOPS FOUND =0
=> OUTER PARALLEL LOOPS FOUND = 2
=>==========================================================
*/
