
#include "uclib.h"
#include <math.h>
#include <stdio.h> 
#include <stdlib.h>
#include <string.h>
#include <windows.h>
#include <time.h>
#include <omp.h>
#include "OSfunc.h"
#include <process.h>
#include "parameter.h"



//[result]:the array of values returned by the user function.(Force/Volume);  [size]:the number of elements in the result array.(int);  [centroid]:網格中心點之座標;  [Velocity]:網格中心點之速度; [Iteration]:the iteration number in CoordReal precision.(int); [Angle]:船體pitch的徑度(艏仰為正); [Cell_id]:圓盤網格之編號; [Drag]:船殼阻力;
void USERFUNCTION_EXPORT bodyforce(Real(*result)[3], int size, CoordReal(*centroid)[3], CoordReal(*Velocity)[3], CoordReal(*Cell_id)[1], CoordReal(*Iteration)[1], CoordReal(*Angle)[1])
{
	
	if (size != 0)
	{
		
		
		
		struct BFKIN bfkin;
		FILE* fp5, * fp3, * fp, * fp1;
		Load_bfkin(&bfkin);
		reitestar = bfkin.itestar % bfkin.ci;
		pitch = 0;
		Output_bfkin(&bfkin, &reitestar, &pitch);

		ite = Iteration[0][0];

		tmp_size = size;
		tmp_Cell_id = (double*)malloc(tmp_size * sizeof(double*));
		for (i = 0; i < tmp_size; i++)
		{
			tmp_Cell_id[i] = (double*)malloc(1 * sizeof(double));
		}
		for (i = 0; i < tmp_size; i++)
		{
			tmp_Cell_id[i][0] = Cell_id[i][0];
		}

		tmp_centroid = (double*)malloc(tmp_size * sizeof(double*));
		for (i = 0; i < tmp_size; i++)
		{
			tmp_centroid[i] = (double*)malloc(3 * sizeof(double));
		}
		for (i = 0; i < tmp_size; i++)
		{
			for (j = 0; j < 3; j++)
			{
				tmp_centroid[i][j] = centroid[i][j];
			}
		}

		tmp_Velocity = (double*)malloc(tmp_size * sizeof(double*));
		for (i = 0; i < tmp_size; i++)
		{
			tmp_Velocity[i] = (double*)malloc(3 * sizeof(double));
		}
		for (i = 0; i < tmp_size; i++)
		{
			for (j = 0; j < 3; j++)
			{
				tmp_Velocity[i][j] = Velocity[i][j];
			}
		}
		debug("------------COPY[][] done!", getpid());


		//-----------------------------------------------------------------------網格參數初始化----------------------
		meshload = 0; minCellid = 999999999, np = 0, threadboy = 0, * newsize, * threadnumber, all_size = 0, non = 0;
		//----------------------------------------------------------------------------------------------------------

		//-----------------------------------------------------------------------抓取螺盤網格------------------------
		if (ite == 0)
		{	
			meshload = 1;//meshload=1 表示已抓取網格
			MESH_STEP_1();
			MESH_STEP_2();
			MESH_STEP_3();
			MESH_STEP_4();
			if (threadboy == 1)
			{
				MESH_STEP_5(&bfkin);
				MESH_STEP_6();
				MESH_STEP_7();
				MESH_STEP_8(&bfkin);
				MESH_STEP_9(&bfkin);
			}

			
		}

		//----------------------------------------------------------------------------------------------------------
		
		
	}
		
	



	
}
