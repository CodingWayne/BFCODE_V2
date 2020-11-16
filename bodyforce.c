
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
		struct PROP prop;
		Load_bfkin(&bfkin);
		reitestar = bfkin.itestar % bfkin.ci;
		pitch = 0;
		//Output_bfkin(&bfkin, &reitestar, &pitch);

		ite = Iteration[0][0];
		int initmesh = 0;
		if (initmesh == 0)
		{
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
			initmesh = 1;
		}

		sprintf(str, ".\\debug\\debug%d.txt", getpid());
		fp = fopen(str, "w");
		fprintf(fp, "=================  Iteration=%d  =================\n", ite);
		fclose(fp);


		//-----------------------------------------------------------------------網格參數初始化----------------------
		meshload = 0; minCellid = 999999999, np = 0, threadboy = 0, * newsize, * threadnumber, all_size = 0, non = 0;
		//----------------------------------------------------------------------------------------------------------

		if (ite == 0 || ite >= bfkin.itestar && ite % bfkin.ci == reitestar && Iteration[0] != 0)
		{
			sprintf(str, ".\\debug\\debug%d.txt", getpid());
			fp = fopen(str, "w");
			fprintf(fp, "=================  bodyforcePre  =================\n", ite);
			fclose(fp);
			Output_bfkin(&bfkin, &reitestar, &pitch);
			threadboy = 0;
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

			//-----------------------------------------------------------------------Body Force開始---------------------
			/*此步驟為每輪Bodyforce開始前必要步驟
			*在每輪Body force 開始之前，1.輸出該輪螺盤網格上的速度 2.以單核心處理把新的速度更新到mesh[][]裡
			*在該輪的Body force結束後 會釋放mesh[][],以及threadboy歸0
			*/
			if (ite >= bfkin.itestar && ite % bfkin.ci == reitestar  && meshload == 0) 
			{
				MESH_STEP_1();	
				MESH_STEP_3();				
				MESH_STEP_4();
				if (threadboy == 1)
				{					
					MESH_STEP_5(&bfkin);
					//讀取minCellid
					fp = fopen(".\\mesh\\minCellid.txt", "r");
					fscanf(fp, "%d", &minCellid);
					fclose(fp);
					//讀取mesh.txt，並更新的網格速度資料
					sprintf(str, ".\\mesh\\mesh.txt");
					fp = fopen(str, "r");
					for (i = 0; i < (bfkin.plane_x_numb * 4 * bfkin.plane_c_numb * bfkin.plane_r_numb); i++) 
					{
						fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf", &mesh[i][0], &mesh[i][1], &mesh[i][2], &non, &non, &non, &mesh[i][6], &mesh[i][7]);
					}
					fclose(fp);
					debug("mesh.txt load Done!", getpid());
					for (i = 0; i < np; i++)
					{
						sprintf(str, ".\\thread\\%d_%d.txt", threadnumber[i], ite);
						fp = fopen(str, "r");

						debug(str, getpid());
						fp3 = fopen("IMD.dat", "w+");
						fgets(str, 200, fp);
						for (j = 0; j < newsize[i]; j++) {
							fscanf(fp, "%lf", &non);
							kk = (int)non - minCellid;
							fscanf(fp, "%lf", &non);
							fscanf(fp, "%lf", &non);
							fscanf(fp, "%lf", &non);
							fscanf(fp, "%lf", &mesh[kk][3]);
							fscanf(fp, "%lf", &mesh[kk][4]);
							fscanf(fp, "%lf", &mesh[kk][5]);
							fprintf(fp3, "%lf\t%lf\t%lf\t%lf\n", mesh[kk][0], mesh[kk][3], mesh[kk][4], mesh[kk][5]);
							if (bfkin.bow_dir == 1)
							{
								mesh[kk][3] = -1 * mesh[kk][3];
								mesh[kk][4] = -1 * mesh[kk][4];
								//mesh[kk][3]=-1*mesh[kk][3];
							}
						}
						fclose(fp);	fclose(fp3);
						sprintf(str, "del .\\thread\\%d_%d.txt", threadnumber[i], ite);
						system(str);
					}
					sprintf(str, "del .\\thread\\allthread%d.txt", threadnumber[0]);
					system(str);
					debug("mesh load done!", getpid());
				}
			}
			//
			//------------------------------------------------------------其他核心讀取mesh.txt---------------------------------
			//
			if (meshload == 0) 
			{
				if (threadboy == 0) 
				{
					mesh = (double**)malloc((bfkin.plane_x_numb * 4 * bfkin.plane_c_numb * bfkin.plane_r_numb) * sizeof(double*));
					for (i = 0; i != (bfkin.plane_x_numb * 4 * bfkin.plane_c_numb * bfkin.plane_r_numb); ++i) 
					{
						mesh[i] = (double*)malloc(8 * sizeof(double));
					}
					debug("mesh[][] done!", getpid());
					//讀取mesh.txt
					sprintf(str, ".\\mesh\\mesh.txt");
					fp = fopen(str, "r");
					debug("mesh.txt open!", getpid());
					for (i = 0; i < (bfkin.plane_x_numb * 4 * bfkin.plane_c_numb * bfkin.plane_r_numb); i++) 
					{
						fscanf(fp, "%lf%lf%lf%lf%lf%lf%lf%lf", &mesh[i][0], &mesh[i][1], &mesh[i][2], &non, &non, &non, &mesh[i][6], &mesh[i][7]);
					}
					fclose(fp);
					debug("mesh[][] load done!", getpid());
				}
			}
			if (ite > 0) {
				fp = fopen(".\\mesh\\minCellid.txt", "r");
				fscanf(fp, "%d", &minCellid);
				fclose(fp);
			}
			debug("=====  Pre     Done  =================", getpid());
			//------------------------------------------------------------讀取prop.geo----------------------------------------
			Load_PROP_GEO(&prop);//Output_PROP_GEO(&prop);			
			//------------------------------------------------------------讀取prop.adm----------------------------------------
			Load_PROP_ADM();//CHECK_PROP_ADM();
			
			

		}

		if (ite > bfkin.itestar && bfkin.Unsteady_offon == 0)
		{
			debug("00000", getpid());
		}

		
		
	}
		
	



	
}
