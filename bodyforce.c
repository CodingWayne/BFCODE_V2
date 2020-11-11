
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
		//-------------------------網格參數初始化----------------------
		meshload = 0; minCellid = 999999999, np = 0, threadboy = 0, * newsize, * threadnumber, all_size = 0, non = 0;
		//------------------------------------------------------------
		//-------------------------抓取羅盤網格------------------------
		if (ite == 0) 
		{
			debug("ite=0", getpid());
			meshload = 1;
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//output thread number & point site
			//依各核心序輸出包含核心之size,centroid[][],Velocity[][]文件
			sprintf(str, ".\\thread\\%d_%d.txt", getpid(), ite); //輸出各核心之size,centroid[][],Velocity[][]
			fp = fopen(str, "w");
			debug("open thread file", getpid());
			fprintf(fp, "%d\n", size);
			for (i = 0; i < size; i++) {
				fprintf(fp, "%lf %lf %lf %lf %lf %lf %lf\n", Cell_id[i][0], centroid[i][0], centroid[i][1], centroid[i][2], Velocity[i][0], Velocity[i][1], Velocity[i][2]);
				if (minCellid > (int)Cell_id[i][0]) {
					minCellid = (int)Cell_id[i][0];
				}
			}
			fclose(fp);

			Sleep(1000);

			debug("output minCellid file", getpid());
			sprintf(str, ".\\mesh\\%d_%d.txt", getpid(), ite); //輸出各核心之minCellid
			fp1 = fopen(str, "w");
			fprintf(fp1, "%d", minCellid);
			fclose(fp1);

			debug("sleeping", getpid());
			Sleep(1000);
			debug("sleep Done", getpid());

			//算np
			np = 0;
			np = count_file_num_in_a_folder(".\\thread");
			debug("np calculation Done!", getpid());

			//get everybody
			//以0開始向上增加來找最小核心序，並以最小核心序核心開啟檔案(allthread.txt)紀錄各核心之核心序和size，避免其他核心開啟此檔，其他核心跳過此步驟
			j = 0;
			for (i = 0; i < 99999999999; i++) {
				sprintf(str, ".\\thread\\%d_%d.txt", i, ite);
				fp = fopen(str, "r");
				if (fp)
				{
					if (i < getpid())
					{  //who is the min thread number then open file(allthread.txt)
						fclose(fp);
						break;
					}
					else
					{
						if (threadboy == 0)
						{
							threadboy = 1;
							sprintf(str, ".\\thread\\allthread%d.txt", getpid());
							fp1 = fopen(str, "w"); //輸出各核心之核心序、size
							newsize = (int*)malloc(np * sizeof(int));
							threadnumber = (int*)malloc(np * sizeof(int));
						}
					}
					Sleep(1000);

					fscanf(fp, "%d", &newsize[j]);
					threadnumber[j] = i; //newsize[j]儲存各核心之size threadnumber[j]=i 儲存各核心之核心序
					fprintf(fp1, "%d %d\n", threadnumber[j], newsize[j]);
					all_size = all_size + newsize[j]; //將各核心之size總和，之後會刪除重複項，使之等於圓盤網格數
					fclose(fp);
					j++;
					if (j == np) {
						fclose(fp1); break;
					}
				}
			}
			debug("Frist thread load done!", getpid());
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			if (threadboy == 1) {
				mesh = (double**)malloc((bfkin.plane_x_numb * 4 * bfkin.plane_c_numb * bfkin.plane_r_numb) * sizeof(double*));
				for (i = 0; i != (bfkin.plane_x_numb * 4 * bfkin.plane_c_numb * bfkin.plane_r_numb); ++i) {
					mesh[i] = (double*)malloc(8 * sizeof(double));


				}
				Sleep(1000);
				debug("malloc done!", getpid());

				//決定最小網格id並輸出
				for (i = 0; i < np; i++) {
					sprintf(str, ".\\mesh\\%d_%d.txt", threadnumber[i], ite);
					fp = fopen(str, "r");
					fscanf(fp, "%d", &numb);
					if (numb < minCellid) {
						minCellid = numb;
					}
					fclose(fp);
					sprintf(str, "del .\\mesh\\%d_%d.txt", threadnumber[i], ite);
					system(str);
				}

				Sleep(1000);

				debug("CPU data load!", getpid());

				sprintf(str, ".\\mesh\\%d_%d.txt", getpid(), ite); //輸出minCellid file
				fp1 = fopen(".\\mesh\\minCellid.txt", "w");
				fprintf(fp1, "%d", minCellid);
				fclose(fp1);

				for (i = 0; i < np; i++) {
					sprintf(str, ".\\thread\\%d_%d.txt", threadnumber[i], ite);
					fp = fopen(str, "r");
					fgets(str, 200, fp);
					for (j = 0; j < newsize[i]; j++) {
						fscanf(fp, "%lf", &non);
						kk = (int)non - minCellid;
						fscanf(fp, "%lf", &mesh[kk][0]);
						fscanf(fp, "%lf", &mesh[kk][1]);
						fscanf(fp, "%lf", &mesh[kk][2]);
						fscanf(fp, "%lf", &mesh[kk][3]);
						fscanf(fp, "%lf", &mesh[kk][4]);
						fscanf(fp, "%lf", &mesh[kk][5]);

					}
					fclose(fp);
					sprintf(str, "del .\\thread\\%d_%d.txt", threadnumber[i], ite);
					system(str);
				}
				Sleep(1000);
				debug("mesh change start!", getpid());

				for (i = 0; i != (bfkin.plane_x_numb * 4 * bfkin.plane_c_numb * bfkin.plane_r_numb); ++i) {

					mesh[i][0] = mesh[i][0] - bfkin.p_xo;
					mesh[i][1] = mesh[i][1] - bfkin.p_yo;
					mesh[i][2] = mesh[i][2] - bfkin.p_zo;


					if (bfkin.bow_dir == 1)

					{
						mesh[i][0] = mesh[i][0] * (-1);
						mesh[i][1] = mesh[i][1] * (-1);
						mesh[i][3] = mesh[i][3] * (-1);
						mesh[i][4] = mesh[i][4] * (-1);

					}

					mesh[i][6] = (double)pow(pow(mesh[i][1], 2.0) + pow(mesh[i][2], 2.0), 0.5);
					if (mesh[i][1] >= 0) {
						mesh[i][7] = (double)-acos(mesh[i][2] / mesh[i][6]); //徑度
					}
					else {
						mesh[i][7] = (double)acos(mesh[i][2] / mesh[i][6]); //徑度
					}
				}
				sprintf(str, "del .\\thread\\allthread%d.txt", threadnumber[0]);
				system(str);
				debug("mesh load done!", getpid());
				Sleep(1000);

				//輸出mesh[ID-1][C123_V123__r_theta]文本
				sprintf(str, ".\\mesh\\mesh.txt");
				fp = fopen(str, "w");
				for (i = 0; i < (bfkin.plane_x_numb * 4 * bfkin.plane_c_numb * bfkin.plane_r_numb); i++) {
					fprintf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf\n", mesh[i][0], mesh[i][1], mesh[i][2], mesh[i][3], mesh[i][4], mesh[i][5], mesh[i][6], mesh[i][7]);
				}
				fclose(fp);
				debug("mesh[ID-1][C123_V123__r_theta] output done!", getpid());


			}
		}
		
		
	}
		
	



	
}
