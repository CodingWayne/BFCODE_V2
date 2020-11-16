#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "uclib.h"
#include "parameter.h"

void Load_bfkin(struct BFKIN* bfkin)
{
	int i;
	FILE* fp;
	fp = fopen("bf.kin", "r+");
	fscanf(fp, "%d%d%d%d%d%d", &bfkin->np, &bfkin->itestar, &bfkin->ci, &bfkin->wake_offon, &bfkin->wake_ci, &bfkin->wake_onite);  //從bf->kin讀入; np:所使用的核心數; itestar:開始計算ci前的先行計算步數; ci:每幾步迭代一次; wake_offon:是否開啟wake輸出; wake_ci:每幾步輸出wake一次; wake_onite:開始輸出wake的步數; 
	fscanf(fp, "%d%d%lf%d%d", &bfkin->self_offon, &bfkin->firstJ_offon, &bfkin->firstJ, &bfkin->bodyforce_offon, &bfkin->P_3);		//self_offon:是否開啟self;	  firstJ_offon:是否開啟自定義第一次J值;  firstJ:自定義的第一次J值;   bodyforce_offon:是否開啟bodyforce;   
	fscanf(fp, "%d%d%d%lf%d%d%d", &bfkin->Unsteady_offon, &bfkin->unsteady_c, &bfkin->UN_ci, &bfkin->theta_one_blade, &bfkin->tnb, &bfkin->disk_offen, &bfkin->UNB);	     // Blade_offon:是否開啟葉片模式 ;  Unsteady_offon:是否開啟Unsteady bodyforce ; unsteady_c: ubf開始進行步數 ; theta_one_blade: 投影螺槳面積  ;  tnb: 每一次迭代螺槳葉面轉動圈數
	fscanf(fp, "%lf%lf%lf%d%d%d", &bfkin->R, &bfkin->x, &bfkin->xp, &bfkin->plane_r_numb, &bfkin->plane_x_numb, &bfkin->plane_c_numb); //R:圓盤半徑; x:圓盤厚度; xp:在x軸上的位置; plane_r_numb:圓盤r方向之網格分割數; plane_x_numb:圓盤x方向之網格分割數; plane_c_numb:圓盤周方向之網格分割數
	fscanf(fp, "%lf%lf%lf%lf%lf", &bfkin->Vs, &bfkin->rho, &bfkin->wr, &bfkin->ar, &bfkin->sfc);					//Vs:船速; rho:密度; WR:興波阻力; AR:波浪附加阻力; SFC:摩擦阻力修正量;
	fscanf(fp, "%d%lf%lf%lf%d", &bfkin->bow_dir, &bfkin->p_xo, &bfkin->p_yo, &bfkin->p_zo, &bfkin->st);	//圓盤中心座標 ; 船艏方向
	bfkin->bfkin_hub = (double*)malloc(bfkin->plane_x_numb * sizeof(double));
	for (i = 0; i != bfkin->plane_x_numb; ++i)fscanf(fp, "%lf", &bfkin->bfkin_hub[i]);  //x方向劃分之hub位置
	bfkin->st = bfkin->st * 1000;
	fclose(fp);
}

void Output_bfkin(struct BFKIN* bfkin, double* reitestar, double* pitch)
{
	int i;
	FILE* fp3;
	fp3 = fopen("Check_bfkin.dat", "w");
	fprintf(fp3, "np=%d\n", bfkin->np);
	fprintf(fp3, "itestart=%d\n", bfkin->itestar);
	fprintf(fp3, "ci=%d\n", bfkin->ci);
	fprintf(fp3, "wake_offon=%d\n", bfkin->wake_offon);
	fprintf(fp3, "wake_ci=%d\n", bfkin->wake_ci);
	fprintf(fp3, "wake_onite=%d\n", bfkin->wake_onite);
	fprintf(fp3, "firstJ_offon=%d\n", bfkin->firstJ_offon);
	fprintf(fp3, "firstJ=%lf\n", bfkin->firstJ);
	fprintf(fp3, "R=%lf\n", bfkin->R);
	fprintf(fp3, "x=%lf\n", bfkin->x);
	fprintf(fp3, "xp=%lf\n", bfkin->xp);
	fprintf(fp3, "P_3=%d\n", bfkin->P_3);
	fprintf(fp3, "Unsteady_offon=%d\n", bfkin->Unsteady_offon);
	fprintf(fp3, "disk_offen=%d\n", bfkin->disk_offen);
	fprintf(fp3, "plane_r_numb=%d\n", bfkin->plane_r_numb);
	fprintf(fp3, "plane_x_numb=%d\n", bfkin->plane_x_numb);
	fprintf(fp3, "plane_c_numb=%d\n", bfkin->plane_c_numb);
	fprintf(fp3, "Vs=%lf\n", bfkin->Vs);
	fprintf(fp3, "rho=%lf\n", bfkin->rho);
	fprintf(fp3, "Wave Resistance=%lf\n", bfkin->wr);
	fprintf(fp3, "Added Resistance=%lf\n", bfkin->ar);
	fprintf(fp3, "SFC=%lf\n", bfkin->sfc);
	fprintf(fp3, "reitestar=%d\n", *reitestar);
	fprintf(fp3, "bfkin_hub={");
	for (i = 0; i != bfkin->plane_x_numb; ++i) {
		fprintf(fp3, " %lf", bfkin->bfkin_hub[i]);
	}
	fprintf(fp3, "}\n");
	fprintf(fp3, "Pitch=%f\n", *pitch);
	fprintf(fp3, "bow_dir=%d\n", bfkin->bow_dir);
	fprintf(fp3, "P_X = %lf\tP_Y = %lf\tP_Z = %lf\n", bfkin->p_xo, bfkin->p_yo, bfkin->p_zo);
	fclose(fp3);
}

void Load_PROP_ADM()
{
	fp = fopen("prop.adm", "r");
	fgets(adm1, 100, fp);
	fgets(adm2, 100, fp);
	fgets(adm3, 100, fp);
	fgets(adm4, 100, fp);
	fgets(adm5, 100, fp);
	fgets(adm6, 100, fp);
	fclose(fp);
	debug("prop.adm load!", getpid());
}
void CHECK_PROP_ADM()
{
	fp3 = fopen("Check_propadm.dat", "w");
	fprintf(fp3, adm1);
	fprintf(fp3, adm2);
	fprintf(fp3, adm3);
	fprintf(fp3, adm4);
	fprintf(fp3, adm5);
	fprintf(fp3, adm6);
	fclose(fp3);
	debug("prop.adm check!", getpid());
}
void Load_PROP_GEO(struct PROP* prop)
{
	sprintf(geoname, "prop.geo");
	fp = fopen(geoname, "r");
	fgets(label, 100, fp);
	fscanf(fp, "%d%d%d%d%d", &prop->NX, &prop->NBLADE, &prop->NC, &prop->MR, &prop->NTMP);
	fscanf(fp, "%lf%d%d%lf%lf%lf", &prop->RHUB, &prop->NHBU, &prop->MHBT, &prop->XHBU, &prop->XHBD, &prop->XHBT);
	fscanf(fp, "%lf%lf%lf%lf%lf%lf%lf", &prop->ADVCO, &prop->RULT, &prop->RHULT, &prop->DCD, &prop->XULT, &prop->DTPROP, &prop->XUWDK);
	prop->NPARAMETER = (double**)malloc(14 * sizeof(double*)); //用以儲存prop.geo第四行之後的數值 
	for (i = 0; i != 14; i++)
	{
		prop->NPARAMETER[i] = (double*)malloc(prop->NX * sizeof(double));
	}
	for (j = 0; j != 14; ++j)
	{
		for (i = 0; i != prop->NX; ++i)
		{
			fscanf(fp, "%lf", &prop->NPARAMETER[j][i]);
		}
					
	}
		
	fclose(fp);
	debug("*.geo load!", getpid());
}

void Output_PROP_GEO(struct PROP* prop)
{
	sprintf(geoname, "propgeo_check.dat");
	fp = fopen(geoname, "w");
	fprintf(fp, "prop\n");
	fprintf(fp, "%d %d %d %d %d \n", prop->NX, prop->NBLADE, prop->NC, prop->MR, prop->NTMP);
	fprintf(fp, "%.3lf %d %d %.4lf %.4lf %.4lf\n", prop->RHUB, prop->NHBU, prop->MHBT, prop->XHBU, prop->XHBD, prop->XHBT);
	fprintf(fp, "%.3lf %.3lf %.3lf %.3lf %.3lf %.3lf %.3lf\n", prop->ADVCO, prop->RULT, prop->RHULT, prop->DCD, prop->XULT, prop->DTPROP, prop->XUWDK);
	for (i = 0; i < 14;i++)
	{
		for (j = 0; j < prop->NX; j++)
		{
			fprintf(fp, "%.5lf\t", prop->NPARAMETER[i][j]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
	debug("*.geo Output!", getpid());
}

void MESH_STEP_1()
{
	/*
	 *STEP[1]
	 *使用getpid輸出各process獨自處理的網格資料(thread_pid_ite.txt)
	 *假設有5個process在分工處理，thread資料夾便會輸出5個(thread_pid_ite.txt)對應5個不同的pid
	 *(thread_pid_ite.txt)=[網格編號,網格中心點X,網格中心點y,網格中心點z,網格中心點速度X,網格中心點速度y,網格中心點速度z]
	 *此時minCellid在各別的process裡記錄最小的網格編號 (minCellid不是矩陣只是單變數，但在各個process裡代表該process的minCellid)
	 */
	FILE* fp;
	sprintf(str, ".\\thread\\%d_%d.txt", getpid(), ite);
	fp = fopen(str, "w");
	debug("STEP[1]-1   open thread file", getpid());
	fprintf(fp, "%d\n", tmp_size);
	for (i = 0; i < tmp_size; i++)
	{
		fprintf(fp, "%lf %lf %lf %lf %lf %lf %lf\n", tmp_Cell_id[i][0], tmp_centroid[i][0], tmp_centroid[i][1], tmp_centroid[i][2], tmp_Velocity[i][0], tmp_Velocity[i][1], tmp_Velocity[i][2]);
		if (minCellid > (int)tmp_Cell_id[i][0])
		{
			minCellid = (int)tmp_Cell_id[i][0];
		}
	}
	fclose(fp);
	debug("STEP[1]-2   output minCellid", getpid());
}

void MESH_STEP_2()
{
	/*
	*STEP[2]
	*輸出各核心之minCellid
	*/
	FILE* fp1;
	sprintf(str, ".\\mesh\\%d_%d.txt", getpid(), ite);
	fp1 = fopen(str, "w");
	fprintf(fp1, "%d", minCellid);
	fclose(fp1);
	debug("STEP[2]     output minCellid file", getpid());
}

void MESH_STEP_3()
{
	/*
	 *STEP[3]
	 *計算thread資料夾下有多少檔案=有多少process在處理網格
	 */
	np = 0;
	np = count_file_num_in_a_folder(".\\thread");
	debug("STEP[3]     np calculation Done!", getpid());
	
}
void MESH_STEP_4()
{
	/*
	 *STEP[4]
	 *
	 *以0開始向上增加來找最小核心序，並以最小核心序核心開啟檔案(allthread.txt)紀錄各核心之核心序和size，避免其他核心開啟此檔，其他核心跳過此步驟
	 *threadboy=1表示已在最小process裡建置矩陣成功並以最小核心開啟各檔案紀錄在矩陣裡
	 */
	j = 0;
	FILE* fp,*fp1;
	for (i = 0; i < 99999999999; i++)
	{
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
			fscanf(fp, "%d", &newsize[j]);
			threadnumber[j] = i; //newsize[j]儲存各核心之size threadnumber[j]=i 儲存各核心之核心序
			fprintf(fp1, "%d %d\n", threadnumber[j], newsize[j]);
			all_size = all_size + newsize[j]; //將各核心之size總和，之後會刪除重複項，使之等於圓盤網格數
			fclose(fp);
			j++;
			if (j == np)
			{
				fclose(fp1); break;
			}
		}
	}
	debug("STEP[4]     Frist thread load done!", getpid());
}
void MESH_STEP_5(struct BFKIN* bfkin)
{
	/*
	 *STEP[5]
	 * 建置2維mesh矩陣  mesh[圓盤網格數][8]
	 */
	mesh = (double**)malloc((bfkin->plane_x_numb * 4 * bfkin->plane_c_numb * bfkin->plane_r_numb) * sizeof(double*));
	for (i = 0; i != (bfkin->plane_x_numb * 4 * bfkin->plane_c_numb * bfkin->plane_r_numb); ++i)
	{
		mesh[i] = (double*)malloc(8 * sizeof(double));
	}
	debug("STEP[5]     MESH[][] malloc done!", getpid());
}
void MESH_STEP_6()
{
	/*
	* STEP[6]
	* 決定最小網格id並輸出,輸出後刪除
	*/
	FILE* fp,* fp1;
	for (i = 0; i < np; i++)
	{
		sprintf(str, ".\\mesh\\%d_%d.txt", threadnumber[i], ite);
		fp = fopen(str, "r");
		fscanf(fp, "%d", &numb);
		if (numb < minCellid)
		{
			minCellid = numb;
		}
		fclose(fp);
		sprintf(str, "del .\\mesh\\%d_%d.txt", threadnumber[i], ite);
		system(str);
	}
	sprintf(str, ".\\mesh\\%d_%d.txt", getpid(), ite); //輸出minCellid file
	fp1 = fopen(".\\mesh\\minCellid.txt", "w");
	fprintf(fp1, "%d", minCellid);
	fclose(fp1);
	debug("STEP[6]     CPU data load!", getpid());
}
void MESH_STEP_7()
{
	/*
	* STEP[7]
	* 放入mesh[][]
	* 依以前紀錄之np,依序開啟各thread輸出之(thread_pid_ite.txt),假設網格數樹有22500個,對應編號須為0~22499以方便日後操作矩陣,因此kk需為22500-minCellid(1)
	* 將(thread_pid_ite.txt)裡各Cellid對應的網格資料放到mesh[id-1][x,y,z,Vx,Vy,Vz,?,?]裡,此時各(thread_pid_ite.txt)裡重複之Cellid資料會覆蓋,可得完整22500個網格資料之矩陣
	* 刪除(thread_pid_ite.txt)
	*/
	FILE* fp;
	for (i = 0; i < np; i++)
	{
		sprintf(str, ".\\thread\\%d_%d.txt", threadnumber[i], ite);
		fp = fopen(str, "r");
		fgets(str, 200, fp);
		for (j = 0; j < newsize[i]; j++)
		{
			fscanf(fp, "%lf", &non);
			kk = (int)non - minCellid;
			fscanf(fp, "%lf", &mesh[kk][0]);// Centroid_X
			fscanf(fp, "%lf", &mesh[kk][1]);// Centroid_Y
			fscanf(fp, "%lf", &mesh[kk][2]);// Centroid_Z
			fscanf(fp, "%lf", &mesh[kk][3]);// Vx
			fscanf(fp, "%lf", &mesh[kk][4]);// Vy
			fscanf(fp, "%lf", &mesh[kk][5]);// Vz
		}
		fclose(fp);
		sprintf(str, "del .\\thread\\%d_%d.txt", threadnumber[i], ite);
		system(str);
	}
	debug("STEP[7]     mesh change start!", getpid());
}
void MESH_STEP_8(struct BFKIN* bfkin)
{
	/*
	*STEP[8]
	* 建立mesh矩陣裡螺盤徑向的距離及徑度 mesh[22500][6]&mesh[22500][7]
	*/
	for (i = 0; i != (bfkin->plane_x_numb * 4 * bfkin->plane_c_numb * bfkin->plane_r_numb); ++i)
	{
		/*
		* p_xo,p_yo,p_zo為螺盤網格中心點,此時移回原點
		*/
		mesh[i][0] = mesh[i][0] - bfkin->p_xo;
		mesh[i][1] = mesh[i][1] - bfkin->p_yo;
		mesh[i][2] = mesh[i][2] - bfkin->p_zo;

		//根據推力方向 修改網格方向
		if (bfkin->bow_dir == 1)
		{
			mesh[i][0] = mesh[i][0] * (-1);
			mesh[i][1] = mesh[i][1] * (-1);
			mesh[i][3] = mesh[i][3] * (-1);
			mesh[i][4] = mesh[i][4] * (-1);
		}
		//mesh[i][6]為螺盤徑向的距離r
		//mesh[i][7]為螺盤徑向的徑度theta
		mesh[i][6] = (double)pow(pow(mesh[i][1], 2.0) + pow(mesh[i][2], 2.0), 0.5);
		if (mesh[i][1] >= 0)
		{
			mesh[i][7] = (double)-acos(mesh[i][2] / mesh[i][6]); //徑度
		}
		else 
		{
			mesh[i][7] = (double)acos(mesh[i][2] / mesh[i][6]); //徑度
		}
	}
	sprintf(str, "del .\\thread\\allthread%d.txt", threadnumber[0]);
	system(str);
	debug("STEP[8]     mesh load done!", getpid());
}
void MESH_STEP_9(struct BFKIN* bfkin)
{
	/*
	*STEP[9]
	*輸出mesh[ID-1][C123_V123__r_theta]文本
	*/
	FILE* fp;
	sprintf(str, ".\\mesh\\mesh.txt");
	fp = fopen(str, "w");
	for (i = 0; i < (bfkin->plane_x_numb * 4 * bfkin->plane_c_numb * bfkin->plane_r_numb); i++)
	{
		fprintf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf\n", mesh[i][0], mesh[i][1], mesh[i][2], mesh[i][3], mesh[i][4], mesh[i][5], mesh[i][6], mesh[i][7]);
	}
	fclose(fp);
	debug("STEP[9]     mesh[ID-1][Cx_Cy_Cz_Vx_Vy_Vz_r_theta] output done!", getpid());
}
