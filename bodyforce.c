
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

		reitestar = bfkin.itestar % bfkin.ci;
		pitch = 1;
		bfkin.R = 1;
		//Output_bfkin(&bfkin, &reitestar, &pitch);
		PI = (double)acos(-1.0);
		ite = Iteration[0][0];
		pitch = (double)Angle[0][0]; //angle[0] 當下的pitch徑度
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
		Load_bfkin(&bfkin);
		sprintf(str, ".\\debug\\debug%d.txt", getpid());
		fp = fopen(str, "w");
		fprintf(fp, "=================  Iteration=%d  =================\n", ite);
		fclose(fp);
		if (ite == 0) 
		{
			system("mkdir debug");
			system("mkdir thread");
			system("mkdir mesh");
			system("mkdir propeller");
			system("mkdir wake");
		}

		//=======================================================================網格參數初始化======================
		meshload = 0; minCellid = 999999999, np = 0, threadboy = 0, * newsize, * threadnumber, all_size = 0, non = 0;
		//==========================================================================================================

		if (ite == 0 || ite >= bfkin.itestar && ite % bfkin.ci == reitestar && Iteration[0] != 0)
		{
			Output_bfkin(&bfkin, &reitestar, &pitch);
			threadboy = 0;
			//=======================================================================抓取螺盤網格========================
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
			if (bfkin.self_offon == 1)
			{
				//=======================================================================Body Force開始=====================
				/*此步驟為每輪Bodyforce開始前必要步驟
				*在每輪Body force 開始之前，1.輸出該輪螺盤網格上的速度 2.以單核心處理把新的速度更新到mesh[][]裡
				*在該輪的Body force結束後 會釋放mesh[][],以及threadboy歸0
				*/
				if (ite >= bfkin.itestar && ite % bfkin.ci == reitestar && meshload == 0)
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
				//============================================================其他核心讀取mesh.txt=================================
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
				if (ite > 0)
				{
					fp = fopen(".\\mesh\\minCellid.txt", "r");
					fscanf(fp, "%d", &minCellid);
					fclose(fp);
				}
				//============================================================讀取prop.geo========================================
				Load_PROP_GEO(&prop);//Output_PROP_GEO(&prop);			
				//============================================================讀取prop.adm========================================
				Load_PROP_ADM();//CHECK_PROP_ADM();
				//============================================================配置指標============================================
				Vx_ave = (double*)malloc(bfkin.plane_r_numb * sizeof(double));
				Vt_ave = (double*)malloc(bfkin.plane_r_numb * sizeof(double));//Vx_ave,Vt_ave,Vr_ave:用以儲存各徑向位置網格的周向平均速度
				Vr_ave = (double*)malloc(bfkin.plane_r_numb * sizeof(double));
				Vx_sum = (double*)malloc(bfkin.plane_r_numb * sizeof(double));//Vx_sum,Vt_sum,Vr_sum:用以儲存各徑向位置網格的周向速度總和
				Vt_sum = (double*)malloc(bfkin.plane_r_numb * sizeof(double));
				Vr_sum = (double*)malloc(bfkin.plane_r_numb * sizeof(double));
				r_ave = (double*)malloc(bfkin.plane_r_numb * sizeof(double)); // r_ave:用以儲存各徑向位置網格的平均半徑大小
				r_sum = (double*)malloc(bfkin.plane_r_numb * sizeof(double)); // r_sum:用以儲存各徑向位置網格的半徑大小總和
				n = (double*)malloc(bfkin.plane_r_numb * sizeof(double));     // n:用以儲存各徑向位置之周向網格的數量
				//使用malloc()函式配置動態的空間 plane_r_numb:number of spanwise  mesh in Disk        MR:number of spanwise panels in patpans11

				fx = (double*)malloc(prop.MR * sizeof(double));         //fx[],fr[],ft[]:patpans11.exe計算之力?
				fr = (double*)malloc(prop.MR * sizeof(double));
				ft = (double*)malloc(prop.MR * sizeof(double));
				bfx = (double*)malloc(prop.MR * sizeof(double));		//bfx[],bfr[],bft[] = fx~r~t[i]/0.5*rho*Vs*Vs*(r_x[i+1]-r_x[i])*R;
				bfr = (double*)malloc(prop.MR * sizeof(double));
				bft = (double*)malloc(prop.MR * sizeof(double));
				r_F = (double*)malloc(prop.MR * sizeof(double));		//r_F[]:用以儲存patpans11.exe計算之力的徑向位置，數量為MR個，(垂向分段為sine spanwise) 
				r_x = (double*)malloc((prop.MR + 1) * sizeof(double));  //r_x[]:用以儲存葉片垂向分段(MR個段)裡，包含頭尾之每點的r/R，數量為MR+1個，(垂向分段為sine spanwise) 
				r_U = (double*)malloc(prop.MR * sizeof(double));		//r_U[]:用以儲存sfpv11.exe計算之誘導速度的徑向位置，數量為MR個，(垂向分段為sine spanwise) 
				Vx_total = (double*)malloc(prop.NX * sizeof(double));   //Vx_total[],Vt_total[],Vr_total[]
				Vt_total = (double*)malloc(prop.NX * sizeof(double));   //用以儲存各徑向位置網格之周向平均速度  SPLINE成 NX 個點後的結果，當作圓盤中心面之入流速度
				Vr_total = (double*)malloc(prop.NX * sizeof(double));
				Ux = (double*)malloc(prop.MR * sizeof(double));			//Ux[],Ut[],Ur[]:用以儲存sfpv11.exe計算之各徑向位置的誘導速度與船速的比值(Vinduced/Vs)，數量為MR個
				Ut = (double*)malloc(prop.MR * sizeof(double));
				Ur = (double*)malloc(prop.MR * sizeof(double));
				Ux_bem = (double*)malloc(prop.NX * sizeof(double));		//Ux_bem[],Ut_bem[],Ur_bem[]
				Ut_bem = (double*)malloc(prop.NX * sizeof(double));     //用以儲存sfpv11.exe計算之各徑向位置的誘導速度與船速的比值(Vinduced/Vs) 
				Ur_bem = (double*)malloc(prop.NX * sizeof(double));		//SPLINE成 NX 個點後的結果，當作圓盤中心面之誘導速度
				debug("=====  Pre     Done  =================", getpid());
				//======================================================計算間距==============================================================
				for (i = 0; i != prop.MR + 1; ++i)
				{
					r_x[i] = prop.RHUB + (1 - prop.RHUB) * (double)sin(0.5 * PI / prop.MR * i);
					//r_x[]:用以儲存葉片垂向分段(MR個段)裡，包含頭尾之每點的r/R，數量為MR+1個，(垂向分段為sine spanwise) 
				}
				//======================================================體積力歸0=============================================================
				if (ite <= bfkin.itestar)
				{
					//第一次執行patpans11前，所放入的單位體積力為0
					for (i = 0; i != size; ++i)
					{
						result[i][0] = 0.0;
						result[i][1] = 0.0;
						result[i][2] = 0.0;
					}
					debug("result 0 input!", getpid());
				}
				/****************************************************************************************************************************************
				*********																														*********
				*********											  開始計算體積														    *********
				*********																										(以單核心運作)	*********
				****************************************************************************************************************************************/
				if (ite >= bfkin.itestar && ite % bfkin.ci == reitestar && Iteration[0] != 0 && threadboy == 1)
				{
					//================================提取total_velocity並做周向平均(徑向速度忽略不計)======================================================

					for (i = 0; i != bfkin.plane_r_numb; ++i)
					{
						Vx_sum[i] = 0.0;
						Vt_sum[i] = 0.0;
						Vr_sum[i] = 0.0;
						r_sum[i] = 0.0;
						n[i] = 0.0;
					}
					//以序提取網格，計算其半徑大小,徑度，以其徑向位置作周向加總
					for (i = 0; i != (bfkin.plane_x_numb * 4 * bfkin.plane_c_numb * bfkin.plane_r_numb); ++i)
					{
						//計算網格所在的半徑大小
						r = mesh[i][6];

						//利用網格所在的x值，計算theta徑度
						theta = mesh[i][7];
						for (j = 0; j != bfkin.plane_r_numb; ++j)
						{
							//判斷是否為圓盤中心線上的網格，並判斷是哪個徑向位置的網格
							//將各徑向位置的網格，周向加總其all_Velocity[][0~2],徑向位置r,數量n
							if (r >= bfkin.R * prop.RHUB + j * bfkin.R * (1 - prop.RHUB) / bfkin.plane_r_numb && r <= bfkin.R * prop.RHUB + (j + 1) * bfkin.R * (1 - prop.RHUB) / bfkin.plane_r_numb && mesh[i][0] >= (bfkin.xp - bfkin.x / 2 / bfkin.plane_r_numb) && mesh[i][0] <= (bfkin.xp + bfkin.x / 2 / bfkin.plane_r_numb))
							{
								Vx_sum[j] = Vx_sum[j] + mesh[i][3] * cos(pitch) - mesh[i][5] * sin(pitch);
								Vr_sum[j] = Vr_sum[j] - mesh[i][4] * sin(theta) + (mesh[i][3] * sin(pitch) + mesh[i][5] * cos(pitch)) * cos(theta);
								Vt_sum[j] = Vt_sum[j] - mesh[i][4] * cos(theta) - (mesh[i][3] * sin(pitch) + mesh[i][5] * cos(pitch)) * sin(theta);
								r_sum[j] = r_sum[j] + r;
								n[j] = n[j] + 1.0;
							}
						}
					}

					debug("Velocity load!", getpid());

					//將各總和值，以各徑向位置的網格數量做周向平均
					for (i = 0; i != bfkin.plane_r_numb; ++i)
					{
						Vx_ave[i] = Vx_sum[i] / n[i];
						Vt_ave[i] = Vt_sum[i] / n[i];
						Vr_ave[i] = Vr_sum[i] / n[i];
						r_ave[i] = r_sum[i] / n[i];
					}
					debug("Velocity Average!", getpid());

					sprintf(str, "Vxrt_ave%d.txt", ite);
					fp = fopen(str, "w");
					for (i = 0; i != bfkin.plane_r_numb; ++i)
					{
						fprintf(fp, "%lf\n", Vx_ave[i]);
						fprintf(fp, "%lf\n", Vt_ave[i]);
						fprintf(fp, "%lf\n", Vr_ave[i]);
						fprintf(fp, "%lf\n", r_ave[i]);
					}
					fclose(fp);
					debug("Velocity Average! Output done!", getpid());

					//==========================================================將各徑向位置之周向平均速度，藉由SPLINE.exe化為 NX 個點===============================================
					//==========================================================點位置由XR:radii of the input geometry parameters for patpans11 決定==============================
					//==========================================================將plane_r_numb個徑向之周向平均速度-----(SPLINE)---->為NX個速度======================================
					fp = fopen("input.dat", "w");
					for (i = 0; i != bfkin.plane_r_numb; ++i)
					{
						fprintf(fp, "%10.6lf %10.6lf\n", r_ave[i] / bfkin.R, Vx_ave[i]);
					}
					fclose(fp);

					fp = fopen("SPLINE.KIN", "w");
					fprintf(fp, "%d %d\ninput.dat\noutputx.dat", bfkin.plane_r_numb, prop.NX);
					fclose(fp);

					fp = fopen("outputx.dat", "w");
					for (i = 0; i != prop.NX; ++i)
					{
						fprintf(fp, "%8.5lf\n", prop.NPARAMETER[0][i]);// NX 個點位置由XR:radii of the input geometry parameters for patpans11 決定
					}
					fclose(fp);

					sprintf(cmd, "SPLINE.exe");
					outexe(cmd);
					Sleep(100);
					fp = fopen("spline.dat", "r");
					for (i = 0; i != prop.NX; ++i)
					{
						fscanf(fp, "%lf%lf", &non, &Vx_total[i]);
					}
					fclose(fp);
					////////////////////////////////////////////////////////////////////////
					fp = fopen("input.dat", "w");
					for (i = 0; i != bfkin.plane_r_numb; ++i)
					{
						fprintf(fp, "%10.6lf %10.6lf\n", r_ave[i] / bfkin.R, Vt_ave[i]);
					}
					fclose(fp);

					fp = fopen("SPLINE.KIN", "w");
					fprintf(fp, "%d %d\ninput.dat\noutputx.dat", bfkin.plane_r_numb, prop.NX);
					fclose(fp);


					sprintf(cmd, "SPLINE.exe");
					outexe(cmd);
					Sleep(100);
					fp = fopen("spline.dat", "r");
					for (i = 0; i != prop.NX; ++i)
					{
						fscanf(fp, "%lf%lf", &non, &Vt_total[i]);
					}
					fclose(fp);
					///////////////////////////////////////////////////////////////////////
					fp = fopen("input.dat", "w");
					for (i = 0; i != bfkin.plane_r_numb; ++i)
					{
						fprintf(fp, "%10.6lf %10.6lf\n", r_ave[i] / bfkin.R, Vr_ave[i]);
					}
					fclose(fp);

					fp = fopen("SPLINE.KIN", "w");
					fprintf(fp, "%d %d\ninput.dat\noutputx.dat", bfkin.plane_r_numb, prop.NX);
					fclose(fp);

					sprintf(cmd, "SPLINE.exe");
					outexe(cmd);
					Sleep(100);
					fp = fopen("spline.dat", "r");
					for (i = 0; i != prop.NX; ++i)
					{
						fscanf(fp, "%lf%lf", &non, &Vr_total[i]);
					}
					fclose(fp);

					debug("Induced Velocity Spline!", getpid());


					/*==================================================================================================
					* ==============================計算誘導速度並提取(徑像速度忽略不計)===================================
					* ==================================================================================================
					*/
					if (ite - bfkin.itestar == 0)
					{
						//第一次放入patpans11時，設誘導速度=0
						for (i = 0; i != prop.NX; ++i)
						{
							Ux_bem[i] = 0.0;
							Ut_bem[i] = 0.0;
							Ur_bem[i] = 0.0;
						}
						debug("First Patpans11 Induced Velocity = 0!", getpid());
					}
					else //除第一輪外誘導速度以SFPV計算
					{
						//讀取前一輪之hXXX，以sfpv11.exe計算誘導速度
						sprintf(filename, "h%d", ite - bfkin.ci);
						fp = fopen("sfpv11.kin", "w");
						fprintf(fp, filename);
						fprintf(fp, "\n0");
						fclose(fp);
						Sleep(100);
						sprintf(cmd, "cmd.exe /c sfpv11.exe<sfpv11.kin");
						outexe(cmd);
						Sleep(100);
						fp = fopen("indvel0.dat", "r");
						fgets(str, 50, fp);
						fgets(str, 50, fp);
						fgets(str, 60, fp);
						for (i = 0; i != prop.MR; ++i)fscanf(fp, "%lf%lf%lf%lf%lf%lf", &r_U[i], &Ux[i], &Ur[i], &Ut[i], &non, &non);
						fclose(fp);

						debug("sfpv calcu Induced Velocity Load!", getpid());
						//將各徑向位置(MR)之誘導速度，藉由SPLINE.exe化為 NX 個點，點位置由XR:radii of the input geometry parameters for patpans11 決定
						//Ux
						fp = fopen("input.dat", "w");
						for (i = 0; i != prop.MR; ++i)fprintf(fp, "%10.6lf %10.6lf\n", r_U[i], Ux[i]);
						fclose(fp);

						fp = fopen("SPLINE.KIN", "w");
						fprintf(fp, "%d %d\ninput.dat\noutputx.dat", prop.MR, prop.NX);
						fclose(fp);

						sprintf(cmd, "SPLINE.exe");
						outexe(cmd);
						Sleep(100);
						fp = fopen("spline.dat", "r");
						for (i = 0; i != prop.NX; ++i)fscanf(fp, "%lf%lf", &non, &Ux_bem[i]);
						fclose(fp);

						//Vt
						fp = fopen("input.dat", "w");
						for (i = 0; i != prop.MR; ++i)fprintf(fp, "%10.6lf %10.6lf\n", r_U[i], Ut[i]);
						fclose(fp);

						fp = fopen("SPLINE.KIN", "w");
						fprintf(fp, "%d %d\ninput.dat\noutputx.dat", prop.MR, prop.NX);
						fclose(fp);

						sprintf(cmd, "SPLINE.exe");
						outexe(cmd);
						Sleep(100);
						fp = fopen("spline.dat", "r");
						for (i = 0; i != prop.NX; ++i)fscanf(fp, "%lf%lf", &non, &Ut_bem[i]);
						fclose(fp);

						//Vr
						fp = fopen("input.dat", "w");
						for (i = 0; i != prop.MR; ++i)fprintf(fp, "%10.6lf %10.6lf\n", r_U[i], Ur[i]);
						fclose(fp);

						fp = fopen("SPLINE.KIN", "w");
						fprintf(fp, "%d %d\ninput.dat\noutputx.dat", prop.MR, prop.NX);
						fclose(fp);

						sprintf(cmd, "SPLINE.exe");
						outexe(cmd);
						Sleep(100);
						fp = fopen("spline.dat", "r");
						for (i = 0; i != prop.NX; ++i)fscanf(fp, "%lf%lf", &non, &Ur_bem[i]);
						fclose(fp);

						debug("sfpv Induced Velocity Spline!", getpid());

					}
					/*==================================================================================================
					* bodyforce_offon = 1為定J計算
					* ==================================================================================================
					*/
					if (bfkin.bodyforce_offon == 1)
					{
						//提取每一輪之前一步的drag
						prop.ADVCO = bfkin.firstJ;
						sprintf(filename, "XYZ_Internal_Table_table_%d.csv", ite);
						fp2 = fopen(filename, "r");
						fgets(str, 80, fp2);
						fscanf(fp2, "%lf,%lf,%lf,%lf", &HULLDRAG, &non, &non, &non);
						fclose(fp2);
						debug("Fixed J -- HULLDRAG load!", getpid());

					}
					else
					{
						/*========================================================================================================
						* 建立patpans11用GEO,ADM檔，先算一次patpans11，所得之hXXX.oup用以得到1-w，之後用於open water test 找kt/j^2 交點
						* ========================================================================================================
						*/
						//建立hXXX.geo，ADVCO採用prop.geo原值，並除入流與船速比外，完全比照prop.geo
						Output_PROP_GEO(&prop, &bfkin);

						//建立hXXX.adm，完全比照prop.adm
						Output_ite_ADM();
						
						//建立patpans11.kin
						sprintf(filename, "h%d", ite);
						fp = fopen("patpans11.kin", "w");
						fprintf(fp, filename);
						fclose(fp);

						//先執行一次patpans11.exe
						sprintf(cmd, "cmd.exe /c patpans11.exe<patpans11.kin");
						outexe(cmd);
						Sleep(100);
						debug("Patpans1  Done!", getpid());

						/*========================================================================================================
						* 計算.geo檔所用J值
						* 第一輪及第二輪以二分法計算[Kt/j^2(常數)]與[螺槳K-Jchart]交點之J值，第三輪後以前兩次J、KT值利用牛頓法計算J值						
						* ========================================================================================================
						*/
						

						if ((ite - (bfkin.itestar - bfkin.ci)) / bfkin.ci == 1)//第一輪
						{
							debug("Round1", getpid());
							if (bfkin.firstJ_offon == 1) 
							{
								prop.ADVCO = bfkin.firstJ;
								debug("First J Done", getpid());
							}
							else 
							{
								//提取先執行一次patpans11所得之hXXX.oup裡的1-w值
								sprintf(filename, "Kt-File.txt");
								fp = fopen(filename, "w");
								fprintf(fp, "h%d.oup", ite);
								fclose(fp);

								sprintf(cmd, "GetKt.exe");
								outexe(cmd);
								Sleep(100);
								fp = fopen("W.dat", "r");
								fscanf(fp, "%lf", &W);
								fclose(fp);

								//計算Ktoj
								Ktoj = (HULLDRAG  + bfkin.wr + bfkin.ar - bfkin.sfc) / (bfkin.rho * pow(bfkin.Vs, 2) * pow(2 * bfkin.R, 2) * pow(W, 2));

								//利用Ktoj，執行KTo.exe，計算ja
								sprintf(filename, "Ktoj.dat");
								fp = fopen(filename, "w");
								fprintf(fp, "%lf", Ktoj);
								fclose(fp);

								sprintf(cmd, "KTo.exe");
								outexe(cmd);
								Sleep(100);
								fp = fopen("ja.dat", "r");
								fscanf(fp, "%lf", &ja);
								fclose(fp);

								//計算J值
								prop.ADVCO = ja / W;
								debug("J Done", getpid());
							}
							debug("第一輪以二分法計算J值或自訂義J值  Done!", getpid());
						}
						else if ((ite - (bfkin.itestar - bfkin.ci)) / bfkin.ci == 2)
						{
							debug("Round2", getpid());
							//提取前一輪hXXX.geo之J值
							sprintf(filename, "h%d.geo", ite - bfkin.ci);
							fp = fopen(filename, "r");
							fgets(str, 80, fp);
							fgets(str, 80, fp);
							fgets(str, 80, fp);
							fscanf(fp, "%lf", &prop.ADVCO);
							fclose(fp);

							//提取前一輪之hXXX.oup裡的CD=0.0035之KT值
							sprintf(filename, "Kt-File.txt");
							fp = fopen(filename, "w");
							fprintf(fp, "h%d.oup", ite - bfkin.ci);
							fclose(fp);

							sprintf(cmd, "GetKt.exe");
							outexe(cmd);
							Sleep(100);
							fp = fopen("KT.dat", "r");
							fscanf(fp, "%lf", &KT);
							fclose(fp);
							debug("*.KT output!", getpid());

							//計算前一輪THRUST
							THRUST = KT * bfkin.rho * pow(bfkin.Vs, 2) * pow(2 * bfkin.R, 2) / pow(prop.ADVCO, 2);

							//輸出此前一輪之Error
							Error = (HULLDRAG + bfkin.wr + bfkin.ar - bfkin.sfc) - THRUST;
							sprintf(filename, "Error%d.dat", ite - bfkin.ci);
							fp3 = fopen(filename, "w");
							fprintf(fp3, "%lf\n", Error);
							fclose(fp3);
							debug("Error output Done!", getpid());


							//提取先執行一次patpans11所得之hXXX.oup裡的1-w值
							sprintf(filename, "Kt-File.txt");
							fp = fopen(filename, "w");
							fprintf(fp, "h%d.oup", ite);
							fclose(fp);

							sprintf(cmd, "GetKt.exe");
							outexe(cmd);
							Sleep(100);
							fp = fopen("W.dat", "r");
							fscanf(fp, "%lf", &W);
							fclose(fp);

							//計算Ktoj
							Ktoj = (HULLDRAG + bfkin.wr + bfkin.ar - bfkin.sfc) / (bfkin.rho * pow(bfkin.Vs, 2) * pow(2 * bfkin.R, 2) * pow(W, 2));

							//利用Ktoj，執行KTo.exe，計算ja
							sprintf(filename, "Ktoj.dat");
							fp = fopen(filename, "w");
							fprintf(fp, "%lf", Ktoj);
							fclose(fp);

							sprintf(cmd, "KTo.exe");
							outexe(cmd);
							Sleep(100);
							fp = fopen("ja.dat", "r");
							fscanf(fp, "%lf", &ja);
							fclose(fp);

							//計算J值
							prop.ADVCO = ja / W;

							debug("第二輪以二分法計算J值  Done!", getpid());
						}
						else if ((ite - (bfkin.itestar - bfkin.ci)) / bfkin.ci > 2)
						{
							debug("Round3", getpid());
							//提取前一輪hXXX.geo之J值
							sprintf(filename, "h%d.geo", ite - bfkin.ci);
							fp = fopen(filename, "r");
							fgets(str, 80, fp);
							fgets(str, 80, fp);
							fgets(str, 80, fp);
							fscanf(fp, "%lf", &prop.ADVCO);
							fclose(fp);

							//提取前二輪hXXX.geo之J值
							sprintf(filename, "h%d.geo", ite - 2 * bfkin.ci);
							fp = fopen(filename, "r");
							fgets(str, 80, fp);
							fgets(str, 80, fp);
							fgets(str, 80, fp);
							fscanf(fp, "%lf", &ADVCO0);
							fclose(fp);

							//提取前一輪之hXXX.oup裡的CD=0.0035之KT值
							sprintf(filename, "Kt-File.txt");
							fp = fopen(filename, "w");
							fprintf(fp, "h%d.oup", ite - bfkin.ci);
							fclose(fp);

							sprintf(cmd, "GetKt.exe");
							outexe(cmd);
							Sleep(100);
							fp = fopen("KT.dat", "r");
							fscanf(fp, "%lf", &KT);
							fclose(fp);
							debug("*.KT output!", getpid());

							//計算前一輪THRUST
							THRUST = KT * bfkin.rho * pow(bfkin.Vs, 2) * pow(2 * bfkin.R, 2) / pow(prop.ADVCO, 2);

							//輸出此前一輪之Error
							Error = (HULLDRAG + bfkin.wr + bfkin.ar - bfkin.sfc) - THRUST;
							sprintf(filename, "Error%d.dat", ite - bfkin.ci);
							fp3 = fopen(filename, "w");
							fprintf(fp3, "%lf\n", Error);
							fclose(fp3);
							debug("Err output Done!", getpid());
							Err = Error;


							//提取前二輪之hXXX.oup裡的CD=0.0035之KT值
							sprintf(filename, "Kt-File.txt");
							fp = fopen(filename, "w");
							fprintf(fp, "h%d.oup", ite - 2 * bfkin.ci);
							fclose(fp);

							sprintf(cmd, "GetKt.exe");
							outexe(cmd);
							Sleep(100);
							fp = fopen("KT.dat", "r");
							fscanf(fp, "%lf", &KT0);
							fclose(fp);
							debug("*.KT0 output!", getpid());

							//計算前二輪THRUST
							THRUST0 = KT0 * bfkin.rho * pow(bfkin.Vs, 2) * pow(2 * bfkin.R, 2) /pow(ADVCO0, 2);

							//提取前二輪之Error
							sprintf(filename, "Error%d.dat", ite - 2 * bfkin.ci);
							fp = fopen(filename, "r");
							fscanf(fp, "%lf", &Err0);
							fclose(fp);
							debug("Err0 load Done!", getpid());

							//牛頓法計算J值
							if (THRUST == THRUST0) {
								prop.ADVCO = ADVCO0;
							}
							else 
							{
								prop.ADVCO = ((HULLDRAG + bfkin.wr + bfkin.ar - bfkin.sfc) - THRUST) / (Err0 - Err) * (prop.ADVCO - ADVCO0) + prop.ADVCO;
							}

							debug("第三輪以牛頓法計算J值計算J值  Done!", getpid());
						}

					}
					//建立hXXX.geo，ADVCO採用--適才之計算值--，並除入流與船速比外，完全比照prop.geo
					Output_PROP_GEO(&prop, &bfkin);

					//建立hXXX.adm，完全比照prop.adm
					Output_ite_ADM();

					/*========================================================================================================
					*                                   輸出周向平均速度、誘導速度、兩相減後速度
					* ========================================================================================================
					*/
					sprintf(filename, "velocity_%d.txt", ite - bfkin.ci);
					fp = fopen(filename, "w");
					fprintf(fp, "r/R  Vx_total  Vr_total  Vt_total  Vx_induce  Vr_induce  Vt_induce  Vx_total-induce  Vr_total-induce  Vt_total-induce\n");
					for (i = 0; i != prop.NX; ++i)fprintf(fp, "%10.6lf %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf\n", prop.NPARAMETER[0][i], Vx_total[i] / bfkin.Vs, Vr_total[i] / bfkin.Vs, Vt_total[i] / bfkin.Vs, Ux_bem[i], Ur_bem[i], Ut_bem[i], Vx_total[i] / bfkin.Vs - Ux_bem[i], Vr_total[i] / bfkin.Vs - Ur_bem[i], Vt_total[i] / bfkin.Vs - Ut_bem[i]);
					fprintf(fp, "%10d", size);
					fclose(fp);

					debug("velocity_XX.txt output!", getpid());
					Sleep(100);
					/*========================================================================================================
					*                                   執行第二次patpans11
					* ========================================================================================================
					*/
					sprintf(filename, "h%d", ite);
					fp = fopen("patpans11.kin", "w");
					fprintf(fp, filename);
					fclose(fp);

					sprintf(cmd, "cmd.exe /c patpans11.exe<patpans11.kin");
					outexe(cmd);
					Sleep(100);
					debug("Patpans2 Done!", getpid());

					//===================================提取hXXX.oup裡的CD=0.0035之KT值，用於計算THRUST========================		
					sprintf(filename, "Kt-File.txt");
					fp = fopen(filename, "w");
					fprintf(fp, "h%d.oup", ite);
					fclose(fp);

					sprintf(cmd, "GetKt.exe");
					outexe(cmd);
					Sleep(100);
					fp = fopen("KT.dat", "r");
					fscanf(fp, "%lf", &KT);
					fclose(fp);

					//====================================================計算THRUST==========================================
					THRUST = KT * bfkin.rho * pow(bfkin.Vs, 2) * pow(2 * bfkin.R, 2) / pow(prop.ADVCO, 2);

					//================================================輸出hulldragXXX.dat=====================================
					sprintf(filename, "hulldrag%d.dat", ite);
					fp = fopen(filename, "w");
					fprintf(fp, "ite=%d\n", ite);
					fprintf(fp, "KT=%lf\n", KT);
					fprintf(fp, "THRUST=%lf\n", THRUST);
					fprintf(fp, "HULLDRAG=%lf\n", HULLDRAG);
					fprintf(fp, "Wave Resistance=%lf\n", bfkin.wr);
					fprintf(fp, "Added Resistance=%lf\n", bfkin.ar);
					fprintf(fp, "SFC=%lf\n", bfkin.sfc);
					fprintf(fp, "J=%lf\n", prop.ADVCO);
					fclose(fp);

					//=================================================更新hulldrag_ite.dat====================================
					fp = fopen("hulldrag_ite.txt", "w");
					fprintf(fp, "%d\n", bfkin.itestar);
					fprintf(fp, "%d\n", ite);
					fclose(fp);

					//讀取受力
					fp = fopen("bf.dat", "r");
					for (i = 0; i != prop.MR + 4; ++i)fgets(str, 80, fp);
					for (i = 0; i != prop.MR; ++i)fscanf(fp, "%lf%lf%lf%lf%lf%lf%lf", &r_F[i], &fx[i], &fr[i], &ft[i], &non, &non, &non);//r/R, F_X, F_Y, F_Z, Q_X, Q_Y, Q_Z
					fclose(fp);
					
					for (i = 0; i != prop.MR; ++i) 
					{
						bfx[i] = fx[i] / (r_x[i + 1] - r_x[i]) * 0.5 * bfkin.rho * bfkin.Vs * bfkin.Vs * bfkin.R;
						bfr[i] = fr[i] / (r_x[i + 1] - r_x[i]) * 0.5 * bfkin.rho * bfkin.Vs * bfkin.Vs * bfkin.R;
						bft[i] = ft[i] / (r_x[i + 1] - r_x[i]) * 0.5 * bfkin.rho * bfkin.Vs * bfkin.Vs * bfkin.R;
					}

					//=====================================================輸出受力=============================================
					sprintf(filename, "force_%d.txt", ite);
					fp = fopen(filename, "w");
					for (i = 0; i != prop.MR; ++i) {
						fprintf(fp, "%lf %lf %lf %lf\n", r_F[i] * bfkin.R, bfx[i], bfr[i], bft[i]);
					}
					fclose(fp);
					/*========================================================================================================
					*											執行chebyToCn.exe取得Cn
					* ========================================================================================================
					*/
					//bfx - 執行chebyToCn.exe取得Cn
					sprintf(filename, "fx%d.dat", ite);
					fp = fopen(filename, "w");
					for (i = 0; i != prop.MR; ++i) { fprintf(fp, "%lf %lf\n", r_F[i] * bfkin.R, bfx[i]); }
					fclose(fp);
					sprintf(filename, "fx.dat");
					fp = fopen(filename, "w");
					for (i = 0; i != prop.MR; ++i) { fprintf(fp, "%f %f\n", r_F[i] * bfkin.R, bfx[i]); }
					fclose(fp);

					fp = fopen("chebyToCn.kin", "w");
					sprintf(str, "fx\n6");
					fprintf(fp, str);
					fclose(fp);

					sprintf(cmd, "chebyToCn.exe");//傅立葉轉換
					outexe(cmd);

					Sleep(100);

					//bfr - 執行chebyToCn.exe取得Cn
					sprintf(filename, "fr%d.dat", ite);
					fp = fopen(filename, "w");
					for (i = 0; i != prop.MR; ++i) { fprintf(fp, "%lf %lf\n", r_F[i] * bfkin.R, bfr[i]); }
					fclose(fp);
					sprintf(filename, "fr.dat");
					fp = fopen(filename, "w");
					for (i = 0; i != prop.MR; ++i) { fprintf(fp, "%lf %lf\n", r_F[i] * bfkin.R, bfr[i]); }
					fclose(fp);

					fp = fopen("chebyToCn.kin", "w");
					sprintf(str, "fr\n6");
					fprintf(fp, str);
					fclose(fp);

					sprintf(cmd, "chebyToCn.exe");
					outexe(cmd);
					Sleep(100);
					//bft - 執行chebyToCn.exe取得Cn
					sprintf(filename, "ft%d.dat", ite);
					fp = fopen(filename, "w");
					for (i = 0; i != prop.MR; ++i) { fprintf(fp, "%lf %lf\n", r_F[i] * bfkin.R, bft[i]); }
					fclose(fp);
					sprintf(filename, "ft.dat");
					fp = fopen(filename, "w");
					for (i = 0; i != prop.MR; ++i) { fprintf(fp, "%lf %lf\n", r_F[i] * bfkin.R, bft[i]); }
					fclose(fp);

					fp = fopen("chebyToCn.kin", "w");
					sprintf(str, "ft\n6");
					fprintf(fp, str);
					fclose(fp);

					sprintf(cmd, "chebyToCn.exe");
					outexe(cmd);
					Sleep(100);
					debug("ChebyToCn Done!", getpid());
					










				}



			}

			if (ite > bfkin.itestar && bfkin.Unsteady_offon == 0)
			{
				debug("00000", getpid());
			}



		}


		




	}
}
