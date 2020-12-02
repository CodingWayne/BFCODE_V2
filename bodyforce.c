
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
#include <io.h>


//[result]:the array of values returned by the user function.(Force/Volume);  [size]:the number of elements in the result array.(int);  [centroid]:���椤���I���y��;  [Velocity]:���椤���I���t��; [Iteration]:the iteration number in CoordReal precision.(int); [Angle]:����pitch���|��(�F������); [Cell_id]:��L���椧�s��; [Drag]:��ߪ��O;
void USERFUNCTION_EXPORT bodyforce(Real(*result)[3], int size, CoordReal(*centroid)[3], CoordReal(*Velocity)[3], CoordReal(*Cell_id)[1], CoordReal(*Iteration)[1], CoordReal(*Angle)[1])
{

	if (size != 0)
	{
		struct BFKIN bfkin;
		struct PROP prop;

		
		pitch = 1;
		bfkin.R = 1;
		//Output_bfkin(&bfkin, &reitestar, &pitch);
		PI = (double)acos(-1.0);
		ite = (int)Iteration[0][0];
		pitch = (double)Angle[0][0]; //angle[0] ��U��pitch�|��
		
		if (ite == 0)
		{
			if (access("debug", 0))system("mkdir debug");
			if (access("thread", 0))system("mkdir thread");
			if (access("mesh", 0))system("mkdir mesh");
			if (access("tmp", 0))system("mkdir tmp");
		}

		Load_bfkin(&bfkin);
		reitestar = bfkin.itestar % bfkin.ci;
		Output_bfkin(&bfkin, &reitestar, &pitch);

		sprintf(str, ".\\debug\\debug%d.txt", getpid());
		fp = fopen(str, "w");
		fprintf(fp, "=================  Iteration=%d  =================\n", ite);
		fclose(fp);
		

		//=======================================================================����Ѽƪ�l��======================
		meshload = 0; minCellid = 999999999, np = 0, threadboy = 0, * newsize, * threadnumber, all_size = 0, non = 0;
		//==========================================================================================================

		/*==========================================================================================================
		*                                       �C�@�B���e�m�@�~      ��J���O
		* ==========================================================================================================
		*/
		if (ite > bfkin.itestar )
		{
			debug("Steady Body Force Star!", getpid());
			/*================================================================================================
			* ��mesh�x�} ���mesh.txt��x,y,z, , ,radian,theta
			* ================================================================================================
			*/
			mesh = (double**)malloc((bfkin.plane_x_numb * 4 * bfkin.plane_c_numb * bfkin.plane_r_numb) * sizeof(double*));
			for (i = 0; i != (bfkin.plane_x_numb * 4 * bfkin.plane_c_numb * bfkin.plane_r_numb); ++i)
			{
				mesh[i] = (double*)malloc(8 * sizeof(double));
			}

			sprintf(str, ".\\mesh\\mesh.txt");
			fp = fopen(str, "r");

			for (i = 0; i < (bfkin.plane_x_numb * 4 * bfkin.plane_c_numb * bfkin.plane_r_numb); i++)
			{
				fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf", &mesh[i][0], &mesh[i][1], &mesh[i][2], &non, &non, &non, &mesh[i][6], &mesh[i][7]);

			}
			fclose(fp);
			debug("mesh load!", getpid());
			/*================================================================================================
			* Ū�Jprop.geo
			* ================================================================================================
			*/
			sprintf(geoname, "prop.geo");
			fp = fopen(geoname, "r");
			fgets(label, 100, fp);
			fscanf(fp, "%d%d%d%d%d", &prop.NX, &prop.NBLADE, &prop.NC, &prop.MR, &prop.NTMP);
			fscanf(fp, "%lf%d%d%lf%lf%lf", &prop.RHUB, &prop.NHBU, &prop.MHBT, &prop.XHBU, &prop.XHBD, &prop.XHBT);
			fscanf(fp, "%lf%lf%lf%lf%lf%lf%lf", &prop.ADVCO, &prop.RULT, &prop.RHULT, &prop.DCD, &prop.XULT, &prop.DTPROP, &prop.XUWDK);
			prop.NPARAMETER = (double**)malloc(14 * sizeof(double*)); //�ΥH�x�sprop.geo�ĥ|�椧�᪺�ƭ� 
			for (i = 0; i != 14; i++)prop.NPARAMETER[i] = (double*)malloc(prop.NX * sizeof(double));
			for (j = 0; j != 14; ++j)for (i = 0; i != prop.NX; ++i)fscanf(fp, "%lf", &prop.NPARAMETER[j][i]);
			fclose(fp);

			debug("prop.geo load!", getpid());

			//r_x:�ΥH�x�s�������V���q(MR�Ӭq)�̡A�]�t�Y�����C�I��r/R (���V���q��sine spanwise) 
			r_x = (double*)malloc((prop.MR + 1) * sizeof(double));

			for (i = 0; i != prop.MR + 1; ++i) 
			{
				r_x[i] = prop.RHUB + (1 - prop.RHUB) * (double)sin(0.5 * PI / prop.MR * i);
			}

			
			//��n�O���|�V���G�Ƭ�6�ӫY��(C0~C5)�����ܼƨ��(T0~T5)�����A�ܼƬ��|�V��mr
			//Ū��Cn��(C0~C5)�A�ƭȤw��chebyToCn.exe��X
			fp = fopen("fxCn.dat", "r");
			for (i = 0; i != 6; ++i)fscanf(fp, "%lf", &Cx[i]);
			fclose(fp);
			fp = fopen("frCn.dat", "r");
			for (i = 0; i != 6; ++i)fscanf(fp, "%lf", &Cr[i]);
			fclose(fp);
			fp = fopen("ftCn.dat", "r");
			for (i = 0; i != 6; ++i)fscanf(fp, "%lf", &Ct[i]);
			fclose(fp);

			debug("Cn load!", getpid());

			/*===================================================================================================================
			*													�C�Ӻ���  ��J���O
			* ===================================================================================================================
			*/
			for (i = 0; i != size; ++i) 
			{
				kk = (int)Cell_id[i][0] - 1;

				//�p�����Ҧb���b�|�j�p
				r = mesh[kk][6];

				//�Q�κ���Ҧb��y�ȡA�p��theta�|��
				theta = mesh[kk][7];


				//�p����ܼƨ��(T0~T5)�A�ܼƬ��|�V��mr
				Tn[0] = 1; //T0(r) = 1
				Tn[1] = (r - (r_x[0] + r_x[1]) / 2.0 * bfkin.R) / ((r_x[prop.MR - 1] + r_x[prop.MR]) / 2.0 * bfkin.R - (r_x[0] + r_x[1]) / 2.0 * bfkin.R) * 2 - 1; //T1(r) = r = (r-(r_x[0]+r_x[1])/2.0*R)/((r_x[MR-1]+r_x[MR])/2.0*R-(r_x[0]+r_x[1])/2.0*R)*2-1
				for (j = 2; j < 6; j++) { Tn[j] = 2 * ((r - (r_x[0] + r_x[1]) / 2.0 * bfkin.R) / ((r_x[prop.MR - 1] + r_x[prop.MR]) / 2.0 * bfkin.R - (r_x[0] + r_x[1]) / 2.0 * bfkin.R) * 2 - 1) * Tn[j - 1] - Tn[j - 2]; } //Tn+1(r) = 2*x*Tn(r) - Tn-1(r) 
				//fprintf(fp1,"%d %f %f %f %f %f %f\n",i,Tn[0],Tn[1],Tn[2],Tn[3],Tn[4],Tn[5]);

				//�p��x��V���ը��O(��n�O) - �qinput��output��V
				sum = 0;
				for (j = 0; j < 6; j++)sum = sum + Tn[j] * Cx[j];
				originfx = sum - 0.5 * Cx[0];

				//�p��r��V���ը��O(��n�O) - �V�~����
				sum = 0;
				for (j = 0; j < 6; j++)sum = sum + Tn[j] * Cr[j];
				originfr = sum - 0.5 * Cr[0];

				//�p��t��V���ը��O(��n�O) - �H���V�J�y���譱�A�q��Z�b�f�ɰw�����A���ɥ�Y�¥k
				sum = 0;
				for (j = 0; j < 6; j++)sum = sum + Tn[j] * Ct[j];
				originft = -(sum - 0.5 * Ct[0]);

				//fprintf(fp2,"%d %f %f %f\n",i,originfx,originfr,originft);


				//�p������n�O
				result[i][0] = (originfx * prop.NBLADE / 2.0 / PI / r / bfkin.x) * cos(pitch) + ((originfr * cos(theta) - originft * sin(theta)) * prop.NBLADE / 2.0 / PI / r / bfkin.x) * sin(pitch); //�p��x��V�����n�O ,NBLADE:������
				result[i][1] = (-originfr * sin(theta) - originft * cos(theta)) * prop.NBLADE / 2.0 / PI / r / bfkin.x;														 //�p��y��V�����n�O ,NBLADE:������	


				if (bfkin.bow_dir == 1)
				{
					result[i][0] = result[i][0] * -1;
					result[i][1] = result[i][1] * -1;
				}

				result[i][2] = -(originfx * prop.NBLADE / 2.0 / PI / r / bfkin.x) * sin(pitch) + ((originfr * cos(theta) - originft * sin(theta)) * prop.NBLADE / 2.0 / PI / r / bfkin.x) * cos(pitch);//�p��z��V�����n�O ,NBLADE:������

			}


			debug("result Steady Force input!", getpid());

			free(r_x);

			if (ite % bfkin.ci != reitestar)
			{
				free(prop.NPARAMETER);//�b������B�J���Τ���

				if (ite == 0) //��if�Τ���A��0�B���|�i"��J���O"�A�|�����i"�@�뱡�p"
				{
					if (threadboy == 1)
					{
						free(newsize); free(threadnumber);

						for (i = 0; i != (bfkin.plane_x_numb * 4 * bfkin.plane_c_numb * bfkin.plane_r_numb); ++i) 
						{
							free(mesh[i]);
						}
						free(mesh);
					}
				}

				if (meshload == 0)
				{
					if (threadboy == 0) 
					{
						//�M���UProcess�����x�}
						for (i = 0; i != (bfkin.plane_x_numb * 4 * bfkin.plane_c_numb * bfkin.plane_r_numb); ++i)
						{
							free(mesh[i]);
						}
						free(mesh);
					}
				}

				//�A�M���@���UProcess�����x�}
				threadboy = 1;
				if (ite >= bfkin.itestar && ite % bfkin.ci == reitestar && meshload == 0) 
				{
					if (threadboy == 1) 
					{
						free(newsize); free(threadnumber);

						for (i = 0; i != (bfkin.plane_x_numb * 4 * bfkin.plane_c_numb * bfkin.plane_r_numb); ++i) 
						{
							free(mesh[i]);
						}
						free(mesh);
					}
				}

			}
			//���n! threadboy�n�k0�H�K����p��
			threadboy = 0;
			debug("[][] Clear!", getpid());
		}
		/*==========================================================================================================
		*                                                     �@�뱡�p
		* ==========================================================================================================
		*/
		if (ite == 0 || ite >= bfkin.itestar && ite % bfkin.ci == reitestar && Iteration[0] != 0)
		{
			
			threadboy = 0;
			Output_bfkin(&bfkin, &reitestar, &pitch);
			//=======================================================================������L����========================
			if (ite == 0)
			{
				meshload = 1;//meshload=1 ��ܤw�������
				/*
				 *STEP[1]
				 *�ϥ�getpid��X�Uprocess�W�۳B�z��������(thread_pid_ite.txt)
				 *���]��5��process�b���u�B�z�Athread��Ƨ��K�|��X5��(thread_pid_ite.txt)����5�Ӥ��P��pid
				 *(thread_pid_ite.txt)=[����s��,���椤���IX,���椤���Iy,���椤���Iz,���椤���I�t��X,���椤���I�t��y,���椤���I�t��z]
				 *����minCellid�b�U�O��process�̰O���̤p������s�� (minCellid���O�x�}�u�O���ܼơA���b�U��process�̥N���process��minCellid)
				 */				
				sprintf(str, ".\\thread\\%d_%d.txt", getpid(), ite);
				fp = fopen(str, "w");
				debug("STEP[1]-1   open thread file", getpid());
				fprintf(fp, "%d\n", size);
				for (i = 0; i < size; i++)
				{
					fprintf(fp, "%lf %lf %lf %lf %lf %lf %lf\n", Cell_id[i][0], centroid[i][0], centroid[i][1], centroid[i][2], Velocity[i][0], Velocity[i][1], Velocity[i][2]);
					if (minCellid > (int)Cell_id[i][0])
					{
						minCellid = (int)Cell_id[i][0];
					}
				}
				fclose(fp);		

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
				//=======================================================================Body Force�}�l=====================
				/*���B�J���C��Bodyforce�}�l�e���n�B�J
				*�b�C��Body force �}�l���e�A1.��X�ӽ����L����W���t�� 2.�H��֤߳B�z��s���t�ק�s��mesh[][]��
				*�b�ӽ���Body force������ �|����mesh[][],�H��threadboy�k0
				*/
				if (ite >= bfkin.itestar && ite % bfkin.ci == reitestar && meshload == 0)
				{
					//output thread number & point site
					//�̦U�֤ߧǿ�X�]�t�֤ߤ�size,centroid[][],Velocity[][]���
					sprintf(str, ".\\thread\\%d_%d.txt", getpid(), ite); //��X�U�֤ߤ�size,centroid[][],Velocity[][]
					fp = fopen(str, "w");
					fprintf(fp, "%d\n", size);
					for (i = 0; i < size; i++) 
					{
						fprintf(fp, "%lf %lf %lf %lf %lf %lf %lf\n", Cell_id[i][0], centroid[i][0], centroid[i][1], centroid[i][2], Velocity[i][0], Velocity[i][1], Velocity[i][2]);
					}
					fclose(fp);

					Sleep(5000);

					MESH_STEP_3();//��np
					MESH_STEP_4();//�H��֤��x�sthread�̪�mesh�x�}

					Sleep(5000);

					if (threadboy == 1)
					{
						MESH_STEP_5(&bfkin);//��mesh�x�}
						//Ū��minCellid
						fp = fopen(".\\mesh\\minCellid.txt", "r");
						fscanf(fp, "%d", &minCellid);
						fclose(fp);
						//Ū��mesh.txt�A�ç�s������t�׸��
						sprintf(str, ".\\mesh\\mesh.txt");
						fp = fopen(str, "r");
						for (i = 0; i < (bfkin.plane_x_numb * 4 * bfkin.plane_c_numb * bfkin.plane_r_numb); i++)
						{
							fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf", &mesh[i][0], &mesh[i][1], &mesh[i][2], &non, &non, &non, &mesh[i][6], &mesh[i][7]);
						}
						fclose(fp);
						debug("threadboy=1 mesh.txt load Done!", getpid());
						//�qthread��Ƨ��̧���s���t�ר�mesh�x�}��
						for (i = 0; i < np; i++)
						{
							sprintf(str, ".\\thread\\%d_%d.txt", threadnumber[i], ite);
							fp = fopen(str, "r");

							debug(str, getpid());
							fp3 = fopen("IMD.dat", "w+");
							fgets(str, 200, fp);
							for (j = 0; j < newsize[i]; j++) 
							{
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
						debug("threadboy=1 mesh load done!", getpid());
					}
				}
				//
				//=======================================��L�֤߫ظmmesh�x�}Ū��mesh.txt(�{�u����֤�threadboy=1��mesh�x�}�Ӥw)=================================
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
						debug("threadboy=0 mesh build done!", getpid());
						//Ū��mesh.txt
						sprintf(str, ".\\mesh\\mesh.txt");
						fp = fopen(str, "r");
						debug("threadboy=0 mesh.txt open!", getpid());
						for (i = 0; i < (bfkin.plane_x_numb * 4 * bfkin.plane_c_numb * bfkin.plane_r_numb); i++)
						{
							fscanf(fp, "%lf%lf%lf%lf%lf%lf%lf%lf", &mesh[i][0], &mesh[i][1], &mesh[i][2], &non, &non, &non, &mesh[i][6], &mesh[i][7]);
						}
						fclose(fp);
						debug("threadboy=0 mesh[][] load done!", getpid());
					}
				}
				if (ite > 0)
				{
					fp = fopen(".\\mesh\\minCellid.txt", "r");
					fscanf(fp, "%d", &minCellid);
					fclose(fp);
				}
				//============================================================Ū��prop.geo========================================
				Load_PROP_GEO(&prop);//Output_PROP_GEO(&prop);			
				//============================================================Ū��prop.adm========================================
				Load_PROP_ADM();//CHECK_PROP_ADM();
				//============================================================�t�m����============================================
				Vx_ave = (double*)malloc(bfkin.plane_r_numb * sizeof(double));
				Vt_ave = (double*)malloc(bfkin.plane_r_numb * sizeof(double));//Vx_ave,Vt_ave,Vr_ave:�ΥH�x�s�U�|�V��m���檺�P�V�����t��
				Vr_ave = (double*)malloc(bfkin.plane_r_numb * sizeof(double));
				Vx_sum = (double*)malloc(bfkin.plane_r_numb * sizeof(double));//Vx_sum,Vt_sum,Vr_sum:�ΥH�x�s�U�|�V��m���檺�P�V�t���`�M
				Vt_sum = (double*)malloc(bfkin.plane_r_numb * sizeof(double));
				Vr_sum = (double*)malloc(bfkin.plane_r_numb * sizeof(double));
				r_ave = (double*)malloc(bfkin.plane_r_numb * sizeof(double)); // r_ave:�ΥH�x�s�U�|�V��m���檺�����b�|�j�p
				r_sum = (double*)malloc(bfkin.plane_r_numb * sizeof(double)); // r_sum:�ΥH�x�s�U�|�V��m���檺�b�|�j�p�`�M
				n = (double*)malloc(bfkin.plane_r_numb * sizeof(double));     // n:�ΥH�x�s�U�|�V��m���P�V���檺�ƶq
				//�ϥ�malloc()�禡�t�m�ʺA���Ŷ� plane_r_numb:number of spanwise  mesh in Disk        MR:number of spanwise panels in patpans11

				fx = (double*)malloc(prop.MR * sizeof(double));         //fx[],fr[],ft[]:patpans11.exe�p�⤧�O?
				fr = (double*)malloc(prop.MR * sizeof(double));
				ft = (double*)malloc(prop.MR * sizeof(double));
				bfx = (double*)malloc(prop.MR * sizeof(double));		//bfx[],bfr[],bft[] = fx~r~t[i]/0.5*rho*Vs*Vs*(r_x[i+1]-r_x[i])*R;
				bfr = (double*)malloc(prop.MR * sizeof(double));
				bft = (double*)malloc(prop.MR * sizeof(double));
				r_F = (double*)malloc(prop.MR * sizeof(double));		//r_F[]:�ΥH�x�spatpans11.exe�p�⤧�O���|�V��m�A�ƶq��MR�ӡA(���V���q��sine spanwise) 
				r_x = (double*)malloc((prop.MR + 1) * sizeof(double));  //r_x[]:�ΥH�x�s�������V���q(MR�Ӭq)�̡A�]�t�Y�����C�I��r/R�A�ƶq��MR+1�ӡA(���V���q��sine spanwise) 
				r_U = (double*)malloc(prop.MR * sizeof(double));		//r_U[]:�ΥH�x�ssfpv11.exe�p�⤧���ɳt�ת��|�V��m�A�ƶq��MR�ӡA(���V���q��sine spanwise) 
				Vx_total = (double*)malloc(prop.NX * sizeof(double));   //Vx_total[],Vt_total[],Vr_total[]
				Vt_total = (double*)malloc(prop.NX * sizeof(double));   //�ΥH�x�s�U�|�V��m���椧�P�V�����t��  SPLINE�� NX ���I�᪺���G�A��@��L���߭����J�y�t��
				Vr_total = (double*)malloc(prop.NX * sizeof(double));
				Ux = (double*)malloc(prop.MR * sizeof(double));			//Ux[],Ut[],Ur[]:�ΥH�x�ssfpv11.exe�p�⤧�U�|�V��m�����ɳt�׻P��t�����(Vinduced/Vs)�A�ƶq��MR��
				Ut = (double*)malloc(prop.MR * sizeof(double));
				Ur = (double*)malloc(prop.MR * sizeof(double));
				Ux_bem = (double*)malloc(prop.NX * sizeof(double));		//Ux_bem[],Ut_bem[],Ur_bem[]
				Ut_bem = (double*)malloc(prop.NX * sizeof(double));     //�ΥH�x�ssfpv11.exe�p�⤧�U�|�V��m�����ɳt�׻P��t�����(Vinduced/Vs) 
				Ur_bem = (double*)malloc(prop.NX * sizeof(double));		//SPLINE�� NX ���I�᪺���G�A��@��L���߭������ɳt��
				debug("=====  Pre     Done  =================", getpid());
				//======================================================�p�ⶡ�Z==============================================================
				for (i = 0; i != prop.MR + 1; ++i)
				{
					r_x[i] = prop.RHUB + (1 - prop.RHUB) * (double)sin(0.5 * PI / prop.MR * i);
					//r_x[]:�ΥH�x�s�������V���q(MR�Ӭq)�̡A�]�t�Y�����C�I��r/R�A�ƶq��MR+1�ӡA(���V���q��sine spanwise) 
				}
				//======================================================��n�O�k0=============================================================
				if (ite <= bfkin.itestar)
				{
					//�Ĥ@������patpans11�e�A�ҩ�J�������n�O��0
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
				*********											  �}�l�p����n														    *********
				*********																										(�H��֤߹B�@)	*********
				****************************************************************************************************************************************/
				if (ite >= bfkin.itestar && ite % bfkin.ci == reitestar && Iteration[0] != 0 && threadboy == 1)
				{
					//================================����total_velocity�ð��P�V����(�|�V�t�ש������p)======================================================

					for (i = 0; i != bfkin.plane_r_numb; ++i)
					{
						Vx_sum[i] = 0.0;
						Vt_sum[i] = 0.0;
						Vr_sum[i] = 0.0;
						r_sum[i] = 0.0;
						n[i] = 0.0;
					}
					//�H�Ǵ�������A�p���b�|�j�p,�|�סA�H��|�V��m�@�P�V�[�`
					for (i = 0; i != (bfkin.plane_x_numb * 4 * bfkin.plane_c_numb * bfkin.plane_r_numb); ++i)
					{
						//�p�����Ҧb���b�|�j�p
						r = mesh[i][6];

						//�Q�κ���Ҧb��x�ȡA�p��theta�|��
						theta = mesh[i][7];
						for (j = 0; j != bfkin.plane_r_numb; ++j)
						{
							//�P�_�O�_����L���߽u�W������A�çP�_�O���Ӯ|�V��m������
							//�N�U�|�V��m������A�P�V�[�`��all_Velocity[][0~2],�|�V��mr,�ƶqn
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

					//�N�U�`�M�ȡA�H�U�|�V��m������ƶq���P�V����
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

					//==========================================================�N�U�|�V��m���P�V�����t�סA�ǥ�SPLINE.exe�Ƭ� NX ���I===============================================
					//==========================================================�I��m��XR:radii of the input geometry parameters for patpans11 �M�w==============================
					//==========================================================�Nplane_r_numb�Ӯ|�V���P�V�����t��-----(SPLINE)---->��NX�ӳt��======================================
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
						fprintf(fp, "%8.5lf\n", prop.NPARAMETER[0][i]);// NX ���I��m��XR:radii of the input geometry parameters for patpans11 �M�w
					}
					fclose(fp);

					sprintf(cmd, "SPLINE.exe");
					outexe(cmd);
					Sleep(1000);
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
					Sleep(1000);
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
					Sleep(1000);
					fp = fopen("spline.dat", "r");
					for (i = 0; i != prop.NX; ++i)
					{
						fscanf(fp, "%lf%lf", &non, &Vr_total[i]);
					}
					fclose(fp);

					debug(" Velocity Spline!", getpid());


					/*==================================================================================================
					* ==============================�p�⻤�ɳt�רô���(�|���t�ש������p)===================================
					* ==================================================================================================
					*/
					if (ite - bfkin.itestar == 0)
					{
						//�Ĥ@����Jpatpans11�ɡA�]���ɳt��=0
						for (i = 0; i != prop.NX; ++i)
						{
							Ux_bem[i] = 0.0;
							Ut_bem[i] = 0.0;
							Ur_bem[i] = 0.0;
						}
						debug("First Patpans11 Induced Velocity = 0!", getpid());
					}
					else //���Ĥ@���~���ɳt�ץHSFPV�p��
					{
						//Ū���e�@����hXXX�A�Hsfpv11.exe�p�⻤�ɳt��
						sprintf(filename, "h%d", ite - bfkin.ci);
						fp = fopen("sfpv11.kin", "w");
						fprintf(fp, filename);
						fprintf(fp, "\n0");
						fclose(fp);
						Sleep(1000);
						sprintf(cmd, "cmd.exe /c sfpv11.exe<sfpv11.kin");
						outexe(cmd);
						Sleep(1000);
						fp = fopen("indvel0.dat", "r");
						fgets(str, 50, fp);
						fgets(str, 50, fp);
						fgets(str, 60, fp);
						for (i = 0; i != prop.MR; ++i)fscanf(fp, "%lf%lf%lf%lf%lf%lf", &r_U[i], &Ux[i], &Ur[i], &Ut[i], &non, &non);
						fclose(fp);

						debug("sfpv calcu Induced Velocity Load!", getpid());
						//�N�U�|�V��m(MR)�����ɳt�סA�ǥ�SPLINE.exe�Ƭ� NX ���I�A�I��m��XR:radii of the input geometry parameters for patpans11 �M�w
						//Ux
						fp = fopen("input.dat", "w");
						for (i = 0; i != prop.MR; ++i)fprintf(fp, "%10.6lf %10.6lf\n", r_U[i], Ux[i]);
						fclose(fp);

						fp = fopen("SPLINE.KIN", "w");
						fprintf(fp, "%d %d\ninput.dat\noutputx.dat", prop.MR, prop.NX);
						fclose(fp);

						sprintf(cmd, "SPLINE.exe");
						outexe(cmd);
						Sleep(1500);
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
						Sleep(1500);
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
						Sleep(3000);
						fp = fopen("spline.dat", "r");
						for (i = 0; i != prop.NX; ++i)fscanf(fp, "%lf%lf", &non, &Ur_bem[i]);
						fclose(fp);

						debug("sfpv Induced Velocity Spline!", getpid());

					}
					/*==================================================================================================
					* bodyforce_offon = 1���wJ�p��
					* ==================================================================================================
					*/
					if (bfkin.bodyforce_offon == 1)
					{
						//�����C�@�����e�@�B��drag
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
						* �إ�patpans11��GEO,ADM�ɡA����@��patpans11�A�ұo��hXXX.oup�ΥH�o��1-w�A����Ω�open water test ��kt/j^2 ���I
						* ========================================================================================================
						*/
						//�إ�hXXX.geo�AADVCO�ĥ�prop.geo��ȡA�ð��J�y�P��t��~�A�������prop.geo
						sprintf(filename, "h%d.geo", ite);
						fp = fopen(filename, "w");
						fprintf(fp, label);
						fprintf(fp, "%d %d %d %d %d\n", prop.NX, prop.NBLADE, prop.NC, prop.MR, prop.NTMP);
						fprintf(fp, "%.4lf %d %d %.4lf %.4lf %.4lf\n", prop.RHUB, prop.NHBU, prop.MHBT, prop.XHBU, prop.XHBD, prop.XHBT);
						fprintf(fp, "%.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf\n", prop.ADVCO, prop.RULT, prop.RHULT, prop.DCD, prop.XULT, prop.DTPROP, prop.XUWDK);
						for (j = 0; j != 7; ++j) {
							for (i = 0; i != prop.NX; ++i)fprintf(fp, "%8.5lf ", prop.NPARAMETER[j][i]);
							fprintf(fp, "\n");
						}
						//�J�y�P��t�񬰶�L�J�y/Vs-���ɳt��/Vs
						for (i = 0; i != prop.NX; ++i)fprintf(fp, "%.5lf ", Vx_total[i] / bfkin.Vs - Ux_bem[i]);
						fprintf(fp, "\n");
						for (i = 0; i != prop.NX; ++i)fprintf(fp, "%.5lf ", Vr_total[i] / bfkin.Vs - Ur_bem[i]);
						fprintf(fp, "\n");
						for (i = 0; i != prop.NX; ++i)fprintf(fp, "%.5lf ", Vt_total[i] / bfkin.Vs - Ut_bem[i]);
						fprintf(fp, "\n");
						//hXXX.geo�̫�4��]��0
						for (j = 0; j != 4; ++j) {
							for (i = 0; i != prop.NX; ++i) { fprintf(fp, "0.00000 "); }
							fprintf(fp, "\n");
						}
						fclose(fp);

						debug("OUT�@GEO  Done!", getpid());
						//�إ�hXXX.adm�A�������prop.adm
						Output_ite_ADM();

						//�إ�patpans11.kin
						sprintf(filename, "h%d", ite);
						fp = fopen("patpans11.kin", "w");
						fprintf(fp, filename);
						fclose(fp);

						//������@��patpans11.exe
						sprintf(cmd, "cmd.exe /c patpans11.exe<patpans11.kin");
						outexe(cmd);
						Sleep(2000);
						debug("Patpans1  Done!", getpid());

						/*========================================================================================================
						* �p��.geo�ɩҥ�J��
						* �Ĥ@���βĤG���H�G���k�p��[Kt/j^2(�`��)]�P[����K-Jchart]���I��J�ȡA�ĤT����H�e�⦸J�BKT�ȧQ�Τ��y�k�p��J��
						* ========================================================================================================
						*/
						sprintf(filename, "XYZ_Internal_Table_table_%d.csv", ite);
						fp2 = fopen(filename, "r");
						fgets(str, 80, fp2);
						fscanf(fp2, "%lf,%lf,%lf,%lf", &HULLDRAG, &non, &non, &non);
						fclose(fp2);
						debug("No fixed J HULLDRAG load!", getpid());

						if ((ite - (bfkin.itestar - bfkin.ci)) / bfkin.ci == 1)//�Ĥ@��
						{
							debug("Round1", getpid());
							if (bfkin.firstJ_offon == 1)
							{
								prop.ADVCO = bfkin.firstJ;
								debug("First J Done", getpid());
							}
							else
							{
								//����������@��patpans11�ұo��hXXX.oup�̪�1-w��
								sprintf(filename, "Kt-File.txt");
								fp = fopen(filename, "w");
								fprintf(fp, "h%d.oup", ite);
								fclose(fp);

								sprintf(cmd, "GetKt.exe");
								outexe(cmd);
								Sleep(1500);
								fp = fopen("W.dat", "r");
								fscanf(fp, "%lf", &W);
								fclose(fp);

								//�p��Ktoj
								Ktoj = (HULLDRAG + bfkin.wr + bfkin.ar - bfkin.sfc) / (bfkin.rho * pow(bfkin.Vs, 2) * pow(2 * bfkin.R, 2) * pow(W, 2));

								//�Q��Ktoj�A����KTo.exe�A�p��ja
								sprintf(filename, "Ktoj.dat");
								fp = fopen(filename, "w");
								fprintf(fp, "%lf", Ktoj);
								fclose(fp);

								sprintf(cmd, "KTo.exe");
								outexe(cmd);
								Sleep(3000);
								fp = fopen("ja.dat", "r");
								fscanf(fp, "%lf", &ja);
								fclose(fp);

								//�p��J��
								prop.ADVCO = ja / W;
								debug("J Done", getpid());
							}
							debug("�Ĥ@���H�G���k�p��J�ȩΦۭq�qJ��  Done!", getpid());
						}
						else if ((ite - (bfkin.itestar - bfkin.ci)) / bfkin.ci == 2)
						{
							debug("Round2", getpid());
							//�����e�@��hXXX.geo��J��
							sprintf(filename, "h%d.geo", ite - bfkin.ci);
							fp = fopen(filename, "r");
							fgets(str, 80, fp);
							fgets(str, 80, fp);
							fgets(str, 80, fp);
							fscanf(fp, "%lf", &prop.ADVCO);
							fclose(fp);

							//�����e�@����hXXX.oup�̪�CD=0.0035��KT��
							sprintf(filename, "Kt-File.txt");
							fp = fopen(filename, "w");
							fprintf(fp, "h%d.oup", ite - bfkin.ci);
							fclose(fp);

							sprintf(cmd, "GetKt.exe");
							outexe(cmd);
							Sleep(1500);
							fp = fopen("KT.dat", "r");
							fscanf(fp, "%lf", &KT);
							fclose(fp);
							debug("*.KT output!", getpid());

							//�p��e�@��THRUST
							THRUST = KT * bfkin.rho * pow(bfkin.Vs, 2) * pow(2 * bfkin.R, 2) / pow(prop.ADVCO, 2);

							//��X���e�@����Error
							Error = (HULLDRAG + bfkin.wr + bfkin.ar - bfkin.sfc) - THRUST;
							sprintf(filename, "Error%d.dat", ite - bfkin.ci);
							fp3 = fopen(filename, "w");
							fprintf(fp3, "%lf\n", Error);
							fclose(fp3);
							debug("Error output Done!", getpid());


							//����������@��patpans11�ұo��hXXX.oup�̪�1-w��
							sprintf(filename, "Kt-File.txt");
							fp = fopen(filename, "w");
							fprintf(fp, "h%d.oup", ite);
							fclose(fp);

							sprintf(cmd, "GetKt.exe");
							outexe(cmd);
							Sleep(1500);
							fp = fopen("W.dat", "r");
							fscanf(fp, "%lf", &W);
							fclose(fp);

							//�p��Ktoj
							Ktoj = (HULLDRAG + bfkin.wr + bfkin.ar - bfkin.sfc) / (bfkin.rho * pow(bfkin.Vs, 2) * pow(2 * bfkin.R, 2) * pow(W, 2));

							//�Q��Ktoj�A����KTo.exe�A�p��ja
							sprintf(filename, "Ktoj.dat");
							fp = fopen(filename, "w");
							fprintf(fp, "%lf", Ktoj);
							fclose(fp);

							sprintf(cmd, "KTo.exe");
							outexe(cmd);
							Sleep(3000);
							fp = fopen("ja.dat", "r");
							fscanf(fp, "%lf", &ja);
							fclose(fp);

							//�p��J��
							prop.ADVCO = ja / W;

							debug("�ĤG���H�G���k�p��J��  Done!", getpid());
						}
						else if ((ite - (bfkin.itestar - bfkin.ci)) / bfkin.ci > 2)
						{
							debug("Round3", getpid());
							//�����e�@��hXXX.geo��J��
							sprintf(filename, "h%d.geo", ite - bfkin.ci);
							fp = fopen(filename, "r");
							fgets(str, 80, fp);
							fgets(str, 80, fp);
							fgets(str, 80, fp);
							fscanf(fp, "%lf", &prop.ADVCO);
							fclose(fp);

							//�����e�G��hXXX.geo��J��
							sprintf(filename, "h%d.geo", ite - 2 * bfkin.ci);
							fp = fopen(filename, "r");
							fgets(str, 80, fp);
							fgets(str, 80, fp);
							fgets(str, 80, fp);
							fscanf(fp, "%lf", &ADVCO0);
							fclose(fp);

							//�����e�@����hXXX.oup�̪�CD=0.0035��KT��
							sprintf(filename, "Kt-File.txt");
							fp = fopen(filename, "w");
							fprintf(fp, "h%d.oup", ite - bfkin.ci);
							fclose(fp);

							sprintf(cmd, "GetKt.exe");
							outexe(cmd);
							Sleep(1500);
							fp = fopen("KT.dat", "r");
							fscanf(fp, "%lf", &KT);
							fclose(fp);
							debug("*.KT output!", getpid());

							//�p��e�@��THRUST
							THRUST = KT * bfkin.rho * pow(bfkin.Vs, 2) * pow(2 * bfkin.R, 2) / pow(prop.ADVCO, 2);

							//��X���e�@����Error
							Error = (HULLDRAG + bfkin.wr + bfkin.ar - bfkin.sfc) - THRUST;
							sprintf(filename, "Error%d.dat", ite - bfkin.ci);
							fp3 = fopen(filename, "w");
							fprintf(fp3, "%lf\n", Error);
							fclose(fp3);
							debug("Err output Done!", getpid());
							Err = Error;


							//�����e�G����hXXX.oup�̪�CD=0.0035��KT��
							sprintf(filename, "Kt-File.txt");
							fp = fopen(filename, "w");
							fprintf(fp, "h%d.oup", ite - 2 * bfkin.ci);
							fclose(fp);

							sprintf(cmd, "GetKt.exe");
							outexe(cmd);
							Sleep(1500);
							fp = fopen("KT.dat", "r");
							fscanf(fp, "%lf", &KT0);
							fclose(fp);
							debug("*.KT0 output!", getpid());

							//�p��e�G��THRUST
							THRUST0 = KT0 * bfkin.rho * pow(bfkin.Vs, 2) * pow(2 * bfkin.R, 2) / pow(ADVCO0, 2);

							//�����e�G����Error
							sprintf(filename, "Error%d.dat", ite - 2 * bfkin.ci);
							fp = fopen(filename, "r");
							fscanf(fp, "%lf", &Err0);
							fclose(fp);
							debug("Err0 load Done!", getpid());

							//���y�k�p��J��
							if (THRUST == THRUST0) {
								prop.ADVCO = ADVCO0;
							}
							else
							{
								prop.ADVCO = ((HULLDRAG + bfkin.wr + bfkin.ar - bfkin.sfc) - THRUST) / (Err0 - Err) * (prop.ADVCO - ADVCO0) + prop.ADVCO;
							}

							debug("�ĤT���H���y�k�p��J�ȭp��J��  Done!", getpid());
						}

					}
					//�إ�hXXX.geo�AADVCO�ĥ�--�A�~���p���--�A�ð��J�y�P��t��~�A�������prop.geo
					sprintf(filename, "h%d.geo", ite);
					fp = fopen(filename, "w");
					fprintf(fp, label);
					fprintf(fp, "%d %d %d %d %d\n", prop.NX, prop.NBLADE, prop.NC, prop.MR, prop.NTMP);
					fprintf(fp, "%.4lf %d %d %.4lf %.4lf %.4lf\n", prop.RHUB, prop.NHBU, prop.MHBT, prop.XHBU, prop.XHBD, prop.XHBT);
					fprintf(fp, "%.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf\n", prop.ADVCO, prop.RULT, prop.RHULT, prop.DCD, prop.XULT, prop.DTPROP, prop.XUWDK);
					for (j = 0; j != 7; ++j) {
						for (i = 0; i != prop.NX; ++i)fprintf(fp, "%8.5lf ", prop.NPARAMETER[j][i]);
						fprintf(fp, "\n");
					}
					//�J�y�P��t�񬰶�L�J�y/Vs-���ɳt��/Vs
					for (i = 0; i != prop.NX; ++i)fprintf(fp, "%.5lf ", Vx_total[i] / bfkin.Vs - Ux_bem[i]);
					fprintf(fp, "\n");
					for (i = 0; i != prop.NX; ++i)fprintf(fp, "%.5lf ", Vr_total[i] / bfkin.Vs - Ur_bem[i]);
					fprintf(fp, "\n");
					for (i = 0; i != prop.NX; ++i)fprintf(fp, "%.5lf ", Vt_total[i] / bfkin.Vs - Ut_bem[i]);
					fprintf(fp, "\n");
					//hXXX.geo�̫�4��]��0
					for (j = 0; j != 4; ++j) {
						for (i = 0; i != prop.NX; ++i) { fprintf(fp, "0.00000 "); }
						fprintf(fp, "\n");
					}
					fclose(fp);

					//�إ�hXXX.adm�A�������prop.adm
					Output_ite_ADM();

					/*========================================================================================================
					*                                   ��X�P�V�����t�סB���ɳt�סB��۴��t��
					* ========================================================================================================
					*/
					sprintf(filename, "velocity_%d.txt", ite - bfkin.ci);
					fp = fopen(filename, "w");
					fprintf(fp, "r/R  Vx_total  Vr_total  Vt_total  Vx_induce  Vr_induce  Vt_induce  Vx_total-induce  Vr_total-induce  Vt_total-induce\n");
					for (i = 0; i != prop.NX; ++i)fprintf(fp, "%10.6lf %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf\n", prop.NPARAMETER[0][i], Vx_total[i] / bfkin.Vs, Vr_total[i] / bfkin.Vs, Vt_total[i] / bfkin.Vs, Ux_bem[i], Ur_bem[i], Ut_bem[i], Vx_total[i] / bfkin.Vs - Ux_bem[i], Vr_total[i] / bfkin.Vs - Ur_bem[i], Vt_total[i] / bfkin.Vs - Ut_bem[i]);
					fprintf(fp, "%10d", size);
					fclose(fp);

					debug("velocity_XX.txt output!", getpid());
					Sleep(1000);
					/*========================================================================================================
					*                                   ����ĤG��patpans11
					* ========================================================================================================
					*/
					sprintf(filename, "h%d", ite);
					fp = fopen("patpans11.kin", "w");
					fprintf(fp, filename);
					fclose(fp);

					sprintf(cmd, "cmd.exe /c patpans11.exe<patpans11.kin");
					outexe(cmd);
					Sleep(2000);
					debug("Patpans2 Done!", getpid());

					//===================================����hXXX.oup�̪�CD=0.0035��KT�ȡA�Ω�p��THRUST========================		
					sprintf(filename, "Kt-File.txt");
					fp = fopen(filename, "w");
					fprintf(fp, "h%d.oup", ite);
					fclose(fp);

					sprintf(cmd, "GetKt.exe");
					outexe(cmd);
					Sleep(1500);
					fp = fopen("KT.dat", "r");
					fscanf(fp, "%lf", &KT);
					fclose(fp);

					//====================================================�p��THRUST==========================================
					THRUST = KT * bfkin.rho * pow(bfkin.Vs, 2) * pow(2 * bfkin.R, 2) / pow(prop.ADVCO, 2);

					//================================================��XhulldragXXX.dat=====================================
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

					//=================================================��shulldrag_ite.dat====================================
					fp = fopen("hulldrag_ite.txt", "w");
					fprintf(fp, "%d\n", bfkin.itestar);
					fprintf(fp, "%d\n", ite);
					fclose(fp);

					//Ū�����O
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
					ite = (int)Iteration[0][0];
					//=====================================================��X���O=============================================
					sprintf(filename, "force_%d.txt", ite);
					fp = fopen(filename, "w");
					for (i = 0; i != prop.MR; ++i) {
						fprintf(fp, "%lf %lf %lf %lf\n", r_F[i] * bfkin.R, bfx[i], bfr[i], bft[i]);
					}
					fclose(fp);
					/*========================================================================================================
					*											����chebyToCn.exe���oCn
					* ========================================================================================================
					*/
					//bfx - ����chebyToCn.exe���oCn
					sprintf(filename, "fx%d.dat", ite);
					fp = fopen(filename, "w");
					for (i = 0; i != prop.MR; ++i) { fprintf(fp, "%lf %lf\n", r_F[i] * bfkin.R, bfx[i]); }
					fclose(fp);
					sprintf(filename, "fx.dat");
					fp = fopen(filename, "w");
					for (i = 0; i != prop.MR; ++i) { fprintf(fp, "%lf %lf\n", r_F[i] * bfkin.R, bfx[i]); }
					fclose(fp);

					fp = fopen("chebyToCn.kin", "w");
					sprintf(str, "fx\n6");
					fprintf(fp, str);
					fclose(fp);

					sprintf(cmd, "chebyToCn.exe");//�ť߸��ഫ
					outexe(cmd);

					Sleep(1000);

					//bfr - ����chebyToCn.exe���oCn
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
					Sleep(1000);
					//bft - ����chebyToCn.exe���oCn
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
					Sleep(1000);
					debug("ChebyToCn Done!", getpid());











				}

				Sleep(1000);
				//�ϥ�free()�k�ٵ��O����				
				for (i = 0; i != 14; i++) 
				{
					free(prop.NPARAMETER[i]);
				}
				free(prop.NPARAMETER);
				free(Vx_ave); free(Vt_ave); free(Vr_ave);
				free(Vx_sum); free(Vt_sum); free(Vr_sum);
				free(r_ave); free(r_sum); free(n);
				free(fx); free(ft); free(fr);
				free(bfx); free(bft); free(bfr);
				free(r_F); free(r_x); free(r_U);
				free(Vx_total); free(Vt_total); free(Vr_total);
				free(Ux); free(Ut); free(Ur);
				free(Ux_bem); free(Ut_bem); free(Ur_bem);

				if (meshload == 0) 
				{
					if (threadboy == 0) 
					{
						for (i = 0; i != (bfkin.plane_x_numb * 4 * bfkin.plane_c_numb * bfkin.plane_r_numb); ++i)
						{
							free(mesh[i]);
						}
						free(mesh);
					}
				}
				debug("Free Done!", getpid());

				Sleep(1000);



			}


		}


		free(bfkin.bfkin_hub);
		if (ite == 0) //�M����0�B��Xmesh.txt��process���x�}
		{
			if (threadboy == 1) 
			{
				free(newsize); free(threadnumber);
				for (i = 0; i != (bfkin.plane_x_numb * 4 * bfkin.plane_c_numb * bfkin.plane_r_numb); ++i) 
				{
					free(mesh[i]);
				}
				free(mesh);

			}
		}


	


		if (ite >= bfkin.itestar && ite % bfkin.ci == reitestar  && meshload == 0)
		{
			if (threadboy == 1) 
			{
				free(newsize); free(threadnumber);
				for (i = 0; i != (bfkin.plane_x_numb * 4 * bfkin.plane_c_numb * bfkin.plane_r_numb); ++i) 
				{
					free(mesh[i]);
				}
				free(mesh);
			}
		}
		


		debug("Free Done!", getpid());






	}
}
