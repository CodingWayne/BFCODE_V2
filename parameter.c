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
	fscanf(fp, "%d%d%d%d%d%d", &bfkin->np, &bfkin->itestar, &bfkin->ci, &bfkin->wake_offon, &bfkin->wake_ci, &bfkin->wake_onite);  //�qbf->kinŪ�J; np:�ҨϥΪ��֤߼�; itestar:�}�l�p��ci�e������p��B��; ci:�C�X�B���N�@��; wake_offon:�O�_�}��wake��X; wake_ci:�C�X�B��Xwake�@��; wake_onite:�}�l��Xwake���B��; 
	fscanf(fp, "%d%d%lf%d%d", &bfkin->self_offon, &bfkin->firstJ_offon, &bfkin->firstJ, &bfkin->bodyforce_offon, &bfkin->P_3);		//self_offon:�O�_�}��self;	  firstJ_offon:�O�_�}�Ҧ۩w�q�Ĥ@��J��;  firstJ:�۩w�q���Ĥ@��J��;   bodyforce_offon:�O�_�}��bodyforce;   
	fscanf(fp, "%d%d%d%lf%d%d%d", &bfkin->Unsteady_offon, &bfkin->unsteady_c, &bfkin->UN_ci, &bfkin->theta_one_blade, &bfkin->tnb, &bfkin->disk_offen, &bfkin->UNB);	     // Blade_offon:�O�_�}�Ҹ����Ҧ� ;  Unsteady_offon:�O�_�}��Unsteady bodyforce ; unsteady_c: ubf�}�l�i��B�� ; theta_one_blade: ��v���խ��n  ;  tnb: �C�@�����N���ո�����ʰ��
	fscanf(fp, "%lf%lf%lf%d%d%d", &bfkin->R, &bfkin->x, &bfkin->xp, &bfkin->plane_r_numb, &bfkin->plane_x_numb, &bfkin->plane_c_numb); //R:��L�b�|; x:��L�p��; xp:�bx�b�W����m; plane_r_numb:��Lr��V��������μ�; plane_x_numb:��Lx��V��������μ�; plane_c_numb:��L�P��V��������μ�
	fscanf(fp, "%lf%lf%lf%lf%lf", &bfkin->Vs, &bfkin->rho, &bfkin->wr, &bfkin->ar, &bfkin->sfc);					//Vs:��t; rho:�K��; WR:���i���O; AR:�i�����[���O; SFC:�������O�ץ��q;
	fscanf(fp, "%d%lf%lf%lf%d", &bfkin->bow_dir, &bfkin->p_xo, &bfkin->p_yo, &bfkin->p_zo, &bfkin->st);	//��L���߮y�� ; ���F��V
	bfkin->bfkin_hub = (double*)malloc(bfkin->plane_x_numb * sizeof(double));
	for (i = 0; i != bfkin->plane_x_numb; ++i)fscanf(fp, "%lf", &bfkin->bfkin_hub[i]);  //x��V������hub��m
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
	prop->NPARAMETER = (double**)malloc(14 * sizeof(double*)); //�ΥH�x�sprop.geo�ĥ|�椧�᪺�ƭ� 
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
	 *�ϥ�getpid��X�Uprocess�W�۳B�z��������(thread_pid_ite.txt)
	 *���]��5��process�b���u�B�z�Athread��Ƨ��K�|��X5��(thread_pid_ite.txt)����5�Ӥ��P��pid
	 *(thread_pid_ite.txt)=[����s��,���椤���IX,���椤���Iy,���椤���Iz,���椤���I�t��X,���椤���I�t��y,���椤���I�t��z]
	 *����minCellid�b�U�O��process�̰O���̤p������s�� (minCellid���O�x�}�u�O���ܼơA���b�U��process�̥N���process��minCellid)
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
	*��X�U�֤ߤ�minCellid
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
	 *�p��thread��Ƨ��U���h���ɮ�=���h��process�b�B�z����
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
	 *�H0�}�l�V�W�W�[�ӧ�̤p�֤ߧǡA�åH�̤p�֤ߧǮ֤߶}���ɮ�(allthread.txt)�����U�֤ߤ��֤ߧǩMsize�A�קK��L�֤߶}�Ҧ��ɡA��L�֤߸��L���B�J
	 *threadboy=1��ܤw�b�̤pprocess�̫ظm�x�}���\�åH�̤p�֤߶}�ҦU�ɮ׬����b�x�}��
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
					fp1 = fopen(str, "w"); //��X�U�֤ߤ��֤ߧǡBsize
					newsize = (int*)malloc(np * sizeof(int));
					threadnumber = (int*)malloc(np * sizeof(int));
				}
			}
			fscanf(fp, "%d", &newsize[j]);
			threadnumber[j] = i; //newsize[j]�x�s�U�֤ߤ�size threadnumber[j]=i �x�s�U�֤ߤ��֤ߧ�
			fprintf(fp1, "%d %d\n", threadnumber[j], newsize[j]);
			all_size = all_size + newsize[j]; //�N�U�֤ߤ�size�`�M�A����|�R�����ƶ��A�Ϥ������L�����
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
	 * �ظm2��mesh�x�}  mesh[��L�����][8]
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
	* �M�w�̤p����id�ÿ�X,��X��R��
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
	sprintf(str, ".\\mesh\\%d_%d.txt", getpid(), ite); //��XminCellid file
	fp1 = fopen(".\\mesh\\minCellid.txt", "w");
	fprintf(fp1, "%d", minCellid);
	fclose(fp1);
	debug("STEP[6]     CPU data load!", getpid());
}
void MESH_STEP_7()
{
	/*
	* STEP[7]
	* ��Jmesh[][]
	* �̥H�e������np,�̧Ƕ}�ҦUthread��X��(thread_pid_ite.txt),���]����ƾ�22500��,�����s������0~22499�H��K���ާ@�x�},�]��kk�ݬ�22500-minCellid(1)
	* �N(thread_pid_ite.txt)�̦UCellid�����������Ʃ��mesh[id-1][x,y,z,Vx,Vy,Vz,?,?]��,���ɦU(thread_pid_ite.txt)�̭��Ƥ�Cellid��Ʒ|�л\,�i�o����22500�Ӻ����Ƥ��x�}
	* �R��(thread_pid_ite.txt)
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
	* �إ�mesh�x�}�����L�|�V���Z���ή|�� mesh[22500][6]&mesh[22500][7]
	*/
	for (i = 0; i != (bfkin->plane_x_numb * 4 * bfkin->plane_c_numb * bfkin->plane_r_numb); ++i)
	{
		/*
		* p_xo,p_yo,p_zo�����L���椤���I,���ɲ��^���I
		*/
		mesh[i][0] = mesh[i][0] - bfkin->p_xo;
		mesh[i][1] = mesh[i][1] - bfkin->p_yo;
		mesh[i][2] = mesh[i][2] - bfkin->p_zo;

		//�ھڱ��O��V �ק�����V
		if (bfkin->bow_dir == 1)
		{
			mesh[i][0] = mesh[i][0] * (-1);
			mesh[i][1] = mesh[i][1] * (-1);
			mesh[i][3] = mesh[i][3] * (-1);
			mesh[i][4] = mesh[i][4] * (-1);
		}
		//mesh[i][6]�����L�|�V���Z��r
		//mesh[i][7]�����L�|�V���|��theta
		mesh[i][6] = (double)pow(pow(mesh[i][1], 2.0) + pow(mesh[i][2], 2.0), 0.5);
		if (mesh[i][1] >= 0)
		{
			mesh[i][7] = (double)-acos(mesh[i][2] / mesh[i][6]); //�|��
		}
		else 
		{
			mesh[i][7] = (double)acos(mesh[i][2] / mesh[i][6]); //�|��
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
	*��Xmesh[ID-1][C123_V123__r_theta]�奻
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
