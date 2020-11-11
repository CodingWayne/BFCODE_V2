#include <stdio.h>
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
