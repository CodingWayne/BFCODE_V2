#ifndef test_H
#define test_H

#include <stdio.h>
#include <stdlib.h>
#include "uclib.h"

struct BFKIN
{
	int np;
	int itestar;
	int ci;
	int wake_offon;
	int wake_ci;
	int wake_onite;

	int self_offon;
	int firstJ_offon;
	double firstJ;
	int bodyforce_offon;
	int P_3;

	int Unsteady_offon;
	int unsteady_c;
	int UN_ci;
	double theta_one_blade;
	int tnb;
	int disk_offen;
	int UNB;

	double R;
	double x;
	double xp;
	int plane_r_numb;
	int plane_x_numb;
	int plane_c_numb;

	double Vs;
	double rho;
	double wr;
	double ar;
	double sfc;

	int bow_dir;
	double p_xo;
	double p_yo;
	double p_zo;
	int st;

	double* bfkin_hub;

};

//---------------------prop.geo�ܼ�------------------------------------------------
char adm1[100], adm2[100], adm3[100], adm4[100], adm5[100], adm6[100];
void Load_PROP_ADM();
void CHECK_PROP_ADM();
char geoname[50];
struct PROP
{
	//NX:number of radii for the input geometry parameters for patpans11; NBLADE:number fo blades;  MR:number of spanwise panels in patpans11
	int NX,NBLADE,NC,MR,NTMP,NHBU,MHBT;
	double RHUB,XHBU,XHBD,XHBT,ADVCO,RULT,RHULT,DCD,XULT,DTPROP,XUWDK;//ADVCO:advance coefficient(Js=Vs/n/D)�A�ΥH�x�s.geo��J�ȡA�Ϋe�@��hXXX.geo��J�ȡA�Ω���y�k��J��
	double** NPARAMETER;//NPARAMETER[][]:�ΥH�x�sprop.geo�ĥ|�椧�᪺�ƭ�;
};
void Load_PROP_GEO(struct PROP* prop);
void Output_PROP_GEO(struct PROP* prop);
//--------------------------------------------------------------------------------
double reitestar, pitch;
FILE* fp5, * fp3, * fp, * fp1;
int ite;
char str[50],label[100];
int i, j;

/*-------------------------�����ܼ�------------------------------------------------------
 *meshload���P�_�ܼơA�짹�����=1
 *minCellid���̤p����s��
 *threadboy���P�_�ܼơA����1�ɪ��ɥγ�@process�}�ҹ�����mesh�ɡA���ɫإ�newsize�x�}��threadnumber�x�}�A���O�}�ҦUthread%d_%d.txt�A��Jpid��threadnumber,size��newsize��
 *allsize=newsize�x�}�U�����[�`�A���O�������L�����(�|�����ƶ�)
 *kk���e=�Ȧs�ܼ�
 *numb,non�Ȭ��O����
 * tmp_size���Ȧssize���ܼ�
 */
int  meshload,minCellid,np, threadboy, * newsize, * threadnumber, all_size ,tmp_size;
double** mesh;
int numb, kk;
double non;
double **tmp_Cell_id,**tmp_centroid,**tmp_Velocity;
void MESH_STEP_1();
void MESH_STEP_2();
void MESH_STEP_3();
void MESH_STEP_4();
void MESH_STEP_5(struct BFKIN* bfkin);
void MESH_STEP_6();
void MESH_STEP_7();
void MESH_STEP_8(struct BFKIN* bfkin);
void MESH_STEP_9(struct BFKIN* bfkin);
//-------------------------------------------------------------------------------------
void Load_bfkin(struct BFKIN* bfkin);
void Output_bfkin(struct BFKIN* bfkin, double* reitestar, double* pitch);


#endif