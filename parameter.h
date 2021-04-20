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
void Output_ite_ADM();
char geoname[50];
struct PROP
{
	//NX:number of radii for the input geometry parameters for patpans11; NBLADE:number fo blades;  MR:number of spanwise panels in patpans11
	int NX, NBLADE, NC, MR, NTMP, NHBU, MHBT;
	double RHUB, XHBU, XHBD, XHBT, ADVCO, RULT, RHULT, DCD, XULT, DTPROP, XUWDK;//ADVCO:advance coefficient(Js=Vs/n/D)�A�ΥH�x�s.geo��J�ȡA�Ϋe�@��hXXX.geo��J�ȡA�Ω���y�k��J��
	double** NPARAMETER;//NPARAMETER[][]:�ΥH�x�sprop.geo�ĥ|�椧�᪺�ƭ�;
};
void Load_PROP_GEO(struct PROP* prop);
double J_Correcction( struct BFKIN* bfkin, struct PROP* prop);
double J_Correcction_coe[5];
void adm_rev( struct BFKIN* bfkin, struct PROP* prop);
//--------------------------------------------------------------------------------
double  pitch;
FILE* fp5, * fp3, * fp, * fp1, * fp2;
int ite, reitestar;
char str[50], label[100], cmd[50], filename[50];
int i, j;
double PI;
double r, theta;//r:�ΥH�x�s����Ҧb���b�|�j�p theta:�ΥH�x�s����Ҧb���|�פj�p R:��L�b�| x:��L�p�� xp;�bx�b�W����m
double HULLDRAG, THRUST, THRUST0, KT, KT0, ADVCO0, Error, W, Err, Err0;
/*HULLDRAG:�Ω��x�shulldrag.csv�̪��̫�@��(�]�N�O�C�@�����e�@�B��drag)
* THRUST:�ΥH�x�s���ı��O,�Ϋe�@�������ı��O
* THRUST0:�ΥH�x�s�e�G�������ı��O
* KT:�ΥH�x�shXXX.oup�̪�CD=0.0035��KT�ȡA�Ϋe�@����hXXX.oup�̪�CD=0.0035��KT�ȡA�Ω�p��THRUST;
* KT0:�ΥH�x�s�e�G����hXXX.oup�̪�CD=0.0035��KT�ȡA�Ω�p��THRUST0;
* ADVCO0:�x�s�e�G��hXXX.geo��J�ȡA�Ω���y�k��J��;
* Error:���Ī��O�P���ı��O���t��;
* W:�ΥH�x�sW.dat�̪��ƭ�(hXXX.oup�̪�1-w);
* Err:�ΥH�x�s�e�@����Error;
* Err0:�ΥH�x�s�e�G����Error;
*/
double Ktoj, ja;
//Cx, Cr, Ct:�ΥH�x�s��n�O�b�|�V���G�Ƭ���6�ӫY��(C0~C5); Tn:�ΥH�x�s��n�O�b�|�V���G�Ƭ���6�ӫY�ƪ����ܼƨ��(T0~T6); sum:�Ω��x�s��n�O�����̤�CnTn�`�M;origin �x�s��n�O
double  Cx[6], Cr[6], Ct[6], Tn[6], sum, originfx, originfr, originft;
double tmp;//tmp�x�sadm���t�ץ���
/*-------------------------�����ܼ�------------------------------------------------------
 *meshload���P�_�ܼơA�짹�����=1
 *minCellid���̤p����s��
 *threadboy���P�_�ܼơA����1�ɪ��ɥγ�@process�}�ҹ�����mesh�ɡA���ɫإ�newsize�x�}��threadnumber�x�}�A���O�}�ҦUthread%d_%d.txt�A��Jpid��threadnumber,size��newsize��
 *allsize=newsize�x�}�U�����[�`�A���O�������L�����(�|�����ƶ�)
 *kk���e=�Ȧs�ܼ�
 *numb,non�Ȭ��O����
 * tmp_size���Ȧssize���ܼ�
 */
int  meshload, minCellid, np, threadboy, * newsize, * threadnumber, all_size;
double** mesh;
int numb, kk;
double non;


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
void Output_bfkin(struct BFKIN* bfkin, int* reitestar, double* pitch);
//-------------------------------------------------------------------------------------
double* Vx_ave, * Vt_ave, * Vr_ave, * Vx_sum, * Vt_sum, * Vr_sum, * r_ave, * r_sum, * n;
double* fx, * fr, * ft, * bfx, * bfr, * bft, * r_F, * r_x, * r_U, * Vx_total, * Vt_total, * Vr_total, * Ux, * Ut, * Ur, * Ux_bem, * Ut_bem, * Ur_bem;
//r_x[]:�ΥH�x�s�������V���q(MR�Ӭq)�̡A�]�t�Y�����C�I��r/R�A�ƶq��MR+1�ӡA(���V���q��sine spanwise) 
//r_F[]:�ΥH�x�spatpans11.exe�p�⤧�O���|�V��m�A�ƶq��MR�ӡA(���V���q��sine spanwise) 
//fx[],fr[],ft[]:patpans11.exe�p�⤧�O?
//bfx[],bfr[],bft[] = fx~r~t[i]/0.5*rho*Vs*Vs*(r_x[i+1]-r_x[i])*R;
//Vx_total[],Vt_total[],Vr_total[]:�ΥH�x�s�U�|�V��m���椧�P�V�����t��  SPLINE�� NX ���I�᪺���G�A��@��L���߭����J�y�t��
//r_U[]:�ΥH�x�ssfpv11.exe�p�⤧���ɳt�ת��|�V��m�A�ƶq��MR�ӡA(���V���q��sine spanwise) 
//Ux[],Ut[],Ur[]:�ΥH�x�ssfpv11.exe�p�⤧�U�|�V��m�����ɳt�׻P��t�����(Vinduced/Vs)�A�ƶq��MR��
//Ux_bem[],Ut_bem[],Ur_bem[]:�ΥH�x�ssfpv11.exe�p�⤧�U�|�V��m�����ɳt�׻P��t�����(Vinduced/Vs) SPLINE�� NX ���I�᪺���G�A��@��L���߭������ɳt��

#endif