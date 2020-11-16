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

//---------------------prop.geo變數------------------------------------------------
char adm1[100], adm2[100], adm3[100], adm4[100], adm5[100], adm6[100];
void Load_PROP_ADM();
void CHECK_PROP_ADM();
char geoname[50];
struct PROP
{
	//NX:number of radii for the input geometry parameters for patpans11; NBLADE:number fo blades;  MR:number of spanwise panels in patpans11
	int NX,NBLADE,NC,MR,NTMP,NHBU,MHBT;
	double RHUB,XHBU,XHBD,XHBT,ADVCO,RULT,RHULT,DCD,XULT,DTPROP,XUWDK;//ADVCO:advance coefficient(Js=Vs/n/D)，用以儲存.geo之J值，或前一輪hXXX.geo之J值，用於牛頓法找J值
	double** NPARAMETER;//NPARAMETER[][]:用以儲存prop.geo第四行之後的數值;
};
void Load_PROP_GEO(struct PROP* prop);
void Output_PROP_GEO(struct PROP* prop);
//--------------------------------------------------------------------------------
double reitestar, pitch;
FILE* fp5, * fp3, * fp, * fp1;
int ite;
char str[50],label[100];
int i, j;

/*-------------------------網格變數------------------------------------------------------
 *meshload為判斷變數，抓完網格後=1
 *minCellid為最小網格編號
 *threadboy為判斷變數，等於1時表此時用單一process開啟對應之mesh檔，此時建立newsize矩陣及threadnumber矩陣，分別開啟各thread%d_%d.txt，輸入pid到threadnumber,size到newsize裡
 *allsize=newsize矩陣各元素加總，但是不等於圓盤網格數(會有重複項)
 *kk為占=暫存變數
 *numb,non僅為記錄用
 * tmp_size為暫存size的變數
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