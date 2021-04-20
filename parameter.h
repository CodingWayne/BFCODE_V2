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
void Output_ite_ADM();
char geoname[50];
struct PROP
{
	//NX:number of radii for the input geometry parameters for patpans11; NBLADE:number fo blades;  MR:number of spanwise panels in patpans11
	int NX, NBLADE, NC, MR, NTMP, NHBU, MHBT;
	double RHUB, XHBU, XHBD, XHBT, ADVCO, RULT, RHULT, DCD, XULT, DTPROP, XUWDK;//ADVCO:advance coefficient(Js=Vs/n/D)，用以儲存.geo之J值，或前一輪hXXX.geo之J值，用於牛頓法找J值
	double** NPARAMETER;//NPARAMETER[][]:用以儲存prop.geo第四行之後的數值;
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
double r, theta;//r:用以儲存網格所在的半徑大小 theta:用以儲存網格所在的徑度大小 R:圓盤半徑 x:圓盤厚度 xp;在x軸上的位置
double HULLDRAG, THRUST, THRUST0, KT, KT0, ADVCO0, Error, W, Err, Err0;
/*HULLDRAG:用於儲存hulldrag.csv裡的最後一項(也就是每一輪之前一步的drag)
* THRUST:用以儲存有效推力,或前一輪之有效推力
* THRUST0:用以儲存前二輪之有效推力
* KT:用以儲存hXXX.oup裡的CD=0.0035之KT值，或前一輪之hXXX.oup裡的CD=0.0035之KT值，用於計算THRUST;
* KT0:用以儲存前二輪之hXXX.oup裡的CD=0.0035之KT值，用於計算THRUST0;
* ADVCO0:儲存前二輪hXXX.geo之J值，用於牛頓法找J值;
* Error:有效阻力與有效推力的差值;
* W:用以儲存W.dat裡的數值(hXXX.oup裡的1-w);
* Err:用以儲存前一輪之Error;
* Err0:用以儲存前二輪之Error;
*/
double Ktoj, ja;
//Cx, Cr, Ct:用以儲存體積力在徑向分佈化為的6個係數(C0~C5); Tn:用以儲存體積力在徑向分佈化為的6個係數的單變數函數(T0~T6); sum:用於儲存體積力公式裡之CnTn總和;origin 儲存體積力
double  Cx[6], Cr[6], Ct[6], Tn[6], sum, originfx, originfr, originft;
double tmp;//tmp儲存adm降速修正值
/*-------------------------網格變數------------------------------------------------------
 *meshload為判斷變數，抓完網格後=1
 *minCellid為最小網格編號
 *threadboy為判斷變數，等於1時表此時用單一process開啟對應之mesh檔，此時建立newsize矩陣及threadnumber矩陣，分別開啟各thread%d_%d.txt，輸入pid到threadnumber,size到newsize裡
 *allsize=newsize矩陣各元素加總，但是不等於圓盤網格數(會有重複項)
 *kk為占=暫存變數
 *numb,non僅為記錄用
 * tmp_size為暫存size的變數
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
//r_x[]:用以儲存葉片垂向分段(MR個段)裡，包含頭尾之每點的r/R，數量為MR+1個，(垂向分段為sine spanwise) 
//r_F[]:用以儲存patpans11.exe計算之力的徑向位置，數量為MR個，(垂向分段為sine spanwise) 
//fx[],fr[],ft[]:patpans11.exe計算之力?
//bfx[],bfr[],bft[] = fx~r~t[i]/0.5*rho*Vs*Vs*(r_x[i+1]-r_x[i])*R;
//Vx_total[],Vt_total[],Vr_total[]:用以儲存各徑向位置網格之周向平均速度  SPLINE成 NX 個點後的結果，當作圓盤中心面之入流速度
//r_U[]:用以儲存sfpv11.exe計算之誘導速度的徑向位置，數量為MR個，(垂向分段為sine spanwise) 
//Ux[],Ut[],Ur[]:用以儲存sfpv11.exe計算之各徑向位置的誘導速度與船速的比值(Vinduced/Vs)，數量為MR個
//Ux_bem[],Ut_bem[],Ur_bem[]:用以儲存sfpv11.exe計算之各徑向位置的誘導速度與船速的比值(Vinduced/Vs) SPLINE成 NX 個點後的結果，當作圓盤中心面之誘導速度

#endif