
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



//[result]:the array of values returned by the user function.(Force/Volume);  [size]:the number of elements in the result array.(int);  [centroid]:���椤���I���y��;  [Velocity]:���椤���I���t��; [Iteration]:the iteration number in CoordReal precision.(int); [Angle]:����pitch���|��(�F������); [Cell_id]:��L���椧�s��; [Drag]:��ߪ��O;
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
		//-----------------------------------------------------------------------����Ѽƪ�l��----------------------
		meshload = 0; minCellid = 999999999, np = 0, threadboy = 0, * newsize, * threadnumber, all_size = 0, non = 0;
		//----------------------------------------------------------------------------------------------------------

		//-----------------------------------------------------------------------������L����------------------------
		if (ite == 0) 
		{
			debug("------------ite=0--------", getpid());
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
			debug("STEP[1]-2   output minCellid file", getpid());			
			/*
			*STEP[2]
			*��X�U�֤ߤ�minCellid
			*/			
			sprintf(str, ".\\mesh\\%d_%d.txt", getpid(), ite);
			fp1 = fopen(str, "w");
			fprintf(fp1, "%d", minCellid);
			fclose(fp1);	
			debug("STEP[2]     output minCellid file", getpid());
			/*
			 *STEP[3]
			 *�p��thread��Ƨ��U���h���ɮ�=���h��process�b�B�z����
			 */
			np = 0;
			np = count_file_num_in_a_folder(".\\thread");
			debug("STEP[3]     np calculation Done!", getpid());
			/*
			 *STEP[4]
			 *
			 *�H0�}�l�V�W�W�[�ӧ�̤p�֤ߧǡA�åH�̤p�֤ߧǮ֤߶}���ɮ�(allthread.txt)�����U�֤ߤ��֤ߧǩMsize�A�קK��L�֤߶}�Ҧ��ɡA��L�֤߸��L���B�J
			 *threadboy=1��ܤw�b�̤pprocess�̫ظm�x�}���\�åH�̤p�֤߶}�ҦU�ɮ׬����b�x�}��
			 */
			j = 0;
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
			/*			 
			 *��threadboy=1��,�}�l�ظm����x�}		
			 *�H�Uif�P�_���̬ҥѳ̤pprocessID����
			 */
			if (threadboy == 1) 
			{
				/*
				 *STEP[5]
				 * �ظm2��mesh�x�}  mesh[��L�����][8]
				 */
				mesh = (double**)malloc((bfkin.plane_x_numb * 4 * bfkin.plane_c_numb * bfkin.plane_r_numb) * sizeof(double*));
				for (i = 0; i != (bfkin.plane_x_numb * 4 * bfkin.plane_c_numb * bfkin.plane_r_numb); ++i) 
				{
					mesh[i] = (double*)malloc(8 * sizeof(double));
				}				
				debug("STEP[5]     MESH[][] malloc done!", getpid());

				/*
				* STEP[6]
				* �M�w�̤p����id�ÿ�X,��X��R��
				*/				
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

				/*
				* STEP[7]
				* ��Jmesh[][]
				* �̥H�e������np,�̧Ƕ}�ҦUthread��X��(thread_pid_ite.txt),���]����ƾ�22500��,�����s������0~22499�H��K���ާ@�x�},�]��kk�ݬ�22500-minCellid(1)
				* �N(thread_pid_ite.txt)�̦UCellid�����������Ʃ��mesh[id-1][x,y,z,Vx,Vy,Vz,?,?]��,���ɦU(thread_pid_ite.txt)�̭��Ƥ�Cellid��Ʒ|�л\,�i�o����22500�Ӻ����Ƥ��x�}
				* �R��(thread_pid_ite.txt)
				*/	
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

				/*
				*STEP[8] 
				* �إ�mesh�x�}�����L�|�V���Z���ή|�� mesh[22500][6]&mesh[22500][7]
				*/
				for (i = 0; i != (bfkin.plane_x_numb * 4 * bfkin.plane_c_numb * bfkin.plane_r_numb); ++i) 
				{
					/*
					* p_xo,p_yo,p_zo�����L���椤���I,���ɲ��^���I
					*/
					mesh[i][0] = mesh[i][0] - bfkin.p_xo;
					mesh[i][1] = mesh[i][1] - bfkin.p_yo;
					mesh[i][2] = mesh[i][2] - bfkin.p_zo;

					//�ھڱ��O��V �ק�����V
					if (bfkin.bow_dir == 1)
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
					else {
						mesh[i][7] = (double)acos(mesh[i][2] / mesh[i][6]); //�|��
					}
				}
				sprintf(str, "del .\\thread\\allthread%d.txt", threadnumber[0]);
				system(str);
				debug("STEP[8]     mesh load done!", getpid());
				/*
				*STEP[9] 
				*��Xmesh[ID-1][C123_V123__r_theta]�奻
				*/
				sprintf(str, ".\\mesh\\mesh.txt");
				fp = fopen(str, "w");
				for (i = 0; i < (bfkin.plane_x_numb * 4 * bfkin.plane_c_numb * bfkin.plane_r_numb); i++) 
				{
					fprintf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf\n", mesh[i][0], mesh[i][1], mesh[i][2], mesh[i][3], mesh[i][4], mesh[i][5], mesh[i][6], mesh[i][7]);
				}
				fclose(fp);
				debug("STEP[9]     mesh[ID-1][Cx_Cy_Cz_Vx_Vy_Vz_r_theta] output done!", getpid());
			}
			
		}
		//----------------------------------------------------------------------------------------------------------
		
		
	}
		
	



	
}
