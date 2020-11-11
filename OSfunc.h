#pragma once
#include "uclib.h"
#include <math.h>
#include <stdio.h> 
#include <stdlib.h>
#include <string.h>
#include <windows.h>
#include <time.h>
#include <omp.h>
#include <process.h>

void outexe(char cmd[]) {
	STARTUPINFO si;
	PROCESS_INFORMATION pi;
	ZeroMemory(&si, sizeof(si));
	si.cb = sizeof(si);
	si.dwFlags = STARTF_USESHOWWINDOW;
	si.wShowWindow = SW_SHOWNORMAL;
	ZeroMemory(&pi, sizeof(pi));
	// Start the child process. 
	CreateProcess(NULL, cmd, NULL, NULL, FALSE, 0, NULL, NULL, &si, &pi);
	// Wait until child process exits.
	WaitForSingleObject(pi.hProcess, INFINITE);
	// Close process and thread handles.  
	CloseHandle(pi.hProcess);
	CloseHandle(pi.hThread);
}

void debug(char txt[], int pid) {
	FILE* fp;
	char str[100];
	sprintf(str, ".\\debug\\debug%d.txt", pid);	//�̦U�֤ߧǶ}�Ҥ��P��debugXXX.txt�A�קK���P�֤߶}�ҦP�@��debug.txt�A�y���֤߼g�ɮɫo����L�֤������ɮת����~
	fp = fopen(str, "a");
	fprintf(fp, "pid = %d, %s\n", pid, txt);
	fclose(fp);
}

int count_file_num_in_a_folder(char target_folder[]) {
	int count = -1;              //�ɮת�counter
	char szDir[128];			 //�nŪ������Ƨ�����}�C 
	WIN32_FIND_DATA FileData;    //���ۥثeŪ���쪺File�����СC
	HANDLE hList;                //���ۭnŪ������Ƨ������СC
	sprintf(szDir, "%s\\*", target_folder);
	if ((hList = FindFirstFile(szDir, &FileData)) == INVALID_HANDLE_VALUE)
		debug("No directories be found.", getpid());
	else {
		while (1) {
			if (!FindNextFile(hList, &FileData)) {
				if (GetLastError() == ERROR_NO_MORE_FILES)
					break;
			}
			count++;
		}
	}
	FindClose(hList);
	return count;
}