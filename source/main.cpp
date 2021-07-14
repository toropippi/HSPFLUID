#define WIN32_LEAN_AND_MEAN		// Exclude rarely-used stuff from Windows headers
#include <windows.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include "hsp3plugin.h"
#include "hsp3struct.h"


void FluidProcess(void);
void GetSpeedPoints(void);

int poissonloop = 8;

int bufsize = 128 * 128;
double* vx_after = new double[128 * 128];
double* vy_after = new double[128 * 128];
double* s = new double[128 * 128];


static double ref_doubleval;						// �Ԓl�̂��߂̕ϐ�
static float ref_floatval;						// �Ԓl�̂��߂̕ϐ�
static int ref_int32val;						// �Ԓl�̂��߂̕ϐ�

static void *reffunc( int *type_res, int cmd )
{
	//		�֐��E�V�X�e���ϐ��̎��s���� (�l�̎Q�Ǝ��ɌĂ΂�܂�)
	//
	//			'('�Ŏn�܂邩�𒲂ׂ�
	//
	if ( *type != TYPE_MARK ) puterror( HSPERR_INVALID_FUNCPARAM );
	if ( *val != '(' ) puterror( HSPERR_INVALID_FUNCPARAM );
	code_next();

	bool fDouble = false;
	bool fFloat = false;
	bool fInt = false;
	
	switch( cmd ) 
	{
	default:
		puterror( HSPERR_UNSUPPORTED_FUNCTION );
	}


	if ( *type != TYPE_MARK ) puterror( HSPERR_INVALID_FUNCPARAM );
	if ( *val != ')' ) puterror( HSPERR_INVALID_FUNCPARAM );
	code_next();

	if (fDouble){
		*type_res = HSPVAR_FLAG_DOUBLE;
		return (void *)&ref_doubleval;
	}

	if (fFloat) {
		*type_res = HSPVAR_FLAG_INT;
		return (void *)&ref_floatval;
	}

	*type_res = HSPVAR_FLAG_INT;
	return (void *)&ref_int32val;
}




static int cmdfunc(int cmd)
{
	code_next();

	switch (cmd) {

	case 0x01:
	{
		FluidProcess();
		break;
	}

	case 0x02:
	{
		GetSpeedPoints();
		break;
	}

	case 0x03:
	{
		poissonloop = code_getdi(8);
		break;
	}

	default:
		puterror(HSPERR_UNSUPPORTED_FUNCTION);
	}
	return RUNMODE_RUN;
}




















static void FluidProcess(void)
{
	PVal* pval1;
	APTR aptr1;	//�z��ϐ��̎擾
	aptr1 = code_getva(&pval1);//	���͕ϐ��̌^�Ǝ��̂̃|�C���^���擾
	HspVarProc* phvp1;
	void* ptr1;
	phvp1 = exinfo->HspFunc_getproc(pval1->flag);	//�^����������HspVarProc�\���̂ւ̃|�C���^
	ptr1 = phvp1->GetPtr(pval1);					//�f�[�^�ipval1�j�̎��Ԃ�����擪�|�C���^���擾�B

	PVal* pval2;
	APTR aptr2;	//�z��ϐ��̎擾
	aptr2 = code_getva(&pval2);//	���͕ϐ��̌^�Ǝ��̂̃|�C���^���擾
	HspVarProc* phvp2;
	void* ptr2;
	phvp2 = exinfo->HspFunc_getproc(pval2->flag);	//�^����������HspVarProc�\���̂ւ̃|�C���^
	ptr2 = phvp2->GetPtr(pval2);					//�f�[�^�ipval1�j�̎��Ԃ�����擪�|�C���^���擾�B

	PVal* pval3;
	APTR aptr3;	//�z��ϐ��̎擾
	aptr3 = code_getva(&pval3);//	���͕ϐ��̌^�Ǝ��̂̃|�C���^���擾
	HspVarProc* phvp3;
	void* ptr3;
	phvp3 = exinfo->HspFunc_getproc(pval3->flag);	//�^����������HspVarProc�\���̂ւ̃|�C���^
	ptr3 = phvp3->GetPtr(pval3);

	PVal* pval4;
	APTR aptr4;	//�z��ϐ��̎擾
	aptr4 = code_getva(&pval4);//	���͕ϐ��̌^�Ǝ��̂̃|�C���^���擾
	HspVarProc* phvp4;
	void* ptr4;
	phvp4 = exinfo->HspFunc_getproc(pval4->flag);	//�^����������HspVarProc�\���̂ւ̃|�C���^
	ptr4 = phvp4->GetPtr(pval4);


	//�^�`�F�b�N
	if (pval1->flag != HSPVAR_FLAG_DOUBLE) 
	{
		MessageBox(NULL, "FluidProcess�̑�P������double�^�ł͂���܂���", "�G���[", 0);
		puterror(HSPERR_UNSUPPORTED_FUNCTION);
		return;
	}
	if (pval2->flag != HSPVAR_FLAG_DOUBLE)
	{
		MessageBox(NULL, "FluidProcess�̑�Q������double�^�ł͂���܂���", "�G���[", 0);
		puterror(HSPERR_UNSUPPORTED_FUNCTION);
		return;
	}
	if (pval3->flag != HSPVAR_FLAG_DOUBLE)
	{
		MessageBox(NULL, "FluidProcess�̑�R������double�^�ł͂���܂���", "�G���[", 0);
		puterror(HSPERR_UNSUPPORTED_FUNCTION);
		return;
	}
	if (pval4->flag != HSPVAR_FLAG_INT)
	{
		MessageBox(NULL, "FluidProcess�̑�S������int�^�ł͂���܂���", "�G���[", 0);
		puterror(HSPERR_UNSUPPORTED_FUNCTION);
		return;
	}


	int lengthx = pval1->len[1];
	int lengthy = pval1->len[2];

	//�z��v�f���`�F�b�N
	if ((lengthx != pval2->len[1]) | (lengthy != pval2->len[2]))
	{
		MessageBox(NULL, "FluidProcess�̑�P�����Ƒ�Q�����̔z��v�f�����قȂ�܂�", "�G���[", 0);
		puterror(HSPERR_UNSUPPORTED_FUNCTION);
		return;
	}
	if ((lengthx != pval3->len[1]) | (lengthy != pval3->len[2]))
	{
		MessageBox(NULL, "FluidProcess�̑�P�����Ƒ�R�����̔z��v�f�����قȂ�܂�", "�G���[", 0);
		puterror(HSPERR_UNSUPPORTED_FUNCTION);
		return;
	}
	if ((lengthx != pval4->len[1]) | (lengthy != pval4->len[2]))
	{
		MessageBox(NULL, "FluidProcess�̑�P�����Ƒ�S�����̔z��v�f�����قȂ�܂�", "�G���[", 0);
		puterror(HSPERR_UNSUPPORTED_FUNCTION);
		return;
	}



	int xysize = lengthx * lengthy;
	if (xysize < 4) 
	{
		MessageBox(NULL, "�z��̗v�f����4�����ł��B", "�G���[", 0);
		puterror(HSPERR_UNSUPPORTED_FUNCTION);
		return;
	}



	double* vx = (double*)ptr1;
	double* vy = (double*)ptr2;
	double* p = (double*)ptr3;
	int* wall = (int*)ptr4;

	
	if (bufsize < xysize) 
	{
		bufsize = xysize * 2;
		vx_after = new double[bufsize];
		vy_after = new double[bufsize];
		s = new double[bufsize];
	}

	

	double u, v;

	//�ǂ̎��O�v�Z
	unsigned int* pwall = new unsigned int[xysize];
	for (int y = 0; y < lengthy; y++)
	{
		for (int x = 0; x < lengthx; x++)
		{
			int ij = y * lengthx + x;
			int i0j = y * lengthx + (x - 1 + lengthx) % lengthx;
			int i1j = y * lengthx + (x + 1) % lengthx;
			int ij0 = (y - 1 + lengthy) % lengthy * lengthx + x;
			int ij1 = (y + 1) % lengthy * lengthx + x;

			unsigned int ans = 0;
			if (wall[i0j] == 1) {
				ans |= 1;
				ans |= 16;//x��
			}
			if (wall[i1j] == 1)ans |= 2;
			if (wall[ij0] == 1) {
				ans |= 4;
				ans |= 32;//y��
			}
			if (wall[ij1] == 1)ans |= 8;


			//x�ǂ�5bit��
			if (wall[ij] == 1)
				ans |= 16;
			//y�ǂ�6bit��
			if (wall[ij] == 1)
				ans |= 32;


			//�����̒l��2bit�ŁB7,8bit��
			ans += wall[ij] * 64;
			pwall[ij] = ans;
		}
	}




	//�ڗ�
	
	for (int y = 0; y < lengthy; y++)
	{
		for (int x = 0; x < lengthx; x++)
		{
			int ij = y * lengthx + x;
			int i0j = y * lengthx + (x - 1 + lengthx) % lengthx;
			int i1j = y * lengthx + (x + 1) % lengthx;
			int ij0 = (y - 1 + lengthy) % lengthy * lengthx + x;
			int ij1 = (y + 1) % lengthy * lengthx + x;
			int i0j1 = (y + 1) % lengthy * lengthx + (x - 1 + lengthx) % lengthx;
			vx_after[ij] = vx[ij];
			if (pwall[ij]&0b00000000000000000000000000010000)continue;

			u = vx[ij];
			v = (vy[i0j] + vy[ij] + vy[i0j1] + vy[ij1]) / 4;


			if (u >= 0.0) 
			{
				if (v >= 0.0) {
					vx_after[ij] = vx[ij] - u * (vx[ij] - vx[i0j]) - v * (vx[ij] - vx[ij0]);
				}
				else 
				{
					vx_after[ij] = vx[ij] - u * (vx[ij] - vx[i0j]) - v * (vx[ij1] - vx[ij]);
				}
			}
			else
			{
				if (v >= 0.0)
				{
					vx_after[ij] = vx[ij] - u * (vx[i1j] - vx[ij]) - v * (vx[ij] - vx[ij0]);
				}
				else 
				{
					vx_after[ij] = vx[ij] - u * (vx[i1j] - vx[ij]) - v * (vx[ij1] - vx[ij]);
				}
			}

			if (vx_after[ij] < -1.2) {
				vx_after[ij] = -1.2;
			}
			else 
			{
				if (vx_after[ij] > 1.2)vx_after[ij] = 1.2;
			}
			/*
			vx_after[ij] -= ((vx[i1j] + vx[ij]) * (vx[i1j] + vx[ij]) - (vx[ij] + vx[i0j]) * (vx[ij] + vx[i0j])
				+ (vx[ij1] + vx[ij]) * (vy[ij1] + vy[i0j1]) - (vx[ij0] + vx[ij]) * (vy[ij] + vy[i0j])
				) * 0.25;
				*/
		}
	}


	for (int y = 0; y < lengthy; y++)
	{
		for (int x = 0; x < lengthx; x++)
		{
			int ij = y * lengthx + x;
			int i0j = y * lengthx + (x - 1 + lengthx) % lengthx;
			int i1j = y * lengthx + (x + 1) % lengthx;
			int ij0 = (y - 1 + lengthy) % lengthy * lengthx + x;
			int ij1 = (y + 1) % lengthy * lengthx + x;
			int i1j0 = (y - 1 + lengthy) % lengthy * lengthx + (x + 1) % lengthx;

			vy_after[ij] = vy[ij];
			if (pwall[ij] & 0b00000000000000000000000000100000)continue;

			u = (vx[ij0] + vx[i1j0] + vx[i1j] + vx[ij]) / 4;
			v = vy[ij];


			if ((u >= 0.0)&(v >= 0.0)) {
				vy_after[ij] = vy[ij] - u * (vy[ij] - vy[i0j]) - v * (vy[ij] - vy[ij0]);
			}

			if ((u < 0.0)&(v >= 0.0)) {
				vy_after[ij] = vy[ij] - u * (vy[i1j] - vy[ij]) - v * (vy[ij] - vy[ij0]);
			}

			if ((u >= 0.0)&(v < 0.0)) {
				vy_after[ij] = vy[ij] - u * (vy[ij] - vy[i0j]) - v * (vy[ij1] - vy[ij]);
			}

			if ((u < 0.0)&(v < 0.0)) {
				vy_after[ij] = vy[ij] - u * (vy[i1j] - vy[ij]) - v * (vy[ij1] - vy[ij]);
			}
			

			if (vy_after[ij] < -1.2) {
				vy_after[ij] = -1.2;
			}
			else
			{
				if (vy_after[ij] > 1.2)vy_after[ij] = 1.2;
			}

			/*
			vy_after[ij] -= ((vy[ij1] + vy[ij]) * (vy[ij1] + vy[ij]) - (vy[ij] + vy[ij0]) * (vy[ij] + vy[ij0])
				+ (vx[i1j0] + vx[i1j]) * (vy[ij] + vy[i1j])- (vx[ij] + vx[ij0]) * (vy[ij] + vy[i0j])
				) * 0.25;
				*/
		}
	}


	//���U�v�Z
	for (int y = 0; y < lengthy; y++)
	{
		for (int x = 0; x < lengthx; x++)
		{
			int ij = y * lengthx + x;
			int i1j = y * lengthx + (x + 1) % lengthx;
			int ij1 = ((y + 1) % lengthy) * lengthx + x;
			s[ij] = (-vx_after[ij] - vy_after[ij] + vx_after[i1j] + vy_after[ij1]);
		}
	}
	
	
	
	//�|�A�\��������
	double pr1, pr2, pr3, pr4;
	for (int lpcnt = 0; lpcnt < poissonloop; lpcnt++)
	{
		for (int y = 0; y < lengthy; y++)
		{
			for (int x = 0; x < lengthx; x++)
			{
				int ij = y * lengthx + x;

				if (pwall[ij] >= 64)continue;

				int i0j = y * lengthx + (x - 1 + lengthx) % lengthx;
				int i1j = y * lengthx + (x + 1) % lengthx;
				int ij0 = ((y - 1 + lengthy) % lengthy) * lengthx + x;
				int ij1 = (y + 1) % lengthy * lengthx + x;
				
				
				if (pwall[ij] & 1)
				{
					pr1 = p[ij];
				}
				else 
				{
					pr1 = p[i0j];
				}

				if (pwall[ij] & 2)
				{
					pr2 = p[ij];
				}
				else
				{
					pr2 = p[i1j];
				}

				if (pwall[ij] & 4)
				{
					pr3 = p[ij];
				}
				else
				{
					pr3 = p[ij0];
				}

				if (pwall[ij] & 8)
				{
					pr4 = p[ij];
				}
				else
				{
					pr4 = p[ij1];
				}

				p[ij] = (1.0 - 1.72) * p[ij] + 1.72 / 4 * (pr1 + pr2 + pr3 + pr4 - s[ij]);
			}
		}
	}
	

	//�C��
	for (int y = 0; y < lengthy; y++)
	{
		for (int x = 0; x < lengthx; x++)
		{
			int ij = y * lengthx + x;
			int i0j = y * lengthx + (x - 1 + lengthx) % lengthx;
			int ij0 = (y - 1 + lengthy) % lengthy * lengthx + x;

			if ((pwall[ij] & 16) == 0) 
			{
				vx[ij] = vx_after[ij] - (p[ij] - p[i0j]);
			}
			if ((pwall[ij] & 32) == 0)
			{
				vy[ij] = vy_after[ij] - (p[ij] - p[ij0]);
			}

		}
	}


}


















static void GetSpeedPoints(void)
{
	PVal* pval1;
	APTR aptr1;	//�z��ϐ��̎擾
	aptr1 = code_getva(&pval1);//	���͕ϐ��̌^�Ǝ��̂̃|�C���^���擾
	HspVarProc* phvp1;
	void* ptr1;
	phvp1 = exinfo->HspFunc_getproc(pval1->flag);	//�^����������HspVarProc�\���̂ւ̃|�C���^
	ptr1 = phvp1->GetPtr(pval1);					//�f�[�^�ipval1�j�̎��Ԃ�����擪�|�C���^���擾�B


	PVal* pval2;
	APTR aptr2;	//�z��ϐ��̎擾
	aptr2 = code_getva(&pval2);//	���͕ϐ��̌^�Ǝ��̂̃|�C���^���擾
	HspVarProc* phvp2;
	void* ptr2;
	phvp2 = exinfo->HspFunc_getproc(pval2->flag);	//�^����������HspVarProc�\���̂ւ̃|�C���^
	ptr2 = phvp2->GetPtr(pval2);					//�f�[�^�ipval1�j�̎��Ԃ�����擪�|�C���^���擾�B


	PVal* pval3;
	APTR aptr3;	//�z��ϐ��̎擾
	aptr3 = code_getva(&pval3);//	���͕ϐ��̌^�Ǝ��̂̃|�C���^���擾
	HspVarProc* phvp3;
	void* ptr3;
	phvp3 = exinfo->HspFunc_getproc(pval3->flag);	//�^����������HspVarProc�\���̂ւ̃|�C���^
	ptr3 = phvp3->GetPtr(pval3);

	PVal* pval4;
	APTR aptr4;	//�z��ϐ��̎擾
	aptr4 = code_getva(&pval4);//	���͕ϐ��̌^�Ǝ��̂̃|�C���^���擾
	HspVarProc* phvp4;
	void* ptr4;
	phvp4 = exinfo->HspFunc_getproc(pval4->flag);	//�^����������HspVarProc�\���̂ւ̃|�C���^
	ptr4 = phvp4->GetPtr(pval4);


	PVal* pval5;
	APTR aptr5;	//�z��ϐ��̎擾
	aptr5 = code_getva(&pval5);//	���͕ϐ��̌^�Ǝ��̂̃|�C���^���擾
	HspVarProc* phvp5;
	void* ptr5;
	phvp5 = exinfo->HspFunc_getproc(pval5->flag);	//�^����������HspVarProc�\���̂ւ̃|�C���^
	ptr5 = phvp5->GetPtr(pval5);

	PVal* pval6;
	APTR aptr6;	//�z��ϐ��̎擾
	aptr6 = code_getva(&pval6);//	���͕ϐ��̌^�Ǝ��̂̃|�C���^���擾
	HspVarProc* phvp6;
	void* ptr6;
	phvp6 = exinfo->HspFunc_getproc(pval6->flag);	//�^����������HspVarProc�\���̂ւ̃|�C���^
	ptr6 = phvp6->GetPtr(pval6);



	//�^�`�F�b�N
	if (pval1->flag != HSPVAR_FLAG_DOUBLE)
	{
		MessageBox(NULL, "��P������double�^�ł͂���܂���", "�G���[", 0);
		puterror(HSPERR_UNSUPPORTED_FUNCTION);
		return;
	}
	if (pval2->flag != HSPVAR_FLAG_DOUBLE)
	{
		MessageBox(NULL, "��Q������double�^�ł͂���܂���", "�G���[", 0);
		puterror(HSPERR_UNSUPPORTED_FUNCTION);
		return;
	}
	if (pval3->flag != HSPVAR_FLAG_DOUBLE)
	{
		MessageBox(NULL, "��R������double�^�ł͂���܂���", "�G���[", 0);
		puterror(HSPERR_UNSUPPORTED_FUNCTION);
		return;
	}
	if (pval4->flag != HSPVAR_FLAG_DOUBLE)
	{
		MessageBox(NULL, "��S������double�^�ł͂���܂���", "�G���[", 0);
		puterror(HSPERR_UNSUPPORTED_FUNCTION);
		return;
	}
	if (pval5->flag != HSPVAR_FLAG_DOUBLE)
	{
		MessageBox(NULL, "��T������double�^�ł͂���܂���", "�G���[", 0);
		puterror(HSPERR_UNSUPPORTED_FUNCTION);
		return;
	}
	if (pval6->flag != HSPVAR_FLAG_DOUBLE)
	{
		MessageBox(NULL, "��U������double�^�ł͂���܂���", "�G���[", 0);
		puterror(HSPERR_UNSUPPORTED_FUNCTION);
		return;
	}


	int lengthx = pval1->len[1];
	int lengthy = pval1->len[2];

	//�z��v�f���`�F�b�N
	if ((lengthx != pval2->len[1]) | (lengthy != pval2->len[2]))
	{
		MessageBox(NULL, "��P�����Ƒ�Q�����̔z��v�f�����قȂ�܂�", "�G���[", 0);
		puterror(HSPERR_UNSUPPORTED_FUNCTION);
		return;
	}

	
	if (lengthx * lengthy < 4)
	{
		MessageBox(NULL, "��P�����܂��͑�Q�����̔z��̗v�f����4�����ł��B", "�G���[", 0);
		puterror(HSPERR_UNSUPPORTED_FUNCTION);
		return;
	}



	//�z��v�f���`�F�b�N
	if ((pval4->len[1] != pval3->len[1]) | (pval4->len[2] != pval3->len[2]))
	{
		MessageBox(NULL, "��R�����Ƒ�S�����̔z��v�f�����قȂ�܂�", "�G���[", 0);
		puterror(HSPERR_UNSUPPORTED_FUNCTION);
		return;
	}
	if ((pval5->len[1] < pval3->len[1]) | (pval5->len[2] < pval3->len[2]))
	{
		MessageBox(NULL, "��R��������T�����̔z��v�f�������Ȃ��ł�", "�G���[", 0);
		puterror(HSPERR_UNSUPPORTED_FUNCTION);
		return;
	}
	if ((pval6->len[1] < pval3->len[1]) | (pval6->len[2] < pval3->len[2]))
	{
		MessageBox(NULL, "��R��������U�����̔z��v�f�������Ȃ��ł�", "�G���[", 0);
		puterror(HSPERR_UNSUPPORTED_FUNCTION);
		return;
	}





	double dlx = lengthx;
	double dly = lengthy;
	int pcount = pval3->len[1];

	double* vx = (double*)ptr1;
	double* vy = (double*)ptr2;
	double* posx = (double*)ptr3;
	double* posy = (double*)ptr4;
	double* out_speedx = (double*)ptr5;
	double* out_speedy = (double*)ptr6;


	for (int i = 0; i < pcount; i++)
	{
		if (!isfinite(posx[i]))break;
		if (!isfinite(posy[i]))break;
		//x�ɂ���
		double xxyyx = posx[i];
		double xxyyy = posy[i] - 0.5;
		if (xxyyx < 0.0) {
			xxyyx = 0.0;
		}
		else 
		{
			if (xxyyx >= dlx) 
			{
				xxyyx = dlx - 0.001;
			}
		}

		if (xxyyy < 0.0) {
			xxyyy = 0.0;
		}
		else
		{
			if (xxyyy >= dly)
			{
				xxyyy = dly - 0.001;
			}
		}

		int ixx = xxyyx;
		int iyy = xxyyy;
		double sxx = xxyyx - ixx;
		double syy = xxyyy - iyy;

		int im1 = ixx + 1;
		if (im1 == lengthx)im1 = 0;

		int jm1 = iyy + 1;
		if (jm1 == lengthy)jm1 = 0;

		iyy *= lengthx;
		jm1 *= lengthx;

		out_speedx[i] =
			(((1.0 - sxx) * vx[ixx + iyy] + sxx * vx[im1 + iyy]) * (1.0 - syy) + ((1.0 - sxx) * vx[ixx + jm1] + sxx * vx[im1 + jm1]) * syy);


		//y�ɂ���
		xxyyx = posx[i] - 0.5;
		xxyyy = posy[i];
		if (xxyyx < 0.0) {
			xxyyx = 0.0;
		}
		else
		{
			if (xxyyx >= dlx)
			{
				xxyyx = dlx - 0.001;
			}
		}

		if (xxyyy < 0.0) {
			xxyyy = 0.0;
		}
		else
		{
			if (xxyyy >= dly)
			{
				xxyyy = dly - 0.001;
			}
		}

		ixx = xxyyx;
		iyy = xxyyy;
		sxx = xxyyx - ixx;
		syy = xxyyy - iyy;

		im1 = ixx + 1;
		if (im1 == lengthx)im1 = 0;

		jm1 = iyy + 1;
		if (jm1 == lengthy)jm1 = 0;

		iyy *= lengthx;
		jm1 *= lengthx;


		out_speedy[i] =
			(((1.0 - sxx) * vy[ixx + iyy] + sxx * vy[im1 + iyy]) * (1.0 - syy) + ((1.0 - sxx) * vy[ixx + jm1] + sxx * vy[im1 + jm1]) * syy);

	}
	
}


















/*------------------------------------------------------------*/

static int termfunc( int option )
{
	//		�I������ (�A�v���P�[�V�����I�����ɌĂ΂�܂�)
	//
	return 0;
}


/*------------------------------------------------------------*/
/*
		interface
*/
/*------------------------------------------------------------*/

int WINAPI DllMain (HINSTANCE hInstance, DWORD fdwReason, PVOID pvReserved)
{
	//		DLL�G���g���[ (��������K�v�͂���܂���)
	//
	return TRUE;
}


EXPORT void WINAPI hsp3cmdinit( HSP3TYPEINFO *info )
{
	//		�v���O�C�������� (���s�E�I��������o�^���܂�)
	//
	hsp3sdk_init( info );		// SDK�̏�����(�ŏ��ɍs�Ȃ��ĉ�����)

	info->reffunc = reffunc;		// �Q�Ɗ֐�(reffunc)�̓o�^
	info->cmdfunc = cmdfunc;		// ���߂̓o�^
	info->termfunc = termfunc;		// �I���֐�(termfunc)�̓o�^
}

/*------------------------------------------------------------*/