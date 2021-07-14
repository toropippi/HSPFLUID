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


static double ref_doubleval;						// 返値のための変数
static float ref_floatval;						// 返値のための変数
static int ref_int32val;						// 返値のための変数

static void *reffunc( int *type_res, int cmd )
{
	//		関数・システム変数の実行処理 (値の参照時に呼ばれます)
	//
	//			'('で始まるかを調べる
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
	APTR aptr1;	//配列変数の取得
	aptr1 = code_getva(&pval1);//	入力変数の型と実体のポインタを取得
	HspVarProc* phvp1;
	void* ptr1;
	phvp1 = exinfo->HspFunc_getproc(pval1->flag);	//型を処理するHspVarProc構造体へのポインタ
	ptr1 = phvp1->GetPtr(pval1);					//データ（pval1）の実態がある先頭ポインタを取得。

	PVal* pval2;
	APTR aptr2;	//配列変数の取得
	aptr2 = code_getva(&pval2);//	入力変数の型と実体のポインタを取得
	HspVarProc* phvp2;
	void* ptr2;
	phvp2 = exinfo->HspFunc_getproc(pval2->flag);	//型を処理するHspVarProc構造体へのポインタ
	ptr2 = phvp2->GetPtr(pval2);					//データ（pval1）の実態がある先頭ポインタを取得。

	PVal* pval3;
	APTR aptr3;	//配列変数の取得
	aptr3 = code_getva(&pval3);//	入力変数の型と実体のポインタを取得
	HspVarProc* phvp3;
	void* ptr3;
	phvp3 = exinfo->HspFunc_getproc(pval3->flag);	//型を処理するHspVarProc構造体へのポインタ
	ptr3 = phvp3->GetPtr(pval3);

	PVal* pval4;
	APTR aptr4;	//配列変数の取得
	aptr4 = code_getva(&pval4);//	入力変数の型と実体のポインタを取得
	HspVarProc* phvp4;
	void* ptr4;
	phvp4 = exinfo->HspFunc_getproc(pval4->flag);	//型を処理するHspVarProc構造体へのポインタ
	ptr4 = phvp4->GetPtr(pval4);


	//型チェック
	if (pval1->flag != HSPVAR_FLAG_DOUBLE) 
	{
		MessageBox(NULL, "FluidProcessの第１引数がdouble型ではありません", "エラー", 0);
		puterror(HSPERR_UNSUPPORTED_FUNCTION);
		return;
	}
	if (pval2->flag != HSPVAR_FLAG_DOUBLE)
	{
		MessageBox(NULL, "FluidProcessの第２引数がdouble型ではありません", "エラー", 0);
		puterror(HSPERR_UNSUPPORTED_FUNCTION);
		return;
	}
	if (pval3->flag != HSPVAR_FLAG_DOUBLE)
	{
		MessageBox(NULL, "FluidProcessの第３引数がdouble型ではありません", "エラー", 0);
		puterror(HSPERR_UNSUPPORTED_FUNCTION);
		return;
	}
	if (pval4->flag != HSPVAR_FLAG_INT)
	{
		MessageBox(NULL, "FluidProcessの第４引数がint型ではありません", "エラー", 0);
		puterror(HSPERR_UNSUPPORTED_FUNCTION);
		return;
	}


	int lengthx = pval1->len[1];
	int lengthy = pval1->len[2];

	//配列要素数チェック
	if ((lengthx != pval2->len[1]) | (lengthy != pval2->len[2]))
	{
		MessageBox(NULL, "FluidProcessの第１引数と第２引数の配列要素数が異なります", "エラー", 0);
		puterror(HSPERR_UNSUPPORTED_FUNCTION);
		return;
	}
	if ((lengthx != pval3->len[1]) | (lengthy != pval3->len[2]))
	{
		MessageBox(NULL, "FluidProcessの第１引数と第３引数の配列要素数が異なります", "エラー", 0);
		puterror(HSPERR_UNSUPPORTED_FUNCTION);
		return;
	}
	if ((lengthx != pval4->len[1]) | (lengthy != pval4->len[2]))
	{
		MessageBox(NULL, "FluidProcessの第１引数と第４引数の配列要素数が異なります", "エラー", 0);
		puterror(HSPERR_UNSUPPORTED_FUNCTION);
		return;
	}



	int xysize = lengthx * lengthy;
	if (xysize < 4) 
	{
		MessageBox(NULL, "配列の要素数が4未満です。", "エラー", 0);
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

	//壁の事前計算
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
				ans |= 16;//x壁
			}
			if (wall[i1j] == 1)ans |= 2;
			if (wall[ij0] == 1) {
				ans |= 4;
				ans |= 32;//y壁
			}
			if (wall[ij1] == 1)ans |= 8;


			//x壁は5bit目
			if (wall[ij] == 1)
				ans |= 16;
			//y壁は6bit目
			if (wall[ij] == 1)
				ans |= 32;


			//自分の値を2bitで。7,8bit目
			ans += wall[ij] * 64;
			pwall[ij] = ans;
		}
	}




	//移流
	
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


	//発散計算
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
	
	
	
	//ポアソン方程式
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
	

	//修正
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
	APTR aptr1;	//配列変数の取得
	aptr1 = code_getva(&pval1);//	入力変数の型と実体のポインタを取得
	HspVarProc* phvp1;
	void* ptr1;
	phvp1 = exinfo->HspFunc_getproc(pval1->flag);	//型を処理するHspVarProc構造体へのポインタ
	ptr1 = phvp1->GetPtr(pval1);					//データ（pval1）の実態がある先頭ポインタを取得。


	PVal* pval2;
	APTR aptr2;	//配列変数の取得
	aptr2 = code_getva(&pval2);//	入力変数の型と実体のポインタを取得
	HspVarProc* phvp2;
	void* ptr2;
	phvp2 = exinfo->HspFunc_getproc(pval2->flag);	//型を処理するHspVarProc構造体へのポインタ
	ptr2 = phvp2->GetPtr(pval2);					//データ（pval1）の実態がある先頭ポインタを取得。


	PVal* pval3;
	APTR aptr3;	//配列変数の取得
	aptr3 = code_getva(&pval3);//	入力変数の型と実体のポインタを取得
	HspVarProc* phvp3;
	void* ptr3;
	phvp3 = exinfo->HspFunc_getproc(pval3->flag);	//型を処理するHspVarProc構造体へのポインタ
	ptr3 = phvp3->GetPtr(pval3);

	PVal* pval4;
	APTR aptr4;	//配列変数の取得
	aptr4 = code_getva(&pval4);//	入力変数の型と実体のポインタを取得
	HspVarProc* phvp4;
	void* ptr4;
	phvp4 = exinfo->HspFunc_getproc(pval4->flag);	//型を処理するHspVarProc構造体へのポインタ
	ptr4 = phvp4->GetPtr(pval4);


	PVal* pval5;
	APTR aptr5;	//配列変数の取得
	aptr5 = code_getva(&pval5);//	入力変数の型と実体のポインタを取得
	HspVarProc* phvp5;
	void* ptr5;
	phvp5 = exinfo->HspFunc_getproc(pval5->flag);	//型を処理するHspVarProc構造体へのポインタ
	ptr5 = phvp5->GetPtr(pval5);

	PVal* pval6;
	APTR aptr6;	//配列変数の取得
	aptr6 = code_getva(&pval6);//	入力変数の型と実体のポインタを取得
	HspVarProc* phvp6;
	void* ptr6;
	phvp6 = exinfo->HspFunc_getproc(pval6->flag);	//型を処理するHspVarProc構造体へのポインタ
	ptr6 = phvp6->GetPtr(pval6);



	//型チェック
	if (pval1->flag != HSPVAR_FLAG_DOUBLE)
	{
		MessageBox(NULL, "第１引数がdouble型ではありません", "エラー", 0);
		puterror(HSPERR_UNSUPPORTED_FUNCTION);
		return;
	}
	if (pval2->flag != HSPVAR_FLAG_DOUBLE)
	{
		MessageBox(NULL, "第２引数がdouble型ではありません", "エラー", 0);
		puterror(HSPERR_UNSUPPORTED_FUNCTION);
		return;
	}
	if (pval3->flag != HSPVAR_FLAG_DOUBLE)
	{
		MessageBox(NULL, "第３引数がdouble型ではありません", "エラー", 0);
		puterror(HSPERR_UNSUPPORTED_FUNCTION);
		return;
	}
	if (pval4->flag != HSPVAR_FLAG_DOUBLE)
	{
		MessageBox(NULL, "第４引数がdouble型ではありません", "エラー", 0);
		puterror(HSPERR_UNSUPPORTED_FUNCTION);
		return;
	}
	if (pval5->flag != HSPVAR_FLAG_DOUBLE)
	{
		MessageBox(NULL, "第５引数がdouble型ではありません", "エラー", 0);
		puterror(HSPERR_UNSUPPORTED_FUNCTION);
		return;
	}
	if (pval6->flag != HSPVAR_FLAG_DOUBLE)
	{
		MessageBox(NULL, "第６引数がdouble型ではありません", "エラー", 0);
		puterror(HSPERR_UNSUPPORTED_FUNCTION);
		return;
	}


	int lengthx = pval1->len[1];
	int lengthy = pval1->len[2];

	//配列要素数チェック
	if ((lengthx != pval2->len[1]) | (lengthy != pval2->len[2]))
	{
		MessageBox(NULL, "第１引数と第２引数の配列要素数が異なります", "エラー", 0);
		puterror(HSPERR_UNSUPPORTED_FUNCTION);
		return;
	}

	
	if (lengthx * lengthy < 4)
	{
		MessageBox(NULL, "第１引数または第２引数の配列の要素数が4未満です。", "エラー", 0);
		puterror(HSPERR_UNSUPPORTED_FUNCTION);
		return;
	}



	//配列要素数チェック
	if ((pval4->len[1] != pval3->len[1]) | (pval4->len[2] != pval3->len[2]))
	{
		MessageBox(NULL, "第３引数と第４引数の配列要素数が異なります", "エラー", 0);
		puterror(HSPERR_UNSUPPORTED_FUNCTION);
		return;
	}
	if ((pval5->len[1] < pval3->len[1]) | (pval5->len[2] < pval3->len[2]))
	{
		MessageBox(NULL, "第３引数より第５引数の配列要素数が少ないです", "エラー", 0);
		puterror(HSPERR_UNSUPPORTED_FUNCTION);
		return;
	}
	if ((pval6->len[1] < pval3->len[1]) | (pval6->len[2] < pval3->len[2]))
	{
		MessageBox(NULL, "第３引数より第６引数の配列要素数が少ないです", "エラー", 0);
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
		//xについて
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


		//yについて
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
	//		終了処理 (アプリケーション終了時に呼ばれます)
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
	//		DLLエントリー (何もする必要はありません)
	//
	return TRUE;
}


EXPORT void WINAPI hsp3cmdinit( HSP3TYPEINFO *info )
{
	//		プラグイン初期化 (実行・終了処理を登録します)
	//
	hsp3sdk_init( info );		// SDKの初期化(最初に行なって下さい)

	info->reffunc = reffunc;		// 参照関数(reffunc)の登録
	info->cmdfunc = cmdfunc;		// 命令の登録
	info->termfunc = termfunc;		// 終了関数(termfunc)の登録
}

/*------------------------------------------------------------*/