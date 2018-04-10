#include <stdio.h> //printf
#include <stdlib.h>
#include "HLS/math.h"
#include "HLS/extendedmath.h"
#include "HLS/hls.h"

// NOTE: &a streams into this component in order
// Can do better in memory usage because we don't need to s and v for 9 * 8 SVD
component void r8mat_svd_linpack_9_8(short m, short n, ihc::stream_in<float> &a, ihc::stream_out<float> &U, ihc::stream_out<float> &S,
	ihc::stream_out<float> &V) 
{
	printf("start svd.....\n");
	// Hardwire this version to 9 * 8 svd 
	const short M = 9;
	const short N = 8;

	short lda;
	short ldu;
	short ldv;
	short job;
	short lwork;

	// Allocate memory
	float a_copy[M * N];
	float u[M * M];
	float s[M * N];
	float v[N * N];
	float e[M + N]; 
	float sdiag[M + N];
	float work[M];
	/*
	Compute the eigenvalues and eigenvectors.
	*/
	job = 11;
	lda = m;
	ldu = m;
	ldv = n;
	/*
	The input matrix is destroyed by the routine.  Since we need to keep
	it around, we only pass a copy to the routine.
	*/
	/*
	for (short j = 0; j < n; j++)
	{
		for (short i = 0; i < m; i++)
		{
			a_copy[i + j * m] = a[i + j * m];
		}
	}
	*/

	// copy in col major order
	for (short i = 0; i < m * n; i++){
		a_copy[i] = a.read();
	}

	/*
	printf("a_copy\n");
	for(short i = 0; i < 9; i++){
		for(short j = 0; j < 8; j++){
			printf("%f\t", a_copy[i  + j * 9]);
		}
		printf("\n");
	}
	printf("\n");
	*/
	//info = dsvdc(a_copy, lda, m, n, sdiag, e, u, ldu, v, ldv, work, job);
	
	float b;
	float c;
	float cs;
	float el;
	float emm1;
	float f;
	float g;
	short i;
	short info;
	short iter;
	short j;
	short jobu;
	short k;
	short kase;
	short kk;
	short l;
	short ll;
	short lls;
	short ls;
	short lu;
	short maxit = 30;
	short mm;
	short mm1;
	short mn;
	short mp1;
	short nct;
	short nctp1;
	short ncu;
	short nrt;
	short nrtp1;
	float scale;
	float shift;
	float sl;
	float sm;
	float smm1;
	float sn;
	float t;
	float t1;
	float test;
	short wantu;
	short wantv;
	float ztest;
	/*
	Determine what is to be computed.
	*/
	info = 0;
	wantu = 0;
	wantv = 0;
	jobu = (job % 100) / 10;

	if (1 < jobu)
	{
		//ncu = std::min(m, n);
		ncu = m < n ? m : n;
	}
	else
	{
		ncu = m;
	}

	if (jobu != 0)
	{
		wantu = 1;
	}

	if ((job % 10) != 0)
	{
		wantv = 1;
	}
	/*
	Reduce A to bidiagonal form, storing the diagonal elements
	in S and the super-diagonal elements in E.
	*/
	/*
	nct = std::min(m - 1, n);
	nrt = std::max(0, std::min(m, n - 2));
	lu = std::max(nct, nrt);
	//*/
	nct = m - 1 < n ? m - 1 : n;
	short tmp = m < n - 2 ? m : n -2;
	nrt = 0 > tmp ? 0 : tmp;
	lu = nct > nrt ? nct : nrt;
	for (l = 1; l <= lu; l++)
	{
		/*
		Compute the transformation for the L-th column and
		place the L-th diagonal in S(L).
		*/
		if (l <= nct)
		{
			//sdiag[l - 1] = dnrm2(m - l + 1, a_copy + l - 1 + (l - 1)*lda, 1);
			// replace dnrm2
			short size = m - l + 1;
			float* x = a_copy + l - 1 + (l - 1)*lda;
			if (size <= 0) {
				sdiag[l - 1] = 0.0;
			}
			else {
				float accum = 0;
				for (short i = 0; i < size; i++) {
					accum += x[i] * x[i];
				}
				sdiag[l - 1] = sqrt(accum);
			}
			///////////////////////////////////

			if (sdiag[l - 1] != 0.0)
			{
				if (a_copy[l - 1 + (l - 1)*lda] < 0.0)
				{
					sdiag[l - 1] = -sdiag[l - 1];
				}
				//dscal(m - l + 1, 1.0 / sdiag[l - 1], a_copy + l - 1 + (l - 1)*lda, 1);
				// Replace dscal
				short size = m - l + 1;
				float sa = 1.0 / sdiag[l - 1];
				float* x = a_copy + l - 1 + (l - 1)*lda;
				if (size > 0) {
					for (short i = 0; i < size; i++) {
						x[i] = sa * x[i];
					}
				}
				/////////////////////////////

				a_copy[l - 1 + (l - 1)*lda] = 1.0 + a_copy[l - 1 + (l - 1)*lda];
			}
			sdiag[l - 1] = -sdiag[l - 1];
		}

		for (j = l + 1; j <= n; j++)
		{
			/*
			Apply the transformation.
			*/
			if (l <= nct && sdiag[l - 1] != 0.0)
			{
				//t = -ddot(m - l + 1, a_copy + l - 1 + (l - 1)*lda, 1, a_copy + l - 1 + (j - 1)*lda, 1)
				//	/ a_copy[l - 1 + (l - 1)*lda];
				// replace: ddot
				short size = m - l + 1;
				float* x = a_copy + l - 1 + (l - 1)*lda;
				float* y = a_copy + l - 1 + (j - 1)*lda;
				t = 0;
				for (short i = 0; i < size; i++) {
					t += x[i] * y[i];
				}
				t = -t / a_copy[l - 1 + (l - 1)*lda];
				////////////////////////////////////

				//daxpy(m - l + 1, t, a_copy + l - 1 + (l - 1)*lda, 1, a_copy + l - 1 + (j - 1)*lda, 1);
				// replace daxpy
				short size_m = m - l + 1;
				float da_m = t;
				float *dx_m = a_copy + l - 1 + (l - 1)*lda;
				float *dy_m = a_copy + l - 1 + (j - 1)*lda;
				for (short i = 0; i < size_m; i++) {
					dy_m[i] = da_m * dx_m[i] + dy_m[i];
				}
				/////////////////////////////////////////////////
			}
			/*
			Place the L-th row of a_copy into E for the
			subsequent calculation of the row transformation.
			*/
			e[j - 1] = a_copy[l - 1 + (j - 1)*lda];
		}
		/*
		Place the transformation in U for subsequent back multiplication.
		*/
		if (wantu && l <= nct)
		{
			for (i = l; i <= m; i++)
			{
				u[i - 1 + (l - 1)*ldu] = a_copy[i - 1 + (l - 1)*lda];
			}
		}

		if (l <= nrt)
		{
			/*
			Compute the L-th row transformation and place the
			L-th superdiagonal in E(L).
			*/
			//e[l - 1] = dnrm2(n - l, e + l, 1);
			// replace dnrm2
			short size = n - l;
			float* x = e + l;
			if (size <= 0) {
				e[l - 1] = 0.0;
			}
			else {
				float accum = 0;
				for (short i = 0; i < size; i++) {
					accum += x[i] * x[i];
				}
				e[l - 1] = sqrt(accum);
			}
			///////////////////////////////////

			if (e[l - 1] != 0.0)
			{
				if (e[l] < 0.0)
				{
					e[l - 1] = -e[l - 1];
				}
				//dscal(n - l, 1.0 / e[l - 1], e + l, 1);
				// Replace dscal
				short size = n - l;
				float sa = 1.0 / e[l - 1];
				float* x = e + l;
				if (size > 0) {
					for (short i = 0; i < size; i++) {
						x[i] = sa * x[i];
					}
				}
				/////////////////////////////

				e[l] = 1.0 + e[l];
			}

			e[l - 1] = -e[l - 1];
			/*
			Apply the transformation.
			*/
			if (l + 1 <= m && e[l - 1] != 0.0)
			{
				for (j = l + 1; j <= m; j++)
				{
					work[j - 1] = 0.0;
				}

				for (j = l + 1; j <= n; j++)
				{
					//daxpy(m - l, e[j - 1], a_copy + l + (j - 1)*lda, 1, work + l, 1);
					// replace daxpy
					short size_m = m - l;
					float da_m = e[j - 1];
					float *dx_m = a_copy + l + (j - 1)*lda;
					float *dy_m = work + l;
					for (short i = 0; i < size_m; i++) {
						dy_m[i] = da_m * dx_m[i] + dy_m[i];
					}
					/////////////////////////////////////////////////
				}

				for (j = l + 1; j <= n; j++)
				{
					//daxpy(m - l, -e[j - 1] / e[l], work + l, 1, a_copy + l + (j - 1)*lda, 1);
					// replace daxpy
					short size_m = m - l;
					float da_m = -e[j - 1] / e[l];
					float *dx_m = work + l;
					float *dy_m = a_copy + l + (j - 1)*lda;
					for (short i = 0; i < size_m; i++) {
						dy_m[i] = da_m * dx_m[i] + dy_m[i];
					}
					/////////////////////////////////////////////////
				}
			}
			/*
			Place the transformation in V for subsequent back multiplication.
			*/
			if (wantv)
			{
				for (j = l + 1; j <= n; j++)
				{
					v[j - 1 + (l - 1)*ldv] = e[j - 1];
				}
			}
		}
	}
	/*
	Set up the final bidiagonal matrix of order MN.
	*/
	//mn = std::min(m + 1, n);
	mn = m + 1 < n ? m + 1 : n;
	nctp1 = nct + 1;
	nrtp1 = nrt + 1;

	if (nct < n)
	{
		sdiag[nctp1 - 1] = a_copy[nctp1 - 1 + (nctp1 - 1)*lda];
	}

	if (m < mn)
	{
		sdiag[mn - 1] = 0.0;
	}

	if (nrtp1 < mn)
	{
		e[nrtp1 - 1] = a_copy[nrtp1 - 1 + (mn - 1)*lda];
	}

	e[mn - 1] = 0.0;
	/*
	If required, generate U.
	*/
	if (wantu)
	{
		for (i = 1; i <= m; i++)
		{
			for (j = nctp1; j <= ncu; j++)
			{
				u[(i - 1) + (j - 1)*ldu] = 0.0;
			}
		}

		for (j = nctp1; j <= ncu; j++)
		{
			u[j - 1 + (j - 1)*ldu] = 1.0;
		}

		for (ll = 1; ll <= nct; ll++)
		{
			l = nct - ll + 1;

			if (sdiag[l - 1] != 0.0)
			{
				for (j = l + 1; j <= ncu; j++)
				{
					//t = -ddot(m - l + 1, u + (l - 1) + (l - 1)*ldu, 1, u + (l - 1) + (j - 1)*ldu, 1)
					//	/ u[l - 1 + (l - 1)*ldu];
					// replace: ddot
					short size = m - l + 1;
					float* x = u + (l - 1) + (l - 1)*ldu;
					float* y = u + (l - 1) + (j - 1)*ldu;
					t = 0;
					for (short i = 0; i < size; i++) {
						t += x[i] * y[i];
					}
					t = -t / u[l - 1 + (l - 1)*ldu];
					////////////////////////////////////
					//daxpy(m - l + 1, t, u + (l - 1) + (l - 1)*ldu, 1, u + (l - 1) + (j - 1)*ldu, 1);
					// replace daxpy
					short size_m = m - l + 1;
					float da_m = t;
					float *dx_m = u + (l - 1) + (l - 1)*ldu;
					float *dy_m = u + (l - 1) + (j - 1)*ldu;
					for (short i = 0; i < size_m; i++) {
						dy_m[i] = da_m * dx_m[i] + dy_m[i];
					}
					/////////////////////////////////////////////////
				}

				//dscal(m - l + 1, -1.0, u + (l - 1) + (l - 1)*ldu, 1);
				// Replace dscal
				short size = m - l + 1;
				float sa = -1.0;
				float* x = u + (l - 1) + (l - 1)*ldu;
				if (size > 0) {
					for (short i = 0; i < size; i++) {
						x[i] = sa * x[i];
					}
				}
				/////////////////////////////
				u[l - 1 + (l - 1)*ldu] = 1.0 + u[l - 1 + (l - 1)*ldu];
				for (i = 1; i <= l - 1; i++)
				{
					u[i - 1 + (l - 1)*ldu] = 0.0;
				}
			}
			else
			{
				for (i = 1; i <= m; i++)
				{
					u[i - 1 + (l - 1)*ldu] = 0.0;
				}
				u[l - 1 + (l - 1)*ldu] = 1.0;
			}
		}
	}
	/*
	If it is required, generate V.
	*/
	if (wantv)
	{
		for (ll = 1; ll <= n; ll++)
		{
			l = n - ll + 1;

			if (l <= nrt && e[l - 1] != 0.0)
			{
				for (j = l + 1; j <= n; j++)
				{
					//t = -ddot(n - l, v + l + (l - 1)*ldv, 1, v + l + (j - 1)*ldv, 1)
					//	/ v[l + (l - 1)*ldv];
					// replace: ddot
					short size = n - l;
					float* x = v + l + (l - 1)*ldv;
					float* y = v + l + (j - 1)*ldv;
					t = 0;
					for (short i = 0; i < size; i++) {
						t += x[i] * y[i];
					}
					t = -t / v[l + (l - 1)*ldv];
					////////////////////////////////////
					//daxpy(n - l, t, v + l + (l - 1)*ldv, 1, v + l + (j - 1)*ldv, 1);
					// replace daxpy
					short size_m = n - l;
					float da_m = t;
					float *dx_m = v + l + (l - 1)*ldv ;
					float *dy_m = v + l + (j - 1)*ldv;
					for (short i = 0; i < size_m; i++) {
						dy_m[i] = da_m * dx_m[i] + dy_m[i];
					}
					/////////////////////////////////////////////////
				}

			}
			for (i = 1; i <= n; i++)
			{
				v[i - 1 + (l - 1)*ldv] = 0.0;
			}
			v[l - 1 + (l - 1)*ldv] = 1.0;
		}
	}
	/*
	Main iteration loop for the singular values.
	*/
	mm = mn;
	iter = 0;

	while (0 < mn)
	{
		/*
		If too many iterations have been performed, set flag and return.
		*/
		if (maxit <= iter)
		{
			info = mn;
			break;
			//return info;
		}
		/*
		This section of the program inspects for
		negligible elements in the sdiag and E arrays.

		On completion the variables KASE and L are set as follows:

		KASE = 1     if sdiag(MN) and E(L-1) are negligible and L < MN
		KASE = 2     if sdiag(L) is negligible and L < MN
		KASE = 3     if E(L-1) is negligible, L < MN, and
		sdiag(L), ..., sdiag(MN) are not negligible (QR step).
		KASE = 4     if E(MN-1) is negligible (convergence).
		*/
		for (ll = 1; ll <= mn; ll++)
		{
			l = mn - ll;

			if (l == 0)
			{
				break;
			}

			test = fabs(sdiag[l - 1]) + fabs(sdiag[l]);
			ztest = test + fabs(e[l - 1]);

			if (ztest == test)
			{
				e[l - 1] = 0.0;
				break;
			}
		}

		if (l == mn - 1)
		{
			kase = 4;
		}
		else
		{
			mp1 = mn + 1;

			for (lls = l + 1; lls <= mn + 1; lls++)
			{
				ls = mn - lls + l + 1;

				if (ls == l)
				{
					break;
				}

				test = 0.0;
				if (ls != mn)
				{
					test = test + fabs(e[ls - 1]);
				}

				if (ls != l + 1)
				{
					test = test + fabs(e[ls - 2]);
				}

				ztest = test + fabs(sdiag[ls - 1]);

				if (ztest == test)
				{
					sdiag[ls - 1] = 0.0;
					break;
				}

			}

			if (ls == l)
			{
				kase = 3;
			}
			else if (ls == mn)
			{
				kase = 1;
			}
			else
			{
				kase = 2;
				l = ls;
			}
		}

		l = l + 1;
		/*
		Deflate negligible sdiag(MN).
		*/
		if (kase == 1)
		{
			mm1 = mn - 1;
			f = e[mn - 2];
			e[mn - 2] = 0.0;

			for (kk = 1; kk <= mm1; kk++)
			{
				k = mm1 - kk + l;
				t1 = sdiag[k - 1];
				//drotg(&t1, &f, &cs, &sn);
				// replace drotg
				float* sa = &t1;
				float* sb = &f;
				//float* c = &cs;
				//float* sdiag = &sn;;
				float r;
				float roe;
				float tmpScale;
				float z;
				if (fabs(*sb) < fabs(*sa))
				{
					roe = *sa;
				}
				else
				{
					roe = *sb;
				}

				tmpScale = fabs(*sa) + fabs(*sb);

				if (tmpScale == 0.0)
				{
					cs = 1.0;
					sn = 0.0;
					r = 0.0;
				}
				else
				{
					r = tmpScale * sqrt((*sa / tmpScale) * (*sa / tmpScale)
						+ (*sb / tmpScale) * (*sb / tmpScale));
					printf("To be sqrted: %f", (*sa / tmpScale) * (*sa / tmpScale)
						+ (*sb / tmpScale) * (*sb / tmpScale));
					printf("r is %f\n", r);
					r = roe >= 0 ? r : -r;
					cs = *sa / r;
					sn = *sb / r;
				}

				if (0.0 < fabs(cs) && fabs(cs) <= sn)
				{
					z = 1.0 / cs;
				}
				else
				{
					z = sn;
				}

				*sa = r;
				*sb = z;
				///////////////////////////////
				sdiag[k - 1] = t1;

				if (k != l)
				{
					f = -sn * e[k - 2];
					e[k - 2] = cs * e[k - 2];
				}

				if (wantv)
				{
					//drot(n, v + 0 + (k - 1)*ldv, 1, v + 0 + (mn - 1)*ldv, 1, cs, sn);
					// replace drot
					short size = n;
					float* x = v + 0 + (k - 1)*ldv;
					float* y = v + 0 + (mn - 1)*ldv;
					//c = cs;
					//sdiag = sn;
					float stemp;
					if (size <= 0) {}
					else {
						for (short i = 0; i < size; i++) {
							stemp = cs * x[i] + sn * y[i];
							y[i] = cs * y[i] - sn * x[i];
							x[i] = stemp;
						}
					}
					///////////////////////////////

				}
			}
		}
		/*
		Split at negligible sdiag(L).
		*/
		else if (kase == 2)
		{
			f = e[l - 2];
			e[l - 2] = 0.0;

			for (k = l; k <= mn; k++)
			{
				t1 = sdiag[k - 1];
				//drotg(&t1, &f, &cs, &sn);
				// replace drotg
				float* sa = &t1;
				float* sb = &f;
				//float* c = &cs;
				//float* sdiag = &sn;;
				float r;
				float roe;
				float tmpScale;
				float z;
				if (fabs(*sb) < fabs(*sa))
				{
					roe = *sa;
				}
				else
				{
					roe = *sb;
				}

				tmpScale = fabs(*sa) + fabs(*sb);

				if (tmpScale == 0.0)
				{
					cs = 1.0;
					sn = 0.0;
					r = 0.0;
				}
				else
				{
					r = tmpScale * sqrt((*sa / tmpScale) * (*sa / tmpScale)
						+ (*sb / tmpScale) * (*sb / tmpScale));
					printf("To be sqrted: %f", (*sa / tmpScale) * (*sa / tmpScale)
						+ (*sb / tmpScale) * (*sb / tmpScale));
					printf("r is %f\n", r);
					r = roe >= 0 ? r : -r;
					cs = *sa / r;
					sn = *sb / r;
				}

				if (0.0 < fabs(cs) && fabs(cs) <= sn)
				{
					z = 1.0 / cs;
				}
				else
				{
					z = sn;
				}

				*sa = r;
				*sb = z;
				////////////////////////////////////////////

				sdiag[k - 1] = t1;
				f = -sn * e[k - 1];
				e[k - 1] = cs * e[k - 1];
				if (wantu)
				{
					//drot(m, u + 0 + (k - 1)*ldu, 1, u + 0 + (l - 2)*ldu, 1, cs, sn);
					// replace drot
					short size = m;
					float* x = u + 0 + (k - 1)*ldu;
					float* y = u + 0 + (l - 2)*ldu;
					//c = cs;
					//sdiag = sn;
					float stemp;
					if (size <= 0) {}
					else {
						for (short i = 0; i < size; i++) {
							stemp = cs * x[i] + sn * y[i];
							y[i] = cs * y[i] - sn * x[i];
							x[i] = stemp;
						}
					}
					///////////////////////////////
				}
			}
		}
		/*
		Perform one QR step.
		*/
		else if (kase == 3)
		{
			/*
			Calculate the shift.
			*/
			scale = fabs(sdiag[mn - 1]);
			if (scale < fabs(sdiag[mn - 2]))
			{
				scale = fabs(sdiag[mn - 2]);
			}
			if (scale < fabs(e[mn - 2]))
			{
				scale = fabs(e[mn - 2]);
			}
			if (scale < fabs(sdiag[l - 1]))
			{
				scale = fabs(sdiag[l - 1]);
			}
			if (scale < fabs(e[l - 1]))
			{
				scale = fabs(e[l - 1]);
			}

			sm = sdiag[mn - 1] / scale;
			smm1 = sdiag[mn - 2] / scale;
			emm1 = e[mn - 2] / scale;
			sl = sdiag[l - 1] / scale;
			el = e[l - 1] / scale;
			b = ((smm1 + sm) * (smm1 - sm) + emm1 * emm1) / 2.0;
			c = (sm * emm1) * (sm * emm1);
			shift = 0.0;

			if (b != 0.0 || c != 0.0)
			{
				shift = sqrt(b * b + c);
				printf("shift is %f\n", shift);
				if (b < 0.0)
				{
					shift = -shift;
				}
				shift = c / (b + shift);
			}

			f = (sl + sm) * (sl - sm) - shift;
			g = sl * el;
			/*
			Chase zeros.
			*/
			mm1 = mn - 1;

			for (k = l; k <= mm1; k++)
			{
				//drotg(&f, &g, &cs, &sn);
				// replace drotg
				float* sa = &f;
				float* sb = &g;
				//float* c = &cs;
				//float* s = &sn;;
				float r;
				float roe;
				float tmpScale;
				float z;
				if (fabs(*sb) < fabs(*sa))
				{
					roe = *sa;
				}
				else
				{
					roe = *sb;
				}

				tmpScale = fabs(*sa) + fabs(*sb);

				if (tmpScale == 0.0)
				{
					cs = 1.0;
					sn = 0.0;
					r = 0.0;
				}
				else
				{
					r = tmpScale * sqrt((*sa / tmpScale) * (*sa / tmpScale)
						+ (*sb / tmpScale) * (*sb / tmpScale));
					printf("To be sqrted: %f", (*sa / tmpScale) * (*sa / tmpScale)
						+ (*sb / tmpScale) * (*sb / tmpScale));
					printf("r is %f\n", r);
					r = roe >= 0 ? r : -r;
					cs = *sa / r;
					sn = *sb / r;
				}

				if (0.0 < fabs(cs) && fabs(cs) <= sn)
				{
					z = 1.0 / cs;
				}
				else
				{
					z = sn;
				}

				*sa = r;
				*sb = z;
				/////////////////////////////////////////

				if (k != l)
				{
					e[k - 2] = f;
				}

				f = cs * sdiag[k - 1] + sn * e[k - 1];
				e[k - 1] = cs * e[k - 1] - sn * sdiag[k - 1];
				g = sn * sdiag[k];
				sdiag[k] = cs * sdiag[k];

				if (wantv)
				{
					//drot(n, v + 0 + (k - 1)*ldv, 1, v + 0 + k * ldv, 1, cs, sn);
					// replace drot
					short size = n;
					float* x = v + 0 + (k - 1)*ldv;
					float* y = v + 0 + k * ldv;
					//c = cs;
					//sdiag = sn;
					float stemp;
					if (size <= 0) {}
					else {
						for (short i = 0; i < size; i++) {
							stemp = cs * x[i] + sn * y[i];
							y[i] = cs * y[i] - sn * x[i];
							x[i] = stemp;
						}
					}
					///////////////////////////////
				}

				//drotg(&f, &g, &cs, &sn);
				// replace drotg
				sa = &f;
				sb = &g;
				//float* c = &cs;
				//float* s = &sn;;
				if (fabs(*sb) < fabs(*sa))
				{
					roe = *sa;
				}
				else
				{
					roe = *sb;
				}

				tmpScale = fabs(*sa) + fabs(*sb);

				if (tmpScale == 0.0)
				{
					cs = 1.0;
					sn = 0.0;
					r = 0.0;
				}
				else
				{
					r = tmpScale * sqrt((*sa / tmpScale) * (*sa / tmpScale)
						+ (*sb / tmpScale) * (*sb / tmpScale));
					printf("To be sqrted: %f", (*sa / tmpScale) * (*sa / tmpScale)
						+ (*sb / tmpScale) * (*sb / tmpScale));
					printf("r is %f\n", r);
					r = roe >= 0 ? r : -r;
					cs = *sa / r;
					sn = *sb / r;
				}

				if (0.0 < fabs(cs) && fabs(cs) <= sn)
				{
					z = 1.0 / cs;
				}
				else
				{
					z = sn;
				}

				*sa = r;
				*sb = z;
				/////////////////////////////////////////

				sdiag[k - 1] = f;
				f = cs * e[k - 1] + sn * sdiag[k];
				sdiag[k] = -sn * e[k - 1] + cs * sdiag[k];
				g = sn * e[k];
				e[k] = cs * e[k];

				if (wantu && k < m)
				{
					//drot(m, u + 0 + (k - 1)*ldu, 1, u + 0 + k * ldu, 1, cs, sn);
					// replace drot
					short size = m;
					float* x = u + 0 + (k - 1)*ldu;
					float* y = u + 0 + k * ldu;
					//c = cs;
					//sdiag = sn;
					float stemp;
					if (size <= 0) {}
					else {
						for (short i = 0; i < size; i++) {
							stemp = cs * x[i] + sn * y[i];
							y[i] = cs * y[i] - sn * x[i];
							x[i] = stemp;
						}
					}
					///////////////////////////////
				}
			}
			e[mn - 2] = f;
			iter = iter + 1;
		}
		/*
		Convergence.
		*/
		else if (kase == 4)
		{
			/*
			Make the singular value nonnegative.
			*/
			if (sdiag[l - 1] < 0.0)
			{
				sdiag[l - 1] = -sdiag[l - 1];
				if (wantv)
				{
					//dscal(n, -1.0, v + 0 + (l - 1)*ldv, 1);
					// Replace dscal
					short size = n;
					float sa = -1.0;
					float* x = v + 0 + (l - 1)*ldv;
					if (size > 0) {
						for (short i = 0; i < size; i++) {
							x[i] = sa * x[i];
						}
					}
					/////////////////////////////
				}
			}
			/*
			Order the singular value.
			*/
			for (; ; )
			{
				if (l == mm)
				{
					break;
				}

				if (sdiag[l] <= sdiag[l - 1])
				{
					break;
				}

				t = sdiag[l - 1];
				sdiag[l - 1] = sdiag[l];
				sdiag[l] = t;

				if (wantv && l < n)
				{
					//dswap(n, v + 0 + (l - 1)*ldv, 1, v + 0 + l * ldv, 1);
					// replace dswap
					short size = n;
					float* x = v + 0 + (l - 1)*ldv;
					float* y = v + 0 + l * ldv;
					for (short i = 0; i < size; i++) {
						float temp = y[i];
						y[i] = x[i];
						x[i] = temp;
					}
					////////////////////////////////////
				}

				if (wantu && l < m)
				{
					//dswap(m, u + 0 + (l - 1)*ldu, 1, u + 0 + l * ldu, 1);
					// replace dswap
					short size = m;
					float* x = u + 0 + (l - 1)*ldu;
					float* y = u + 0 + l * ldu;
					for (short i = 0; i < size; i++) {
						float temp = y[i];
						y[i] = x[i];
						x[i] = temp;
					}
					////////////////////////////////////
				}

				l = l + 1;
			}
			iter = 0;
			mn = mn - 1;
		}
	}


	///////////////////////////////////////////////////////////////////////
	if (info != 0)
	{
		printf("\n");
		printf("R8MAT_SVD_LINPACK - Failure!\n");
		printf("  The SVD could not be calculated.\n");
		printf("  LINPACK routine DSVDC returned a_copy nonzero\n");
		printf("  value of the error flag, INFO = %hd\n", info);
		return;
	}

	// In our case we don't need S for 9 * 8 decomposition
	/*
	Make the MxN matrix S from the diagonal values in SDIAG.
	*/
	
	for (j = 0; j < n; j++)
	{
		for (i = 0; i < m; i++)
		{
			if (i == j)
			{
				s[i + j * m] = sdiag[i];
			}
			else
			{
				s[i + j * m] = 0.0;
			}
		}
	}

	/*
    printf("U\n");
	for(short i = 0; i < 9; i++){
		for(short j = 0; j < 9; j++){
			printf("%f\t", u[i + j * 9]);
		}
		printf("\n");
	}
	printf("\n");

	printf("S\n");
	for(short i = 0; i < 9; i++){
		for(short j = 0; j < 8; j++){
			printf("%f\t", s[i + j * 9]);
		}
		printf("\n");
	}
	printf("\n");

	printf("V\n");
	for(short i = 0; i < 8; i++){
		for(short j = 0; j < 8; j++){
			printf("%f\t", v[i  + j * 8]);
		}
		printf("\n");
	}
	printf("\n");
	printf("Start streaming out\n");
	*/

	// Stream out U, S and V matrix in row major order
	// col major order would be straight forward U.wirte = u[i*N + j]
	for(short i = 0; i < M * M; i++){
		U.write(u[i]);
	}
	for(short i = 0; i < M * N; i++){
		S.write(s[i]);
	}
	for(short i = 0; i < N * N; i++){
		V.write(v[i]);
	}
	return;
}