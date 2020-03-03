#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "diff_tensor.h"

/*--------------------------------------------------------------------------*/
/*                                                                          */
/*           COHERENCE-ENHANCING ANISOTROPIC DIFFUSION FILTERING            */
/*                                                                          */
/*                       (Joachim Weickert, 6/2000)                         */
/*                                                                          */
/*--------------------------------------------------------------------------*/

/* 
 features:
 - explicit scheme
 - presmoothing at noise scale:  convolution-based, Neumann b.c.
 - presmoothing at integration scale: convolution-based, Dirichlet b.c.
*/

/* ------------------------------------------------------------------------ */

void PA_trans

	(float a11,   /* coeffs of (2*2)-matrix */
	 float a12,   /* coeffs of (2*2)-matrix */
	 float a22,   /* coeffs of (2*2)-matrix */
	 float *c,	/* 1. comp. of 1. eigenvector, output */
	 float *s,	/* 2. comp. of 1. eigenvector, output */
	 float *lam1, /* larger  eigenvalue, output */
	 float *lam2) /* smaller eigenvalue, output */

/*
 Principal axis transformation. 
*/

{
	float help, norm; /* time savers */

	/* eigenvalues */
	help = sqrt(pow(a22 - a11, 2.0) + 4 * a12 * a12);
	*lam1 = (a11 + a22 + help) / 2.0;
	*lam2 = (a11 + a22 - help) / 2.0;

	/* eigenvectors */
	*c = 2.0 * a12;
	*s = a22 - a11 + help;

	/* normalized eigenvectors */
	norm = sqrt(*c * *c + *s * *s);
	if (norm >= 0.0000001)
	{
		*c = *c / norm;
		*s = *s / norm;
	}
	else
	{
		*c = 1.0;
		*s = 0.0;
	}
	return;

} /* PA_trans */

/* ----------------------------------------------------------------------- */

void PA_backtrans

	(float c,	/* 1. comp. of 1. eigenvector */
	 float s,	/* 2. comp. of 1. eigenvector */
	 float lam1, /* 1. eigenvalue */
	 float lam2, /* 2. eigenvalue */
	 float *a11, /* coeff. of (2*2)-matrix, output */
	 float *a12, /* coeff. of (2*2)-matrix, output */
	 float *a22) /* coeff. of (2*2)-matrix, output */

/*
  Principal axis backtransformation of a symmetric (2*2)-matrix. 
 A = U * diag(lam1, lam2) * U_transpose with U = (v1 | v2)     
 v1 = (c, s) is first eigenvector
*/

{

	*a11 = c * c * lam1 + s * s * lam2;
	*a22 = lam1 + lam2 - *a11; /* trace invariance */
	*a12 = c * s * (lam1 - lam2);

	return;

} /* PA_backtrans */

/*--------------------------------------------------------------------------*/

void diff_tensor

	(float C,	  /* contrast parameter */
	 float alpha, /* linear diffusion fraction */
	 long nx,	  /* image dimension in x direction */
	 long ny,	  /* image dimension in y direction */
	 float **dxx, /* in: structure tensor el., out: diff. tensor el. */
	 float **dxy, /* in: structure tensor el., out: diff. tensor el. */
	 float **dyy) /* in: structure tensor el., out: diff. tensor el. */

/*
 Calculates the diffusion tensor of CED by means of the structure tensor.
*/

{
	long i, j;		  /* loop variables */
	float beta;		  /* time saver */
	float c, s;		  /* specify first eigenvector */
	float mu1, mu2;   /* eigenvalues of structure tensor */
	float lam1, lam2; /* eigenvalues of diffusion tensor */

	beta = 1.0 - alpha;

	for (i = 1; i <= nx; i++)
		for (j = 1; j <= ny; j++)
		{
			/* principal axis transformation */
			PA_trans( dxx[i][j], dxy[i][j], dyy[i][j],
	 				  &c, &s, &mu1, &mu2 );

			/* calculate eigenvalues */
			lam1 = alpha;
			if (mu1==mu2)
				lam2 = alpha;
			else
				lam2 = alpha + (1-alpha) * exp( -C / pow(mu1-mu2,2) );

			/* principal axis backtransformation */
			PA_backtrans ( c, s, lam1, lam2,
						   &dxx[i][j], &dxy[i][j], &dyy[i][j]);
		}

	return;

} /* diff_tensor */

/*--------------------------------------------------------------------------*/
