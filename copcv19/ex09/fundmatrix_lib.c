/*****************************************************************************/
/*                                                                           */
/*                   Copyright 08/2006 by Dr. Andres Bruhn                   */
/*     Faculty of Mathematics and Computer Science, Saarland University,     */
/*                           Saarbruecken, Germany.                          */
/*                                                                           */
/*****************************************************************************/

#ifndef FUNDMATRIX_LIB_INCLUDED
#define FUNDMATRIX_LIB_INCLUDED

/* ----------------------------------------------------------------------- */

void PA_trans_n(
	/************************************************************/
	double **B,		/* in  : matrix; only lower triangle needed; unchanged      */
	int n,			/* in  : B is n * n matrix                                  */
	double delta,   /* in  : machine precision: lowest delta with 1 + delta > 1 */
	double eps,		/* in  : absolute precision of the calculated eigenvalues   */
	double *lambda, /* out : eigenvalues, ordered with decreasing modulus       */
	double **v		/* out : orthonormal eigenvectors as columns                */
					/************************************************************/

)

/*
 Cyclic Jacobi method for determining the eigenvalues and eigenvectors
 of a symmetric matrix.
 Ref.:  H.R. Schwarz: Numerische Mathematik. Teubner, Stuttgart, 1988.
        pp. 243-246.
*/

{
	/*****************************************************/
	int i, j, p, q;		/* loop variables                                    */
	double sum;			/* for summing up                                    */
	double theta, t, c; /* intermediate results                              */
	double r, s, g, h;  /* intermediate results                              */
	int k;				/* for switching                                     */
	double help;		/* for switching                                     */
	double **a;			/* working copy of B                                 */
						/*****************************************************/

	/* ---- allocate storage ---- */
	// ALLOC_MATRIX (1, n+1, n+1, &a);
	a = new double *[n + 1];
	for (int i = 0; i < 10; i++)
		a[i] = new double[n + 1];

	/* ---- initializations ---- */
	/* a := B */
	for (i = 1; i <= n; i++)
		for (j = 1; j <= n; j++)
			a[i][j] = B[i][j];

	/* v := unit matrix */
	for (i = 1; i <= n; i++)
	{
		for (j = 1; j <= n; j++)
			v[i][j] = 0.0;
		v[i][i] = 1;
	}

	/* ---- loop ---- */

	do
	{
		/* check whether accuracy is reached */
		sum = 0.0;
		for (i = 2; i <= n; i++)
			for (j = 1; j <= i - 1; j++)
				sum = sum + a[i][j] * a[i][j];

		if (sum + sum > eps * eps)
		/* accuracy not reached yet, new cycle */
		{
			for (p = 1; p <= n - 1; p++)
				for (q = p + 1; q <= n; q++)
					if (fabs(a[q][p]) >= eps * eps)
					{
						theta = (a[q][q] - a[p][p]) / (2.0 * a[q][p]);
						t = 1.0;
						if (fabs(theta) > delta)
							t = 1.0 / (theta + theta / fabs(theta) *
												   sqrt(theta * theta + 1.0));
						c = 1.0 / sqrt(1.0 + t * t);
						s = c * t;
						r = s / (1.0 + c);
						a[p][p] = a[p][p] - t * a[q][p];
						a[q][q] = a[q][q] + t * a[q][p];
						a[q][p] = 0.0;
						for (j = 1; j <= p - 1; j++)
						{
							g = a[q][j] + r * a[p][j];
							h = a[p][j] - r * a[q][j];
							a[p][j] = a[p][j] - s * g;
							a[q][j] = a[q][j] + s * h;
						}
						for (i = p + 1; i <= q - 1; i++)
						{
							g = a[q][i] + r * a[i][p];
							h = a[i][p] - r * a[q][i];
							a[i][p] = a[i][p] - s * g;
							a[q][i] = a[q][i] + s * h;
						}
						for (i = q + 1; i <= n; i++)
						{
							g = a[i][q] + r * a[i][p];
							h = a[i][p] - r * a[i][q];
							a[i][p] = a[i][p] - s * g;
							a[i][q] = a[i][q] + s * h;
						}
						for (i = 1; i <= n; i++)
						{
							g = v[i][q] + r * v[i][p];
							h = v[i][p] - r * v[i][q];
							v[i][p] = v[i][p] - s * g;
							v[i][q] = v[i][q] + s * h;
						}
					} /* if */
		}			  /* if */
	}				  /* do */
	while (sum + sum > eps * eps);

	for (i = 1; i <= n; i++)
		lambda[i] = a[i][i];

	/* ---- order eigenvalues and eigenvectors ---- */

	for (i = 1; i <= n - 1; i++)
	{
		k = i;
		for (j = i + 1; j <= n; j++)
			if (fabs(lambda[j]) > fabs(lambda[k]))
				k = j;
		if (k != i)
		{
			/* switch eigenvalue i and k */
			help = lambda[k];
			lambda[k] = lambda[i];
			lambda[i] = help;
			/* switch eigenvector i and k */
			for (j = 1; j <= n; j++)
			{
				help = v[j][k];
				v[j][k] = v[j][i];
				v[j][i] = help;
			}
		} /* if */
	}	 /* for */

	/* ---- disallocate storage ---- */

	// FREE_MATRIX(1, n+1, n+1, a);
	for (int i = 0; i < 10; i++)
		delete[] a[i];
	delete[] a;

	return;
}

/*---------------------------------------------------------------------------*/

void multiply_matrix_vector(
	double **A,
	double *x,
	double *b,
	int dim)
{
	double *tmp; /* avoid problems when x=b */
	ALLOC_VECTOR_DOUBLE(1, 1 + dim, &tmp);

	for (int i = 1; i <= dim; i++)
	{
		tmp[i] = 0.0;

		for (int j = 0; j <= dim; j++)
		{
			tmp[i] += A[i][j] * x[j];
		}
	}

	for (int i = 1; i <= dim; i++)
	{
		b[i] = tmp[i];
	}

	FREE_VECTOR_DOUBLE(1, 1 + dim, tmp);
}

/*---------------------------------------------------------------------------*/

void multiply_matrix_matrix(
	double **A,
	double **B,
	double **C,
	int dim)
{
	double **tmp; /* avoid problems when A=C or B=C */
	ALLOC_MATRIX_DOUBLE(1, 1 + dim, 1 + dim, &tmp);

	for (int i = 1; i <= dim; i++)
		for (int k = 1; k <= dim; k++)
		{
			tmp[i][k] = 0.0;

			for (int j = 1; j <= dim; j++)
			{
				tmp[i][k] += A[i][j] * B[j][k];
			}
		}

	for (int i = 1; i <= dim; i++)
		for (int k = 1; k <= dim; k++)
		{
			C[i][k] = tmp[i][k];
		}

	FREE_MATRIX_DOUBLE(1, 1 + dim, 1 + dim, tmp);
}

/*---------------------------------------------------------------------------*/

void transpose_matrix(
	double **A,
	double **B,
	int dim)
{
	double **tmp; /* avoid problems when A=B */
	ALLOC_MATRIX_DOUBLE(1, 1 + dim, 1 + dim, &tmp);

	for (int i = 1; i <= dim; i++)
		for (int j = 1; j <= dim; j++)
		{
			tmp[i][j] = A[j][i];
		}

	for (int i = 1; i <= dim; i++)
		for (int j = 1; j <= dim; j++)
		{
			B[i][j] = tmp[i][j];
		}
	FREE_MATRIX_DOUBLE(1, 1 + dim, 1 + dim, tmp);
}

/*---------------------------------------------------------------------------*/

void compute_fundamental_matrix_TLS(
	/*****************************************************/
	float **u,		   /* in     : x-component of computed flow field       */
	float **v,		   /* in     : y-component of computed flow field       */
	int nx,			   /* in     : size in x-direction                      */
	int ny,			   /* in     : size in y-direction                      */
	int bx,			   /* in     : boundary size in x-direction             */
	int by,			   /* in     : boundary size in y-direction             */
	int normalisation, /* in     : use normalisation?                       */
	double **F		   /* out    : Fundamental matrix                       */
					   /*****************************************************/
)
{
	/********************************************************/
	int i, j, k, l;		   /* loop variables                                       */
	double ***m1;		   /* point in the first image (projective coordinates)    */
	double ***m2;		   /* point in the second image (projective coordinates)   */
	double *s;			   /* constraint vector for correspondence between m1/m2   */
	double **StS;		   /* S^tS -> TLS matrix for eigenvalue problem            */
	double *e_val;		   /* eigenvalues of S^tS                                  */
	double **e_vec;		   /* eigenvectors of S^tS                                 */
	double fnorm;		   /* normalisation factor                                 */
	double **T1;		   /* transformation matrix in first image                 */
	double **T2;		   /* transformation matrix in second image                */
	double **Fn;		   /* fundamental matrix for normalised coordinates        */
	double x1, x2, y1, y2; /* shift parameters                                     */
	double s1, s2;		   /* scaling parameters                                   */
						   /********************************************************/

	/* alloc memory */
	ALLOC_CUBIX_DOUBLE(2, nx + 2 * bx, ny + 2 * by, 3 + 1, &m1, &m2);
	ALLOC_VECTOR_DOUBLE(2, 9 + 1, &s, &e_val);
	ALLOC_MATRIX_DOUBLE(2, 9 + 1, 9 + 1, &StS, &e_vec);
	ALLOC_MATRIX_DOUBLE(3, 3 + 1, 3 + 1, &T1, &T2, &Fn);

	/* Initialise TLS matrix S^tS with zero */
	for (k = 1; k <= 9; k++)
		for (l = 1; l <= 9; l++)
		{
			StS[k][l] = 0.0;
		}

	/* Initialise transformation matrices with zero */
	for (k = 1; k <= 3; k++)
		for (l = 1; l <= 3; l++)
		{
			T1[k][l] = 0.0;
			T2[k][l] = 0.0;
		}

	/* Determine corresponding points */
	for (i = bx; i < nx + bx; i++)
		for (j = by; j < ny + by; j++)
		{
			/* point in left image (m1, projective coordinates) */
			/* --- TODO (Exercise 9.1): insert your code here -> compute m1 --- */
			m1[i][j][1] = i-bx;
			m1[i][j][2] = j-by;
			m1[i][j][3] = 1;

			/* !!! ATTENTION !!! 
			we allocate vectors with one additional entry so that their index runs 
			from 1,..,n */
			/* ------------------------------------------------- */

			/* corresponding point in right image (m2, projective coordinates) */
			/* --- TODO (Exercise 9.1): insert your code here -> compute m2 --- */
			m2[i][j][1] = m1[i][j][1]  + u[i][j];
			m2[i][j][2] = m1[i][j][2]  + v[i][j];
			m2[i][j][3] = 1;
			/* ------------------------------------------------- */
		}

	/* If normalisation is used, determine */
	/* transformation matrices T1 and T2 */
	if (normalisation)
	{
		/* compute shift parameters (mean x/y coordinate)*/
		x1 = x2 = y1 = y2 = 0.0;
		/* --- TODO (Exercise 9.3): insert your code here -> compute x1,x2,y1,y2 */
		for (i = bx; i < nx + bx; i++)
			for (j = by; j < ny + by; j++){
				x1 += 1/(nx*ny) * m1[i][j][1];
				y1 += 1/(nx*ny) * m1[i][j][2];
				x2 += 1/(nx*ny) * m2[i][j][1];
				y2 += 1/(nx*ny) * m2[i][j][2];
			}

		/* compute scaling parameters */
		s1 = s2 = 0.0;
		/* --- TODO (Exercise 9.3): insert your code here -> compute s1,s2 */
		for (i = bx; i < nx + bx; i++)
			for (j = by; j < ny + by; j++){
				s1 += 1/(nx*ny) * sqrt( pow( m1[i][j][1] - x1, 2) + pow( m1[i][j][2] - y1, 2) + 1 );
				s2 += 1/(nx*ny) * sqrt( pow( m2[i][j][1] - x2, 2) + pow( m2[i][j][2] - y2, 2) + 1 );
			}
		s1 = sqrt(3)/s1;
		s2 = sqrt(3)/s2;

		/* set up transformation matrices */
		T1[1][1] = T1[2][2] = s1;
		T2[1][1] = T2[2][2] = s2;

		T1[1][3] = -s1 * x1;
		T1[2][3] = -s1 * y1;

		T2[1][3] = -s2 * x2;
		T2[2][3] = -s2 * y2;

		T1[3][3] = T2[3][3] = 1.0;

		for (i = bx; i < nx + bx; i++)
			for (j = by; j < ny + by; j++)
			{
				/* --- TODO (Exercise 9.3): insert your code here -> compute transformed coordinates */
				multiply_matrix_vector( T1, m1[i][j], m1[i][j], 3 );
				multiply_matrix_vector( T2, m2[i][j], m2[i][j], 3 );
			}
	}

	/* Determine entries of S^tS matrix */
	for (i = bx; i < nx + bx; i++)
		for (j = by; j < ny + by; j++)
		{
			/* constraint vector s */
			/* --- TODO (Exercise 9.1): insert your code here -> compute the constraint s_i 
	             (later used for s_i^\top f) --- */
			/* ------------------------------------------------- */
			s[1] = m2[i][j][1] * m1[i][j][1];
			s[2] = m2[i][j][1] * m1[i][j][2];
			s[3] = m2[i][j][1] * m1[i][j][3];
			s[4] = m2[i][j][2] * m1[i][j][1];
			s[5] = m2[i][j][2] * m1[i][j][2];
			s[6] = m2[i][j][2] * m1[i][j][3];
			s[7] = m2[i][j][3] * m1[i][j][1];
			s[8] = m2[i][j][3] * m1[i][j][2];
			s[9] = m2[i][j][3] * m1[i][j][3];

			/* update matrix S^tS by new constraint s_i */
			for (k = 1; k <= 9; k++)
				for (l = 1; l <= 9; l++)
				{
					/* --- TODO (Exercise 9.1): insert your code here -> compute 
    				S^\top S + s_i s_i^\top --- */
					/* ------------------------------------------------- */
					StS[k][l] += s[k] * s[l];
				}
		}

	/* compute eigenvalues and eigenvectors */
	PA_trans_n(StS, 9, 1e-14, 1e-12, e_val, e_vec);

	for (k = 1; k <= 9; k++)
		printf("\n %f", e_val[k]);

	/* use smallest eigenvalue as solution */
	F[1][1] = e_vec[1][9];
	F[1][2] = e_vec[2][9];
	F[1][3] = e_vec[3][9];

	F[2][1] = e_vec[4][9];
	F[2][2] = e_vec[5][9];
	F[2][3] = e_vec[6][9];

	F[3][1] = e_vec[7][9];
	F[3][2] = e_vec[8][9];
	F[3][3] = e_vec[9][9];

	/* If normalisation is used, apply the */
	/* backtransformation to the fundamental matrix */
	if (normalisation)
	{
		/* --- TODO (Exercise 9.3): insert your code here -> compute */
		/* the backtransformed fundamental matrix */
		transpose_matrix( T2, T2, 3 );
		multiply_matrix_matrix( T2, F, F, 3 );
		multiply_matrix_matrix( F, T1, F, 3 );
	}

	/* normalise matrix */
	fnorm = F[2][3];

	for (k = 1; k <= 3; k++)
		for (l = 1; l <= 3; l++)
		{
			F[k][l] /= fnorm;
		}

	/* print out fundamental matrix */
	printf("\n\nThe estimated fundamental matrix reads:");
	printf("\n %lf %lf %lf", F[1][1], F[1][2], F[1][3]);
	printf("\n %lf %lf %lf", F[2][1], F[2][2], F[2][3]);
	printf("\n %lf %lf %lf", F[3][1], F[3][2], F[3][3]);
	printf("\n ");

	// /* free memory */
	FREE_CUBIX_DOUBLE(2, nx + 2 * bx, ny + 2 * by, 3 + 1, m1, m2);
	FREE_VECTOR_DOUBLE(2, 9 + 1, s, e_val);

	FREE_MATRIX_DOUBLE(2, 9 + 1, 9 + 1, StS, e_vec);
	FREE_MATRIX_DOUBLE(3, 3 + 1, 3 + 1, T1, T2, Fn);

	return;
}

/*--------------------------------------------------------------------------*/

#endif
