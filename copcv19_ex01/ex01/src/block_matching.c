/*****************************************************************************/
/*                                                                           */
/*                   Copyright 08/2006 by Dr. Andres Bruhn                   */
/*     Faculty of Mathematics and Computer Science, Saarland University,     */
/*                           Saarbruecken, Germany.                          */
/*                                                                           */
/*****************************************************************************/

#ifndef OF_BLOCK_MATCHING_INCLUDED
#define OF_BLOCK_MATCHING_INCLUDED

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "alloc_mem_linear.c"
#include "alloc_mem_linear_mult.c"
#include "funct_lib.c"

/* ---------------------------------------------------------------------- */

void matching_cost

	(
		/********************************************/
		float **f1, /* in  :  image 1                           */
		float **f2, /* in  :  image 2                           */
		int i,		/* in  :  x-coordinate of window center     */
		int j,		/* in  :  y-coordinate of window center     */
		int n,		/* in  :  window size 2*n+1                 */
		int sx,		/* in  :  shift of second image in x-dir.   */
		int sy,		/* in  :  shift of second image in y-dir.   */
		float &cost /* out :  cost map                          */
		/********************************************/
	)

/* 
 computes matching cost for two windows
*/

{
	/* 
	compute SSD matching cost between a reference window of size (2*n+1)x(2*n+1) 
	in the first frame at location i,j and a candidate window at location 
	i+sx, i+sy in the second image
	*/

	cost = 0;
	for (int di = -n; di++; di <= n)
	{
		for (int dj = -n; dj++; dj <= n)
		{
			cost += pow( f1[i+di][j+dj] - f2[i+sx+di][j+sy+dj] ,2) ;
		}
	}

	return;

} /* matching_cost */

/* ---------------------------------------------------------------------- */

void compute_cost_map

	(
		/********************************************/
		float **f1,		   /* in     : 1st image                       */
		float **f2,		   /* in     : 2nd image                       */
		int i,			   /* in     : x-coordinate of window center   */
		int j,			   /* in     : y-coordinate of window center   */
		int m_window_size, /* in     : window size for cross-corr      */
		int m_max_disp,	/* in     : maximum displacement            */
		float **cost_map   /* out    : cost map for the pixel i,j      */
						   /********************************************/
	)

/* 
 Computes the block matching costs between neigbourhoods of a pixel i,j
in the first frame and all shifted neighbourhoods in the second image.
*/

{
	/**********************************************/
	int k, l; /* loop variables                             */
	/**********************************************/

	/* ---- compute cost map ---- */
	for (k = -m_max_disp; k <= m_max_disp; k++)
		for (l = -m_max_disp; l <= m_max_disp; l++)
			matching_cost(f1, f2, i, j, m_window_size, k, l,
						  cost_map[k + m_max_disp][l + m_max_disp]);

	return;

} /* compute_cost_map */

/* ---------------------------------------------------------------------- */

void subpixel_refinement

	(
		/********************************************/
		float c_1,   /* in  :  left neigbour of min cost         */
		float c0,	/* in  :  min cost                          */
		float c1,	/* in  :  right neighbour of min cost       */
		float &delta /* out :  subpixel correction               */
					 /********************************************/
	)

/*
computes subpixel correction 
*/

{
	delta = (c_1 - c1) / ( 2 * (c1 -2*c0 +c_1) );
	return;
}

/* ---------------------------------------------------------------------- */

void pick_min(
	/**********************************************/
	float **cost_map, /* in     : cost map for pixel i,j            */
	float &u,		  /* out    : x-component of displacement field */
	float &v,		  /* out    : y-component of displacement field */
	int nx,			  /* in     : size in x-direction               */
	int ny,			  /* in     : size in y-direction               */
	int bx,			  /* in     : boundary size in x-direction      */
	int by,			  /* in     : boundary size in y-direction      */
	int m_max_disp,   /* in     : maximum displacement              */
	int m_subpixel	/* in     : flag for subpixel refinement      */
	/**********************************************/
)

/* 
 Picks for a pixel the displacement with the smallest matching cost
*/

{
	/**********************************************/
	int k, l;				/* loop variables                             */
	float min_cost;			/* current minimal cost                       */
	int min_cost_disp_x;	/* corresponding displacement in x-direction  */
	int min_cost_disp_y;	/* corresponding displacement in y-direction  */
	float delta_x, delta_y; /* subpixel correction in x- and y-direction  */
							/**********************************************/

	/*  initialise minimum search */
	min_cost = cost_map[m_max_disp][m_max_disp];
	min_cost_disp_x = 0;
	min_cost_disp_y = 0;

	/* perform minimum search */
	for (k = -m_max_disp; k <= m_max_disp; k++)
		for (l = -m_max_disp; l <= m_max_disp; l++)
			if (cost_map[k + m_max_disp][l + m_max_disp] < min_cost)
			{
				/* save current minimal cost and corresponding displacement */
				min_cost = cost_map[k + m_max_disp][l + m_max_disp];
				min_cost_disp_x = k;
				min_cost_disp_y = l;
			}

	/* if no subpixel refinement is desired */
	if (m_subpixel == 0)
	{
		u = (float)min_cost_disp_x;
		v = (float)min_cost_disp_y;
	}

	/* if subpixel refinement is desired */
	if (m_subpixel == 1)
	{
		/* perform subpixel refinement in x-direction */
		/* (if maximum is not at cost map border) */
		if ((min_cost_disp_x > -m_max_disp) && (min_cost_disp_x < m_max_disp))
		{
			subpixel_refinement(
				cost_map[min_cost_disp_x + m_max_disp - 1][min_cost_disp_y + m_max_disp],
				cost_map[min_cost_disp_x + m_max_disp][min_cost_disp_y + m_max_disp],
				cost_map[min_cost_disp_x + m_max_disp + 1][min_cost_disp_y + m_max_disp],
				delta_x);

			u = (float)min_cost_disp_x + delta_x;
		}
		else
		{
			u = (float)min_cost_disp_x;
		}

		/* perform subpixel refinement in y-direction */
		/* (if maximum is not at cost map border) */
		if ((min_cost_disp_y > -m_max_disp) && (min_cost_disp_y < m_max_disp))
		{
			subpixel_refinement(
				cost_map[min_cost_disp_x + m_max_disp][min_cost_disp_y + m_max_disp - 1],
				cost_map[min_cost_disp_x + m_max_disp][min_cost_disp_y + m_max_disp],
				cost_map[min_cost_disp_x + m_max_disp][min_cost_disp_y + m_max_disp + 1],
				delta_y);

			v = (float)min_cost_disp_y + delta_y;
		}
		else
		{
			v = (float)min_cost_disp_y;
		}
	}

	return;
}

/* ------------------------------------------------------------------------- */

void BLOCK_MATCHING(
	/*****************************************************/
	float **f1,		   	/* in     : 1st image                                */
	float **f2,		   	/* in     : 2nd image                                */
	float **u,		   	/* out    : x-component of displacement field        */
	float **v,		   	/* out    : y-component of displacement field        */
	int nx,			   	/* in     : size in x-direction                      */
	int ny,			   	/* in     : size in y-direction                      */
	int bx,			   	/* in     : boundary size in x-direction             */
	int by,			   	/* in     : boundary size in y-direction             */
	int m_window_size, 	/* in     : window size for block matching           */
	int m_max_disp,		/* in     : maximum disparity                        */
	int m_subpixel,		/* in     : flag for subpixel refinement             */
	int m_occlusion,   	/* in     : flag for occlusion detection             */
	float m_occ_thres  	/* in     : error threshold for occlusion detection  */
	/*****************************************************/
)

/* computes displacement field with block matching */

{

	/*****************************************************/
	int i, j, k, l;			/* loop variables                                    */
	float **cost_map;		/* cost map for the pixel i,j                        */
	float **u_fw, **v_fw;   /* forward flow field                                */
	float **u_bw, **v_bw;   /* backward flow field                               */
	float delta_u, delta_v; /* differences of forward and backward flow field in */
							/* x- and y-direction (w.r.t. to u and v)            */
	float mag_delta_u_v;	/* magnitude of vector (delta u, delta v)^T          */
	int ii, jj;				/* corresponding positions of i and j in f2          */
	/*****************************************************/

	/* ---- alloc memory ---- */
	ALLOC_MATRIX(1, 2 * m_max_disp + 1, 2 * m_max_disp + 1, &cost_map);
	ALLOC_MATRIX(2, nx + 2 * bx, ny + 2 * by, &u_fw, &v_fw);
	ALLOC_MATRIX(2, nx + 2 * bx, ny + 2 * by, &u_bw, &v_bw);

	/* ---- if no occlusion detection is desired ---- */
	if (m_occlusion == 0)
	{
		/* initialise forward displacement field with zero */
		for (i = bx; i < bx + nx; i++)
			for (j = by; j < by + ny; j++)
			{
				u_fw[i][j] = 0;
				v_fw[i][j] = 0;
			}

		/* ---- for each pixel that is not to close to the boundary ---- */
		for (i = bx + m_window_size + m_max_disp; i < bx + nx - m_window_size - m_max_disp; i++)
			for (j = by + m_window_size + m_max_disp; j < by + ny - m_window_size - m_max_disp; j++)
			{
				/* compute the forward cost map for the pixels i,j */
				compute_cost_map(f1, f2, i, j, m_window_size, m_max_disp, cost_map);

				/* pick forward displacement with minimum matching cost */
				pick_min(cost_map, u_fw[i][j], v_fw[i][j], nx, ny, bx, by,
						 m_max_disp, m_subpixel);
			}
	}

	/* ---- if occlusion detection is desired ---- */
	if (m_occlusion == 1)
	{
		/* initialise forward and backward displacement field with zero */
		for (i = bx; i < bx + nx; i++)
			for (j = by; j < by + ny; j++)
			{
				u_fw[i][j] = 0;
				v_fw[i][j] = 0;
				u_bw[i][j] = 0;
				v_bw[i][j] = 0;
			}

		/* ---- for each pixel that is not too close to the boundary ---- */
		for (i = bx + m_window_size + m_max_disp; i < bx + nx - m_window_size - m_max_disp; i++)
			for (j = by + m_window_size + m_max_disp; j < by + ny - m_window_size - m_max_disp; j++)
			{
				/* compute the forward cost map for the pixels i,j */
				compute_cost_map(f1, f2, i, j, m_window_size, m_max_disp, cost_map);

				/* pick forward displacement with minimum matching cost */
				pick_min(cost_map, u_fw[i][j], v_fw[i][j], nx, ny, bx, by,
						 m_max_disp, m_subpixel);

				/* compute the backward cost map for the pixels i,j */
				compute_cost_map(f2, f1, i, j, m_window_size, m_max_disp, cost_map);

				/* pick backward displacement with minimum matching cost */
				pick_min(cost_map, u_bw[i][j], v_bw[i][j], nx, ny, bx, by,
						 m_max_disp, m_subpixel);
			}

		/* ---- for each pixel that is not too close to the boundary ---- */
		for (i = bx + m_window_size + m_max_disp; i < bx + nx - m_window_size - m_max_disp; i++)
			for (j = by + m_window_size + m_max_disp; j < by + ny - m_window_size - m_max_disp; j++)
			{
				/* determine corresponding location of (i,j) in f2 */
				ii = (int)(i + u_fw[i][j]);
				jj = (int)(j + v_fw[i][j]);

				/* check if this location (ii,jj) is in the image domain */
				if ((ii >= bx) && (ii < nx + bx) && (jj >= by) && (jj < ny + by))
				{
					float consistency = sqrt( pow( u_fw[i][j] + u_bw[ii][jj] ,2) + pow( v_fw[i][j] + v_bw[ii][jj] ,2) );
					if (consistency < m_occ_thres){
						u_fw[i][j] = 0;
					}
				}
			}
	}

	/* use forward displacements as solution */
	copy_matrix_2d(u_fw, u, nx, ny, bx, by);
	copy_matrix_2d(v_fw, v, nx, ny, bx, by);

	/* ---- free memory */
	FREE_MATRIX(1, 2 * m_max_disp + 1, 2 * m_max_disp + 1, cost_map);
	FREE_MATRIX(2, nx + 2 * bx, ny + 2 * by, u_fw, v_fw);
	FREE_MATRIX(2, nx + 2 * bx, ny + 2 * by, u_bw, v_bw);
}
/* ------------------------------------------------------------------------- */

#endif
