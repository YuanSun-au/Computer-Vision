#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*--------------------------------------------------------------------------*/
/*                                                                          */
/*                       Hough Transform for Circles                        */
/*                                                                          */
/*  (Copyright Andres Bruhn and Sebastian Volz, 10/2012 , changed 10/2015)  */
/*                                                                          */
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/

void alloc_vector

	(float **vector, /* vector */
	 long n)		 /* size */

/* allocates storage for a vector of size n */

{
	*vector = (float *)malloc(n * sizeof(float));
	if (*vector == NULL)
	{
		printf("alloc_vector: not enough storage available\n");
		exit(1);
	}
	return;
}

/*--------------------------------------------------------------------------*/

void alloc_matrix

	(float ***matrix, /* matrix */
	 long nx,		  /* size in x direction */
	 long ny)		  /* size in y direction */

/* allocates storage for matrix of size nx * ny */

{
	long i;

	*matrix = (float **)malloc(nx * sizeof(float *));
	if (*matrix == NULL)
	{
		printf("alloc_matrix: not enough storage available\n");
		exit(1);
	}
	for (i = 0; i < nx; i++)
	{
		(*matrix)[i] = (float *)malloc(ny * sizeof(float));
		if ((*matrix)[i] == NULL)
		{
			printf("alloc_matrix: not enough storage available\n");
			exit(1);
		}
	}
	return;
}

/*--------------------------------------------------------------------------*/

void alloc_cubix

	(float ****cubix, /* cubix */
	 long nx,		  /* size in x direction */
	 long ny,		  /* size in y direction */
	 long nz)		  /* size in z direction */

/* allocates storage for cubix of size nx * ny * nz */

{
	long i, j;

	*cubix = (float ***)malloc(nx * sizeof(float **));
	if (*cubix == NULL)
	{
		printf("alloc_cubix: not enough storage available\n");
		exit(1);
	}
	for (i = 0; i < nx; i++)
	{
		(*cubix)[i] = (float **)malloc(ny * sizeof(float *));
		if ((*cubix)[i] == NULL)
		{
			printf("alloc_cubix: not enough storage available\n");
			exit(1);
		}

		for (j = 0; j < ny; j++)
		{
			(*cubix)[i][j] = (float *)malloc(nz * sizeof(float));
			if ((*cubix)[i][j] == NULL)
			{
				printf("alloc_cubix: not enough storage available\n");
				exit(1);
			}
		}
	}
	return;
}

/*--------------------------------------------------------------------------*/

void disalloc_vector

	(float *vector, /* vector */
	 long n)		/* size */

/* disallocates storage for a vector of size n */

{
	free(vector);
	return;
}

/*--------------------------------------------------------------------------*/

void disalloc_matrix

	(float **matrix, /* matrix */
	 long nx,		 /* size in x direction */
	 long ny)		 /* size in y direction */

/* disallocates storage for matrix of size nx * ny */

{
	long i;
	for (i = 0; i < nx; i++)
		free(matrix[i]);
	free(matrix);
	return;
}

/*--------------------------------------------------------------------------*/

void disalloc_cubix

	(float ***cubix, /* cubix */
	 long nx,		 /* size in x direction */
	 long ny,		 /* size in y direction */
	 long nz)		 /* size in z direction */

/* disallocates storage for cubix of size nx * ny * nz */

{
	long i, j;
	for (i = 0; i < nx; i++)
		for (j = 0; j < ny; j++)
			free(cubix[i][j]);
	for (i = 0; i < nx; i++)
		free(cubix[i]);
	free(cubix);
	return;
}

/*--------------------------------------------------------------------------*/

void alloc_vector_long

	(long **vector, /* vector */
	 long n)		/* size */

/* allocates storage for a vector of size n */

{
	*vector = (long *)malloc(n * sizeof(long));
	if (*vector == NULL)
	{
		printf("alloc_vector: not enough storage available\n");
		exit(1);
	}
	return;
}

/*--------------------------------------------------------------------------*/

void alloc_matrix_long

	(long ***matrix, /* matrix */
	 long nx,		 /* size in x direction */
	 long ny)		 /* size in y direction */

/* allocates storage for matrix of size nx * ny */

{
	long i;

	*matrix = (long **)malloc(nx * sizeof(long *));
	if (*matrix == NULL)
	{
		printf("alloc_matrix: not enough storage available\n");
		exit(1);
	}
	for (i = 0; i < nx; i++)
	{
		(*matrix)[i] = (long *)malloc(ny * sizeof(long));
		if ((*matrix)[i] == NULL)
		{
			printf("alloc_matrix: not enough storage available\n");
			exit(1);
		}
	}
	return;
}

/*--------------------------------------------------------------------------*/

void disalloc_vector_long

	(long *vector, /* vector */
	 long n)	   /* size */

/* disallocates storage for a vector of size n */

{
	free(vector);
	return;
}

/*--------------------------------------------------------------------------*/

void disalloc_matrix_long

	(long **matrix, /* matrix */
	 long nx,		/* size in x direction */
	 long ny)		/* size in y direction */

/* disallocates storage for matrix of size nx * ny */

{
	long i;
	for (i = 0; i < nx; i++)
		free(matrix[i]);
	free(matrix);
	return;
}

/*--------------------------------------------------------------------------*/

void alloc_circle_list

	(float ****c_list, /* list of points on circles of different radii */
	 long r_max,	   /* maximum radius */
	 long r_min)	   /* minimum radius */

/*
  allocate memory for circle list
*/

{
	long i, j, r; /* loop variables */
	long s1, s2;  /* auxiliary variable */
	float help;   /* auxiliary variable */
	long count;   /* number of points to be drawn for a certain circle */
	float **tmp;  /* aux image to draw circles and count points */

	/* allocate memory */
	alloc_matrix(&tmp, 2 * r_max + 1 + 2, 2 * r_max + 1 + 2);

	/* ---- determine number of points on circle with r_max ---- */
	r = r_max;

	/* set image to black */
	for (i = 1; i <= 2 * r + 1; i++)
		for (j = 1; j <= 2 * r + 1; j++)
		{
			tmp[i][j] = 0;
		}

	/* draw circle of radius r_max */

	/* left to right sweep */
	for (i = 1; i <= 2 * r + 1; i++)
	{
		/* solve circle equation for y^2 */
		help = r * r - (i - r - 1) * (i - r - 1);

		/* continue if no solution exists */
		if (help >= 0)
			help = round(sqrt(help));
		else
			continue;

		/* otherwise draw both solutions */
		s1 = r + 1 + help;
		tmp[i][s1] = 255.0;

		s2 = r + 1 - help;
		tmp[i][s2] = 255.0;
	}

	/* top to bottom sweep */
	for (j = 1; j <= 2 * r + 1; j++)
	{
		/* solve circle equation for x^2 */
		help = r * r - (j - r - 1) * (j - r - 1);

		/* continue if no solution exists */
		if (help >= 0)
			help = round(sqrt(help));
		else
			continue;

		/* otherwise draw both solutions */
		s1 = r + 1 + help;
		tmp[s1][j] = 255.0;

		s2 = r + 1 - help;
		tmp[s2][j] = 255.0;
	}

	/* count points on circle with r_max */
	count = 0;

	for (i = 1; i <= 2 * r + 1; i++)
		for (j = 1; j <= 2 * r + 1; j++)
		{
			if (tmp[i][j] == 255.0)
			{
				count++;
			}
		}

	/* free memory */
	disalloc_matrix(tmp, 2 * r_max + 1 + 2, 2 * r_max + 1 + 2);

	/* allocate memory for circle list */
	alloc_cubix(c_list, r_max + 1, count + 1, 2);

	return;
} /* alloc_circle_list */

/*--------------------------------------------------------------------------*/

void fill_circle_list

	(float ***c_list, /* list of points on circles of different radii */
	 long r_max,	  /* maximum radius */
	 long r_min)	  /* minimum radius */

/*
  compute relative coordinates of all pixels for each circle with a radius
  between r_min and r_max
*/

{
	long i, j, r; /* loop variables */
	long s1, s2;  /* auxiliary variable */
	float help;   /* auxiliary variable */
	long count;   /* number of points to be drawn for a certain circle */
	float **tmp;  /* aux image to draw circles and count points */

	/* allocate memory */
	alloc_matrix(&tmp, 2 * r_max + 1 + 2, 2 * r_max + 1 + 2);

	/* ---- compute points on circles for different r ---- */
	for (r = r_min; r <= r_max; r++)
	{
		/* set image to black */
		for (i = 1; i <= 2 * r + 1; i++)
			for (j = 1; j <= 2 * r + 1; j++)
			{
				tmp[i][j] = 0;
			}

		/* draw circle of radius r_max */

		/* left to right sweep */
		for (i = 1; i <= 2 * r + 1; i++)
		{
			/* solve circle equation for y^2 */
			help = r * r - (i - r - 1) * (i - r - 1);

			/* continue if no solution exists */
			if (help >= 0)
				help = round(sqrt(help));
			else
				continue;

			/* otherwise draw both solutions */
			s1 = r + 1 + help;
			tmp[i][s1] = 255.0;

			s2 = r + 1 - help;
			tmp[i][s2] = 255.0;
		}

		/* top to bottom sweep */
		for (j = 1; j <= 2 * r + 1; j++)
		{
			/* solve circle equation for x^2 */
			help = r * r - (j - r - 1) * (j - r - 1);

			/* continue if no solution exists */
			if (help >= 0)
				help = round(sqrt(help));
			else
				continue;

			/* otherwise draw both solutions */
			s1 = r + 1 + help;
			tmp[s1][j] = 255.0;

			s2 = r + 1 - help;
			tmp[s2][j] = 255.0;
		}

		count = 0;

		for (i = 1; i <= 2 * r + 1; i++)
			for (j = 1; j <= 2 * r + 1; j++)
			{

				if (tmp[i][j] == 255.0)
				{
					/* store relative coordinates of pixels */
					c_list[r][count + 1][0] = i - r - 1;
					c_list[r][count + 1][1] = j - r - 1;

					/* count points on circle with radius r */
					count++;
				}
			}

		/* store pixel number instead of zeroth pixel's coordinates */
		c_list[r][0][0] = count;
		c_list[r][0][1] = count;
	}

	/* free memory */
	disalloc_matrix(tmp, 2 * r_max + 1 + 2, 2 * r_max + 1 + 2);

	return;

} /* fill_circle_list */

/*--------------------------------------------------------------------------*/

void free_circle_list

	(float ***c_list, /* list of points on circles of different radii */
	 long r_max,	  /* maximum radius */
	 long r_min)	  /* minimum radius */

/*
  free memory for circle list
*/

{
	/* free memory for circle list */
	disalloc_cubix(c_list, r_max + 1, (c_list[r_max][0][0]) + 1, 2);

	return;

} /* free_circle_list */

/*--------------------------------------------------------------------------*/

void vote_circle

	(float **u,		  /* voting space,changed */
	 float ***c_list, /* list of points on circles of different radii */
	 long r_max,	  /* maximum radius */
	 long r_min,	  /* minimum radius */
	 long nx,		  /* pixel number in x-direction */
	 long ny,		  /* pixel number in y-direction */
	 long cx,		  /* x-coordinate of center of circle */
	 long cy,		  /* y-coordinate of center of circle */
	 long r)		  /* radius */

/*
 vote for circle
*/

{
	long k, i, j; /* loop variables */
	float help;   /* auxiliary variable */
	float k_max;  /* number of points */

	/* check if vote for desired radius is allowed */
	if ((r < r_min) || (r > r_max))
	{
		printf("\n radius not in the interval [r_min,r_max] !!!");
		exit(0);
	}

	/* read out number of points for circle of radius r */
	k_max = c_list[r][0][0];

	/* for each point on the circle of radius r */
	for (k = 1; k <= k_max; k++)
	{
		/* determine coordinate */
		i = cx + c_list[r][k][0];
		j = cy + c_list[r][k][1];

		/* vote if inside the image */
		if ((i >= 1) && (i <= nx) && (j >= 1) && (j <= ny))
			u[i][j] = u[i][j] + 1;
	}

	return;

} /* vote_circle */

/*--------------------------------------------------------------------------*/

void vote_hough

	(float ***h,	  /* Hough voting space */
	 float **u_mag,   /* edge image  */
	 float ***c_list, /* list of points on circles of different radii */
	 long r_max,	  /* maximum radius */
	 long r_min,	  /* minimum radius */
	 long nx,		  /* pixel number in x-direction */
	 long ny		  /* pixel number in y-direction */
	)

/*
 vote for each point that is an edge
*/

{
	long i, j, r; /* loop variables */

	for (i = 1; i <= nx; i++)
		for (j = 1; j <= ny; j++)
		{
			for (r = r_min; r <= r_max; r++)
			{
				vote_circle (u_mag, c_list, r_max, r_min, nx, ny, i, j, r) ;
			}
		}
	return;

} /* vote_hough */

/*--------------------------------------------------------------------------*/

void normalise_hough

	(float ***h,	  /* Hough voting space, changed */
	 float ***c_list, /* list of points on circles of different radii */
	 long r_max,	  /* maximum radius */
	 long r_min,	  /* minimum radius */
	 long nx,		  /* pixel number in x-direction */
	 long ny		  /* pixel number in y-direction */
	)

/* 
 normalises entries in Hough space by number of possible votes
*/

{
	long i, j, r; /* loop variables */

	for (r = r_min; r <= r_max; r++)
		for (i = 1; i <= nx; i++)
			for (j = 1; j <= ny; j++)
			{
				/* normalise hough space votes by numer of possible votes */
				h[r][i][j] = h[r][i][j] / c_list[r][0][0];
			}

	return;

} /* normalise_hough */

/*--------------------------------------------------------------------------*/

void threshold_hough

	(float ***h,   /* Hough voting space, changed */
	 long r_max,   /* maximum radius */
	 long r_min,   /* minimum radius */
	 long nx,	  /* pixel number in x-direction */
	 long ny,	  /* pixel number in y-direction */
	 float t_hough /* threshold for hough transform */
	)

/*
 computes thresholded Hough space
*/

{
	long i, j, r; /* loop variables */

	for (r = r_min; r <= r_max; r++)
		for (i = 1; i <= nx; i++)
			for (j = 1; j <= ny; j++)
			{
				/* if more than t_hough percent of votes */
				if (h[r][i][j] < t_hough)
					h[r][i][j] = 0.0;
			}

	return;

} /* threshold_hough */

/*--------------------------------------------------------------------------*/

void denormalise_hough

	(float ***h,	  /* Hough voting space, changed */
	 float ***c_list, /* list of points on circles of different radii */
	 long r_max,	  /* maximum radius */
	 long r_min,	  /* minimum radius */
	 long nx,		  /* pixel number in x-direction */
	 long ny		  /* pixel number in y-direction */
	)

/* 
 denormalises entries in Hough space by number of possible votes
*/

{
	long i, j, r; /* loop variables */

	for (r = r_min; r <= r_max; r++)
		for (i = 1; i <= nx; i++)
			for (j = 1; j <= ny; j++)
			{
				/* denormalise hough space votes by numer of possible votes */
				h[r][i][j] = h[r][i][j] * c_list[r][0][0];
			}

	return;

} /* denormalise_hough */

/*--------------------------------------------------------------------------*/

void nonmaxima_suppression_hough

	(float ***h,  /* normalised Hough voting space, changed */
	 long r_max,  /* maximum radius */
	 long r_min,  /* minimum radius */
	 long nx,	 /* pixel number in x-direction */
	 long ny,	 /* pixel number in y-direction */
	 long r_dist, /* minimum distance between circle radii */
	 long i_dist, /* minimum distance between circle centers */
	 long j_dist  /* minimum distance between circle centers */
	)

/*
 perform nonmaxima suppression in Hough space
*/

{
	long i, j, r;			  /* loop variables */
	long k, l, s;			  /* loop variables */
	long k_min, k_max;		  /* bounding box neighbourhood window */
	long l_min, l_max;		  /* bounding box neighbourhood window */
	long s_min, s_max;		  /* bounding box neighbourhood window */
	float max;				  /* local maximum */
	long max_k, max_l, max_s; /* coordinates of maximum */

	for (r = r_min; r <= r_max; r++)
		for (i = 1; i <= nx; i++)
			for (j = 1; j <= ny; j++)
			{
				if (h[r][i][j] > 0)
				{
					// determine neighbourhood window around point (r,i,j)
					// crop window at image boundaries
					s_min = r - r_dist;
					if (s_min < r_min)
						s_min = r_min;
					s_max = r + r_dist;
					if (s_max > r_max)
						s_max = r_max;
					k_min = i - i_dist;
					if (k_min < 1)
						k_min = 1;
					k_max = i + i_dist;
					if (k_max > nx)
						k_max = nx;
					l_min = j - i_dist;
					if (l_min < 1)
						l_min = 1;
					l_max = j + i_dist;
					if (l_max > ny)
						l_max = ny;

					// determine maximum in window and delete entries */
					max = h[r][i][j];
					max_s = r;
					max_k = i;
					max_l = j;

					for (s = s_min; s <= s_max; s++)
						for (k = k_min; k <= k_max; k++)
							for (l = l_min; l <= l_max; l++)
							{
								if (h[s][k][l] > max)
								{
									max = h[s][k][l];
									max_s = s;
									max_k = k;
									max_l = l;
								}
								h[s][k][l] = 0.0;
							}

					// restore maximum entry
					h[max_s][max_k][max_l] = max;
				}
			}

	return;

} /* nonmaxima_suppression_hough */

/*--------------------------------------------------------------------------*/

void count_hough

	(float ***h,  /* Hough voting space, changed */
	 long r_max,  /* maximum radius */
	 long r_min,  /* minimum radius */
	 long nx,	 /* pixel number in x-direction */
	 long ny,	 /* pixel number in y-direction */
	 long *number /* return number of objects */
	)

/*
 count non-zero entries in Hough space
*/

{
	long i, j, r; /* loop variables */

	*number = 0;
	for (r = r_min; r <= r_max; r++)
		for (i = 1; i <= nx; i++)
			for (j = 1; j <= ny; j++)
			{
				/* if more than t_hough percent of votes */
				if (h[r][i][j] > 0.0)
					(*number)++;
			}

	printf("\n ... found %ld independent components! \n\n", *number);

	return;

} /* count_hough */

/*--------------------------------------------------------------------------*/

void draw_circle

	(float **u,		  /* image, changed */
	 float ***c_list, /* list of points on circles of different radii */
	 long r_max,	  /* maximum radius */
	 long r_min,	  /* minimum radius */
	 long nx,		  /* pixel number in x-direction */
	 long ny,		  /* pixel number in y-direction */
	 long cx,		  /* x-coordinate of center of circle */
	 long cy,		  /* y-coordinate of center of circle */
	 long r)		  /* radius */

/*
 draw circle
*/

{
	long k, i, j; /* loop variables */
	float help;   /* auxiliary variable */
	float k_max;  /* number of points */

	/* check if vote for desired radius is allowed */
	if ((r < r_min) || (r > r_max))
	{
		printf("\n radius not in the interval [r_min,r_max] !!!");
		exit(0);
	}

	/* for each point on the circle of radius r */
	k_max = c_list[r][0][0];

	for (k = 1; k <= k_max; k++)
	{
		/* determine coordinate */
		i = cx + c_list[r][k][0];
		j = cy + c_list[r][k][1];

		/* vote if inside the image */
		if ((i >= 1) && (i <= nx) && (j >= 1) && (j <= ny))
			u[i][j] = 0.0;
	}

	return;

} /* draw_circle */

/*--------------------------------------------------------------------------*/

void draw_circle_image

	(
		float **u,		 /* original image, changed */
		float ***h,		 /* Hough voting space */
		float ***c_list, /* list of points on circles of different radii */
		long r_max,		 /* maximum radius */
		long r_min,		 /* minimum radius */
		long nx,		 /* pixel number in x-direction */
		long ny			 /* pixel number in y-direction */
	)

/*
 draws circles in image based on Hough space
*/

{
	long i, j, r; /* loop variables */

	for (r = r_min; r <= r_max; r++)
		for (i = 1; i <= nx; i++)
			for (j = 1; j <= ny; j++)
			{
				if (h[r][i][j] > 0.0)
				{
					/* draw circle with width of 3 pixels */
					draw_circle(u, c_list, r_max, r_min, nx, ny, i, j, r);
					if (r + 1 <= r_max)
						draw_circle(u, c_list, r_max, r_min, nx, ny, i, j, r + 1);
					if (r + 2 <= r_max)
						draw_circle(u, c_list, r_max, r_min, nx, ny, i, j, r + 2);
				}
			}

	return;

} /* draw_circle_image */

/*--------------------------------------------------------------------------*/

void draw_full_circle_interior

	(float **u,		 /* image, changed */
	 long nx,		 /* pixel number in x-direction */
	 long ny,		 /* pixel number in y-direction */
	 long cx,		 /* x-coordinate of center of circle */
	 long cy,		 /* y-coordinate of center of circle */
	 long r,		 /* radius */
	 long greyvalue) /* greyvalue */

/*
 draw circle
*/

{
	long i, j;			 /* loop variables */
	long xs, ys, xe, ye; /* corners of enclosing rectangle */
	float help;			 /* auxiliary variable */
	float k_max;		 /* number of points */

	/* compute reduced radius for interior */
	r = r;
	if (r <= 1)
		r = 1;

	/* determine corners of enclosing rectangle */
	xs = -r;
	ys = -r;
	xe = r;
	ye = r;

	if ((cx - r - 1) < 0)
		xs = -r - (cx - r - 1);
	if ((cy - r - 1) < 0)
		ys = -r - (cy - r - 1);
	if ((cx + r + 1) > nx + 1)
		xe = r - (cx + r - nx);
	if ((cy + r + 1) > ny + 1)
		ye = r - (cy + r - ny);

	/* draw circle */
	for (i = xs; i <= xe; i++)
		for (j = ys; j <= ye; j++)
		{
			if ((i * i + j * j) <= (r * r))
				u[cx + i][cy + j] = greyvalue;
		}

	return;

} /* draw_full_circle */

/*--------------------------------------------------------------------------*/

void draw_full_circle_image

	(
		float **u,		 /* original image, changed */
		float ***h,		 /* Hough voting space */
		float ***c_list, /* list of points on circles of different radii */
		long r_max,		 /* maximum radius */
		long r_min,		 /* minimum radius */
		long nx,		 /* pixel number in x-direction */
		long ny,		 /* pixel number in y-direction */
		long number		 /* numer of objects */
	)

/*
 draws circles in image based on Hough space
*/

{
	long i, j, r;   /* loop variables */
	long greyvalue; /* greyvalue */
	long count;		/* counting variale */

	/* initialise counter */
	count = 0;

	/* draw circles with corresponding grey value */
	for (r = r_min; r <= r_max; r++)
		for (i = 1; i <= nx; i++)
			for (j = 1; j <= ny; j++)
			{
				if (h[r][i][j] > 0.0)
				{
					/* use different grey value for each circle */
					greyvalue = 255.0 / (number + 20) * (count + 10);
					/* draw full circles */
					draw_full_circle_interior(u, nx, ny, i, j, r, greyvalue);
					count++;
				}
			}

	return;

} /* draw_full_circle_image */

/*--------------------------------------------------------------------------*/

int main()

{
	char row[80];			  /* for reading data */
	char in1[80];			  /* for reading data */
	char in2[80];			  /* for reading data */
	char out1[80];			  /* for reading data */
	char out2[80];			  /* for reading data */
	char tmp1[80];			  /* for copying information */
	char tmp2[80];			  /* for copying information */
	float **f;				  /* original image */
	float **u_mag;			  /* gradient magnitude */
	float **b;				  /* binary image for counting coins */
	float **g;				  /* output image */
	float ***h;				  /* Hough voting space */
	long i, j, k;			  /* loop variables */
	long nx, ny;			  /* image size in x, y direction */
	float hx, hy;			  /* grid size in x, y direction */
	FILE *inimage, *outimage; /* input file, output file */
	float t_hough;			  /* Hough threshold */
	long r_max;				  /* maximum radius in Hough space */
	long r_min;				  /* minimum radius in Hough space */
	long r_dist;			  /* minimum distance between radii */
	long s_dist;			  /* minimum distance betwen centers */
	unsigned char byte;		  /* for data conversion */
	long number;			  /* number of objects */
	float ***c_list;		  /* list of points on circles of different radii */

	printf("\n");
	printf("HOUGH TRANSFORM FOR CIRCLES\n\n");
	printf("****************************************************************\n");
	printf("\n");
	printf("    Copyright 2015 by Andres Bruhn and Sebastian Volz   \n");
	printf("    Institute for Visualization and Interactive Systems \n");
	printf("    University of Stuttgart, Germany\n");
	printf("\n");
	printf("    All rights reserved. Unauthorized usage,\n");
	printf("    copying, hiring, and selling prohibited.\n");
	printf("\n");
	printf("    Send bug reports to\n");
	printf("    bruhn@vis.uni-stuttgart.de\n");
	printf("\n");
	printf("****************************************************************\n\n");

	/* -------------- R E A D   I N ------------- */

	/* ---- read input image (pgm format P5) ---- */

	/* read image name */
	printf("original image:                                         ");
	gets(in1);
	printf("edge map:                                               ");
	gets(in2);

	/* open pgm file and read header */
	inimage = fopen(in1, "rb");
	fgets(row, 80, inimage);
	fgets(row, 80, inimage);
	while (row[0] == '#')
		fgets(row, 80, inimage);
	sscanf(row, "%ld %ld", &nx, &ny);
	fgets(row, 80, inimage);

	/* allocate storage */
	alloc_matrix(&f, nx + 2, ny + 2);
	alloc_matrix(&u_mag, nx + 2, ny + 2);
	alloc_matrix(&b, nx + 2, ny + 2);
	alloc_matrix(&g, nx + 2, ny + 2);

	/* read image data */
	for (j = 1; j <= ny; j++)
		for (i = 1; i <= nx; i++)
			f[i][j] = (float)getc(inimage);
	fclose(inimage);

	/* read edge image */
	inimage = fopen(in2, "rb");
	fgets(row, 80, inimage);
	fgets(row, 80, inimage);
	while (row[0] == '#')
	{
		strcpy(tmp1, tmp2);
		strcpy(tmp2, row);
		fgets(row, 80, inimage);
	}
	sscanf(row, "%ld %ld", &nx, &ny);
	fgets(row, 80, inimage);

	/* read image data */
	for (j = 1; j <= ny; j++)
		for (i = 1; i <= nx; i++)
			u_mag[i][j] = (float)getc(inimage);
	fclose(inimage);

	/* ---- read other parameters ---- */

	printf("maximum circle radius:                                 ");
	gets(row);
	sscanf(row, "%ld", &r_max);
	printf("minimum circle radius:                                 ");
	gets(row);
	sscanf(row, "%ld", &r_min);
	printf("Hough threshold t_hough [0..1]:                        ");
	gets(row);
	sscanf(row, "%f", &t_hough);
	printf("output image (circles):                                ");
	gets(out1);
	printf("output image (rescaled class map):                     ");
	gets(out2);
	printf("\n");

	/* -------------- C O M P U T A T I O N S ------------- */

	/* allocate memory for Hough space */
	alloc_cubix(&h, r_max + 1, nx + 2, ny + 2);

	/* ---- define grid size ---- */
	hx = 1.0;
	hy = 1.0;

	/* allocate memory for circle list */
	alloc_circle_list(&c_list, r_max, r_min);

	/* ---- compute for which points a circle of radius r votes ---- */
	fill_circle_list(c_list, r_max, r_min);

	/* ---- initialise Hough space (voting space) ---- */
	for (k = 0; k <= r_max; k++)
		for (j = 1; j <= ny; j++)
			for (i = 1; i <= nx; i++)
				h[k][i][j] = 0;

	/* ---- vote for each point that is an edge ---- */
	vote_hough(h, u_mag, c_list, r_max, r_min, nx, ny);

	/* normalise hough space */
	normalise_hough(h, c_list, r_max, r_min, nx, ny);

	/* threshold hough space */
	threshold_hough(h, r_max, r_min, nx, ny, t_hough);

	/* denormalise hough space */
	denormalise_hough(h, c_list, r_max, r_min, nx, ny);

	/* perform non-maxima suppression in Hough space */
	r_dist = r_min;
	s_dist = r_min;
	nonmaxima_suppression_hough(h, r_max, r_min, nx, ny, r_dist, s_dist, s_dist);

	/* count number of circles in Hough space */
	count_hough(h, r_max, r_min, nx, ny, &number);

	/* copy original image */
	for (j = 1; j <= ny; j++)
		for (i = 1; i <= nx; i++)
			g[i][j] = f[i][j];

	/* mark all found circles in copied image */
	draw_circle_image(g, h, c_list, r_max, r_min, nx, ny);

	/* initialise binary image for counting coins */
	for (j = 1; j <= ny; j++)
		for (i = 1; i <= nx; i++)
			b[i][j] = 255.0;

	/* draw all found circles in binary image */
	draw_full_circle_image(b, h, c_list, r_max, r_min, nx, ny, number);

	/* free memory for circle list */
	free_circle_list(c_list, r_max, r_min);

	/* -------------- W R I T E   O U T ------------- */

	/* ---- write output image (pgm format P5) ---- */

	/* open file and write header (incl. filter parameters) */
	outimage = fopen(out1, "wb");
	fprintf(outimage, "P5 \n");
	fprintf(outimage, "# Hough Transform for Circles \n");
	fprintf(outimage, "# initial image:         %s\n", in1);
	fprintf(outimage, "%s", tmp1);
	fprintf(outimage, "%s", tmp2);
	fprintf(outimage, "# maximum radius:        %ld\n", r_max);
	fprintf(outimage, "# minimum radius:        %ld\n", r_min);
	fprintf(outimage, "# Hough threshold:       %f\n", t_hough);
	fprintf(outimage, "%ld %ld \n255\n", nx, ny);

	/* write image data and close file */
	for (j = 1; j <= ny; j++)
		for (i = 1; i <= nx; i++)
		{
			if (g[i][j] < 0.0)
				byte = (unsigned char)(0.0);
			else if (g[i][j] > 255.0)
				byte = (unsigned char)(255.0);
			else
				byte = (unsigned char)(g[i][j]);
			fwrite(&byte, sizeof(unsigned char), 1, outimage);
		}
	fclose(outimage);
	printf("output image %s successfully written\n\n", out1);

	/* ---- write output image (pgm format P5) ---- */

	/* open file and write header (incl. filter parameters) */
	outimage = fopen(out2, "wb");
	fprintf(outimage, "P5 \n");
	fprintf(outimage, "# Rescaled Class Map \n");
	fprintf(outimage, "# initial image:         %s\n", in1);
	fprintf(outimage, "%s", tmp1);
	fprintf(outimage, "%s", tmp2);
	fprintf(outimage, "# maximum radius:        %ld\n", r_max);
	fprintf(outimage, "# minimum radius:        %ld\n", r_min);
	fprintf(outimage, "# Hough threshold:       %f\n", t_hough);
	fprintf(outimage, "%ld %ld \n255\n", nx, ny);

	/* write image data and close file */
	for (j = 1; j <= ny; j++)
		for (i = 1; i <= nx; i++)
		{
			if (b[i][j] < 0.0)
				byte = (unsigned char)(0.0);
			else if (b[i][j] > 255.0)
				byte = (unsigned char)(255.0);
			else
				byte = (unsigned char)(b[i][j]);
			fwrite(&byte, sizeof(unsigned char), 1, outimage);
		}
	fclose(outimage);
	printf("output image %s successfully written\n\n", out2);

	/* ---- disallocate storage ---- */
	disalloc_matrix(f, nx + 2, ny + 2);
	disalloc_matrix(u_mag, nx + 2, ny + 2);
	disalloc_matrix(b, nx + 2, ny + 2);
	disalloc_matrix(g, nx + 2, ny + 2);
	disalloc_cubix(h, r_max + 1, nx + 2, ny + 2);

	return (0);
}
