/*****************************************************************************/
/*                                                                           */
/*                   Copyright 08/2006 by Dr. Andres Bruhn,                  */
/*     Faculty of Mathematics and Computer Science, Saarland University,     */
/*                           Saarbruecken, Germany.                          */
/*                                                                           */
/*****************************************************************************/


#ifndef OF_HORN_SCHUNCK_INCLUDED
#define OF_HORN_SCHUNCK_INCLUDED

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "alloc_mem_linear.c"
#include "alloc_mem_linear_mult.c"
#include "funct_lib.c"


/* ------------------------------------------------------------------------- */

void horn_schunck_jacobi
(
                        /*****************************************************/
float **J_11,           /* in     : entry 11 of the motion tensor            */
float **J_22,           /* in     : entry 22 of the motion tensor            */
float **J_33,           /* in     : entry 33 of the motion tensor            */
float **J_12,           /* in     : entry 12 of the motion tensor            */
float **J_13,           /* in     : entry 13 of the motion tensor            */
float **J_23,           /* in     : entry 23 of the motion tensor            */
float **u,              /* in+out : x-component of displacement field        */
float **v,              /* in+out : y-component of displacement field        */
int   nx,               /* in     : size in x-direction                      */
int   ny,               /* in     : size in y-direction                      */
int   bx,               /* in     : boundary size in x-direction             */
int   by,               /* in     : boundary size in y-direction             */
float hx,               /* in     : grid spacing in x-direction              */
float hy,               /* in     : grid spacing in y-direction              */
float alpha             /* in     : smoothness weight                        */
                        /*****************************************************/
)

/*
 Computes one Jacobi iteration
*/

{
                        /*****************************************************/
int   i,j;              /* loop variables                                    */
float hx_2,hy_2;        /* time saver variables                              */
float xp,xm,yp,ym;      /* neighbourhood weights                             */
float sum;              /* central weight                                    */      
float **u_old;          /* x-component from old time step                    */
float **v_old;          /* y-component from old time step                    */
                        /*****************************************************/


/* define time saver variables */
hx_2=alpha/(hx*hx);
hy_2=alpha/(hy*hy);

/* alloc memory */
ALLOC_MATRIX(2, nx+2*bx, ny+2*by, 
	     &u_old,
	     &v_old);

/* copy previous result */
copy_matrix_2d(u,u_old,nx,ny,bx,by);
copy_matrix_2d(v,v_old,nx,ny,bx,by);

/* set boundaries zero */
set_bounds_2d(u_old,nx,ny,bx,by,0.0);
set_bounds_2d(v_old,nx,ny,bx,by,0.0);


/* Jacobi iteration */
for(i=bx;i<nx+bx;i++)
    for(j=by;j<ny+by;j++)
    {
	/* compute weights */
	
        /* TODO
          ----- fill in your code for computing the neighbourhood weights
                here - do not forget the correct boundary conditions  ----
	    */
        xm = hx_2 * (i>bx);
        xp = hx_2 * (i<nx+bx-1);
        ym = hy_2 * (j>by);
        yp = hy_2 * (j<ny+by-1);

		/* compute the sum of weights */
        sum = xp + xm + yp + ym;
		
		/* ------------------------------------------------------------ */

	/* perform iteration */

        /* TODO
         ----- fill in your code for the HS Jacobi iteration here ----
        */
       u[i][j] = ( -J_13[i][j] - ( J_12[i][j] *  v_old[i][j] - xm*u_old[i-1][j]
                                                             - xp*u_old[i+1][j]
                                                             - ym*u_old[i][j-1]
                                                             - yp*u_old[i][j+1] ) ) 
                    / ( J_11[i][j] + sum );

        v[i][j] = ( -J_23[i][j] - ( J_12[i][j] * u_old[i][j] - xm*v_old[i-1][j]
                                                             - xp*v_old[i+1][j]
                                                             - ym*v_old[i][j-1]
                                                             - yp*v_old[i][j+1] ) ) 
                    / ( J_22[i][j] + sum );
        		
		/* ------------------------------------------------------------ */

	}

/* free memory */
FREE_MATRIX(2, nx+2*bx, ny+2*by, 
	    u_old,
	    v_old);
}


/* ------------------------------------------------------------------------- */

void compute_motion_tensor
(
                        /*****************************************************/
float **f1,             /* in     : 1st image                                */
float **f2,             /* in     : 2nd image                                */
int   nx,               /* in     : size in x-direction                      */
int   ny,               /* in     : size in y-direction                      */
int   bx,               /* in     : boundary size in x-direction             */
int   by,               /* in     : boundary size in y-direction             */
float hx,               /* in     : grid spacing in x-direction              */
float hy,               /* in     : grid spacing in y-direction              */
float **J_11,           /* out    : entry 11 of the motion tensor            */
float **J_22,           /* out    : entry 22 of the motion tensor            */
float **J_33,           /* out    : entry 33 of the motion tensor            */
float **J_12,           /* out    : entry 12 of the motion tensor            */
float **J_13,           /* out    : entry 13 of the motion tensor            */
float **J_23            /* out    : entry 23 of the motion tensor            */
                        /*****************************************************/
)

/*
 Computes the motion tensor entries from the given image pair
*/

{
                        /*****************************************************/
int   i,j;              /* loop variables                                    */
float fx,fy,ft;         /* image derivatives                                 */
float hx_1,hy_1;        /* time saver variables                              */
                        /*****************************************************/


/* define time saver variables */
hx_1=1.0/(2.0*hx);
hy_1=1.0/(2.0*hy);

/* mirror boundaries */
mirror_bounds_2d(f1,nx,ny,bx,by);
mirror_bounds_2d(f2,nx,ny,bx,by);

/* compute motion tensor entries entries */        
for(i=bx;i<nx+bx;i++)
    for(j=by;j<ny+by;j++)
    {
	/* compute derivatives */

        /* TODO
            ----- fill in your code for computing the derivatives here ----
        */
       
       fx = 0.5 * ( (f2[i+1][j] - f2[i-1][j])*hx_1 + (f1[i+1][j] - f1[i-1][j])*hx_1 );
       fy = 0.5 * ( (f2[i][j+1] - f2[i][j-1])*hy_1 + (f1[i][j+1] - f1[i][j-1])*hy_1 );
       ft = f2[i][j] - f1[i][j];
       
		
		/* ------------------------------------------------------------ */

	/* set up motion tensor */

        /* TODO
         ----- fill in your code for the motion tensor entries here  ----
        */
		
       J_11[i][j] = fx*fx;
       J_22[i][j] = fy*fy;
       J_33[i][j] = ft*ft;
       J_12[i][j] = fx*fy;
       J_13[i][j] = fx*ft;
       J_23[i][j] = fy*ft;
       
		/* ------------------------------------------------------------ */
	
    }

}


/* ------------------------------------------------------------------------- */


void HORN_SCHUNCK
(
                        /*****************************************************/
float **f1,             /* in     : 1st image                                */
float **f2,             /* in     : 2nd image                                */
float **u,              /* out    : x-component of displacement field        */
float **v,              /* out    : y-component of displacement field        */
int   nx,               /* in     : size in x-direction                      */
int   ny,               /* in     : size in y-direction                      */
int   bx,               /* in     : boundary size in x-direction             */
int   by,               /* in     : boundary size in y-direction             */
float hx,               /* in     : grid spacing in x-direction              */
float hy,               /* in     : grid spacing in y-direction              */
float m_alpha,          /* in     : smoothness weight                        */
int   n_iter            /* in     : number of iterations                     */
                        /*****************************************************/
)

/* computes optic flow with Horn/Schunck */

{

                        /*****************************************************/
int   i,j;              /* loop variables                                    */
float **J_11;           /* entry 11 of the motion tensor                     */
float **J_22;           /* entry 22 of the motion tensor                     */
float **J_33;           /* entry 33 of the motion tensor                     */
float **J_12;           /* entry 12 of the motion tensor                     */
float **J_13;           /* entry 13 of the motion tensor                     */
float **J_23;           /* entry 23 of the motion tensor                     */
                        /*****************************************************/
  


/* ---- alloc memory ---- */
ALLOC_MATRIX (6, nx+2*bx,  ny+2*by, 
	      &J_11,
	      &J_22,
	      &J_33,
	      &J_12,
	      &J_13,
	      &J_23);


/* ---- initialise displacement field with zero ---- */
for (i=bx; i<bx+nx; i++)
    for (j=by; j<by+ny; j++)	 
    {
	u[i][j]=0;
	v[i][j]=0;
    }


/* ---- compute motion tensor ---- */
compute_motion_tensor(f1, f2, nx, ny, bx, by, hx, hy,
		      J_11, J_22, J_33, J_12, J_13, J_23);



/* ---- perform Jacobi iterations ---- */
for(i=1;i<=n_iter;i++)
 {
     horn_schunck_jacobi(J_11, J_22, J_33, J_12, J_13, J_23, 
			 u, v, nx, ny, bx, by, hx, hy, m_alpha);
 }



/* ---- free memory ---- */
FREE_MATRIX (6, nx+2*bx,  ny+2*by, 
	     J_11,
	     J_22,
	     J_33,
	     J_12,
	     J_13,
	     J_23);


}
/* ------------------------------------------------------------------------- */

#endif
