/*****************************************************************************/
/*                                                                           */
/*                   Copyright 08/2006 by Dr. Andres Bruhn,                  */
/*     Faculty of Mathematics and Computer Science, Saarland University,     */
/*                           Saarbruecken, Germany.                          */
/*                                                                           */
/*****************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <GL/glut.h>

/*---------------------------------------------------------------------------*/
/* include own libraries */

#include "alloc_mem_linear.c"
#include "alloc_mem_linear_mult.c"
#include "io_lib.c"
#include "bounds_lib.c"
#include "matrix_lib.c"
#include "mg_trans_lib.c"
#include "funct_lib.c"

/*---------------------------------------------------------------------------*/

                         /****************************************************/
char g_image1[80];       /* name of left image                               */
char g_image2[80];       /* name of right image                              */
char g_matrix_t[80];     /* name of ground truth fundamental matrix          */
char g_matrix_e[80];     /* name of estimated fundamental matrix             */
                         /*                                                  */
float **g_f1;            /* left image                                       */
float **g_f2;            /* right image                                      */
float ***g_f1_c;         /* colour copy of left image                        */
float ***g_f2_c;         /* colour copy of right image                       */
                         /*                                                  */
double **g_F_t;           /* ground truth fundamental matrix                  */
double **g_F_e;           /* estimated fundamental matrix                     */
                         /*                                                  */
float ***g_p6;           /* P6 image                                         */
                         /*                                                  */
float g_g_pixel_ratio;   /* pixel zoom factor                                */
                         /*                                                  */
int  g_nx, g_ny;         /* size of images in x- and y-direction             */
int  g_bx, g_by;         /* size of image boundary in x- and y-direction     */
int  g_px, g_py;         /* reference point in left image                    */
                         /*                                                  */
long g_position;         /* position in stream                               */ 
                         /*                                                  */
GLubyte *g_pixels;       /* pixel aray for Open-GL                           */
                         /****************************************************/


/*---------------------------------------------------------------------------*/

void drawGlutScene_from_image
(  
                         /****************************************************/
    float   ***p6,       /* in  : RGB image                                  */
    GLubyte *pixels,     /* use : display array                              */
    float   magnify,     /* in  : scaling factor                             */
    int     nx,          /* in  : size in x-direction                        */
    int     ny,          /* in  : size in y-direction                        */
    int     bx,          /* in  : boundary size in x-direction               */
    int     by           /* in  : boundary size in y-direction               */
                         /****************************************************/
)

/* visualises image with Open-GL */

{

    int odd;      /* flag for odd image size */
    int counter;  /* pixel index counter */
    int i,j;      /* loop variable */
 
   
    /* check if image size is odd */
    if (nx%4==0) odd=0;
    else odd=1;

    /* set pixel counter zero */
    counter=0;

    /* prepare Open-GL */    
    glViewport(0, 0, nx, ny);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, nx , 0, ny, -1.0, 1.0);
    glMatrixMode(GL_MODELVIEW);
    glDisable(GL_DITHER);
    glPixelZoom((GLfloat)magnify,(GLfloat)magnify);
      
    /* draw pixels in pixel array */
    for(i=by; i < ny+by;i++)
      {
	for(j=bx; j < nx+bx;j++)
	  {	   
	    pixels[counter++]=(GLubyte)(p6[j][ny+2*by-1-i][0]);
	    pixels[counter++]=(GLubyte)(p6[j][ny+2*by-1-i][1]);
	    pixels[counter++]=(GLubyte)(p6[j][ny+2*by-1-i][2]);
	  }
	if (odd==1) counter+=2;
      }

    /* draw pixels in draw buffer */
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0.0, 0.0, 0.0, 1.0);
    glRasterPos3f(0, 0, 0.0);
    glDrawPixels(nx,ny,GL_RGB,GL_UNSIGNED_BYTE,pixels);

    /* swap draw and display buffer */
    glutSwapBuffers();
    
    return;
}


/*---------------------------------------------------------------------------*/

void  draw_rectangle_3x3
(
                   /**********************************************************/
float  ***f,       /* in+out  : colour image                                 */
int    nx,         /* in      : size in x-direction                          */
int    ny,         /* in      : size in y-direction                          */
int    bx,         /* in      : boundary size in x-direction                 */
int    by,         /* in      : boundary size in y-direction                 */
float  px,         /* in      : x-coordinate of point without boundary       */
float  py,         /* in      : y-coordinate of point without boundary       */
float  *c          /* in      : RGB colour of line                           */
                   /**********************************************************/
)

/* draws 3x3 colour rectangle around a point */

{
                   /**********************************************************/
int   i;           /* loop variable                                          */
int   s,t;         /* pixel position with boundary                           */
                   /**********************************************************/


/* determine pixel position with boundary */
 
s=(int)round(px)+bx;
t=(int)round(py)+by;

 
/* draw rectangle around reference point */
for(i=0;i<3;i++)
 {
     f[s-1][t-1][i]=c[i];
     f[s-1][t  ][i]=c[i];
     f[s-1][t+1][i]=c[i];
     f[s  ][t-1][i]=c[i];
     f[s  ][t+1][i]=c[i];
     f[s+1][t-1][i]=c[i];
     f[s+1][t  ][i]=c[i];
     f[s+1][t+1][i]=c[i];
 }

 return;
}


/*---------------------------------------------------------------------------*/

void draw_line
(
                   /**********************************************************/
float  ***f,       /* in+out  : colour image                                 */
int    nx,         /* in      : size in x-direction                          */
int    ny,         /* in      : size in y-direction                          */
int    bx,         /* in      : boundary size in x-direction                 */
int    by,         /* in      : boundary size in y-direction                 */
float  *l,         /* in      : line normal                                  */
float  *c          /* in      : RGB colour of line                           */
                   /**********************************************************/
)

/* rasterize line with given line normal */

{

                   /**********************************************************/
int   i,j;         /* loop variables                                         */
float x,y;         /* pixel position without boundary                        */
int   s,t;         /* pixel position with boundary                           */
                   /**********************************************************/
  


/* ---- rasterize line from left ro right ---------------------------------- */

for(i=bx;i<nx+bx;i++)
 {
     x=i-bx;
     y=(l[1]*x+l[3])/(-l[2]);
     
     s=i;
     t=(int)round(y+by);
     
     if ((t>=by)&&(t<ny+by))
     {
	 f[s][t][0]=c[0];
	 f[s][t][1]=c[1];
	 f[s][t][2]=c[2];
     }
 }
 
/* ---- rasterize line from top to bottom ---------------------------------- */

for(j=by;j<ny+by;j++)
 {
     y=j-by;
     x=(l[2]*y+l[3])/(-l[1]);
     
     t=j;
     s=(int)round(x+bx);
     
     if ((s>=bx)&&(s<nx+bx))
     {
	 f[s][t][0]=c[0];
	 f[s][t][1]=c[1];
	 f[s][t][2]=c[2];
     }
 }

 return; 
}


/*---------------------------------------------------------------------------*/

void create_image_with_epipolar_lines
(
                   /**********************************************************/
float  **f1,       /* in  : left image                                       */
float  **f2,       /* in  : rigth image                                      */
float  ***f1_c,    /* in  : modified left image                              */
float  ***f2_c,    /* in  : modified rigth image                             */
float  ***p6,      /* out : RGB image                                        */
int    nx,         /* in  : size in x-direction                              */
int    ny,         /* in  : size in y-direction                              */
int    bx,         /* in  : boundary size in x-direction                     */
int    by,         /* in  : boundary size in y-direction                     */
int    px,         /* in  : x-coordinate of reference point                  */
int    py,         /* in  : y-coordinate of reference point                  */
double **F_t,      /* in  : ground truth fundamental matrix                  */
double **F_e       /* in  : estimated fundamental matrix                     */
                   /**********************************************************/
)

/* draw epipolar liner for ground truth and fundamental matrix */
   
{
                   /**********************************************************/
int   i,j;         /* loop variables                                         */
float l[4];        /* epipolar normal (projective coordinates)               */
float m[4];        /* pixel postion (projective coordinates)                 */
float c[3];        /* RGB colour                                             */
                   /**********************************************************/


/* copy both frames in colour images */
for(i=bx;i<nx+bx;i++)
    for(j=by;j<ny+by;j++)
    {
	f1_c[i][j][0]=f1[i][j];	
	f1_c[i][j][1]=f1[i][j];	
	f1_c[i][j][2]=f1[i][j];	
	
	f2_c[i][j][0]=f2[i][j];	
	f2_c[i][j][1]=f2[i][j];	
	f2_c[i][j][2]=f2[i][j];	
    }
    
/* if reference point is in the left frame */
if(px<nx)
 {       
     /* compute point in left frame (projective coorindates) */
     m[1]=px;
     m[2]=py;
     m[3]=1;
            
     /* draw a white rectangle around reference point */
     c[0]=255;
     c[1]=255;
     c[2]=255;     

     draw_rectangle_3x3(f1_c,nx,ny,bx,by,m[1],m[2],c);
     
     /* determine corresponding epipolar with F_t line in right frame */
     l[1]=F_t[1][1]*m[1]+F_t[1][2]*m[2]+F_t[1][3]*m[3];
     l[2]=F_t[2][1]*m[1]+F_t[2][2]*m[2]+F_t[2][3]*m[3];
     l[3]=F_t[3][1]*m[1]+F_t[3][2]*m[2]+F_t[3][3]*m[3];
     
     /* draw epipolar line for F_t in green */
     c[0]=0;
     c[1]=255;
     c[2]=0;  

     draw_line(f2_c,nx,ny,bx,by,l,c);
     
     /* determine corresponding epipolar with F_e line in right frame */
     l[1]=F_e[1][1]*m[1]+F_e[1][2]*m[2]+F_e[1][3]*m[3];
     l[2]=F_e[2][1]*m[1]+F_e[2][2]*m[2]+F_e[2][3]*m[3];
     l[3]=F_e[3][1]*m[1]+F_e[3][2]*m[2]+F_e[3][3]*m[3];
     
     /* draw epipolar line for F_e in red */
     c[0]=255;
     c[1]=0;
     c[2]=0;
  
     draw_line(f2_c,nx,ny,bx,by,l,c);
 }
 
 
/* if reference point is in the right frame */
if(px>=nx)
 {       
     /* compute point in right frame (projective coordinates */
     m[1]=px-nx;
     m[2]=py;
     m[3]=1;
     
     /* draw white rectangle around reference point */
     c[0]=255;
     c[1]=255;
     c[2]=255;  

     draw_rectangle_3x3(f2_c,nx,ny,bx,by,m[1],m[2],c);
     
     /* determine corresponding epipolar with F_t^T line in left frame */
     l[1]=F_t[1][1]*m[1]+F_t[2][1]*m[2]+F_t[3][1]*m[3];
     l[2]=F_t[1][2]*m[1]+F_t[2][2]*m[2]+F_t[3][2]*m[3];
     l[3]=F_t[1][3]*m[1]+F_t[2][3]*m[2]+F_t[3][3]*m[3];
     
     /* draw epipolar line for F_t in green */
     c[0]=0;
     c[1]=255;
     c[2]=0;  

     draw_line(f1_c,nx,ny,bx,by,l,c);
     
     /* determine corresponding epipolar with F_e line in right frame */
     l[1]=F_e[1][1]*m[1]+F_e[2][1]*m[2]+F_e[3][1]*m[3];
     l[2]=F_e[1][2]*m[1]+F_e[2][2]*m[2]+F_e[3][2]*m[3];
     l[3]=F_e[1][3]*m[1]+F_e[2][3]*m[2]+F_e[3][3]*m[3];
     
     /* draw epipolar line for F_e in red */
     c[0]=255;
     c[1]=0;
     c[2]=0;  

     draw_line(f1_c,nx,ny,bx,by,l,c);
 }
 

/* copy images in combined image (for output) */
for(i=bx;i<nx+bx;i++)
    for(j=by;j<ny+by;j++)
    {
	p6[i][j][0]=f1_c[i][j][0];	
	p6[i][j][1]=f1_c[i][j][1];	
	p6[i][j][2]=f1_c[i][j][2];	
    }
 
for(i=bx;i<nx+bx;i++)
    for(j=by;j<ny+by;j++)
    {
	p6[i+nx][j][0]=f2_c[i][j][0];
	p6[i+nx][j][1]=f2_c[i][j][1];
	p6[i+nx][j][2]=f2_c[i][j][2];	
    } 

 return;
}


/*---------------------------------------------------------------------------*/


void handleDoNothing()
{     
}

/*---------------------------------------------------------------------------*/


void handleDrawLine()
{    
    /* create image */
    create_image_with_epipolar_lines(g_f1,g_f2,g_f1_c,g_f2_c,g_p6,
				     g_nx,g_ny,g_bx,g_by,g_px,g_py,
				     g_F_t,g_F_e);    
}

/*---------------------------------------------------------------------------*/

void handleDraw()
{     
    /* draw flowfield */
    drawGlutScene_from_image(g_p6,g_pixels,g_g_pixel_ratio,
			     2*g_nx,g_ny,g_bx,g_by);
    /* refresh display */
    glutPostRedisplay();
}

/*---------------------------------------------------------------------------*/

void handleKeyboardspecial(int key, int x, int y)
{   
    switch(key) {	
	case GLUT_KEY_DOWN:
	    g_py++;
	    if (g_py>=g_ny) g_py=g_ny-1;
	    break;

	case GLUT_KEY_UP:
	    g_py--;
	    if (g_py<0) g_py=0;
	    break;

	case GLUT_KEY_RIGHT:
	    g_px++;
	    if (g_px>=2*g_nx) g_px=2*g_nx-1;
	    break;

	case GLUT_KEY_LEFT:
	    g_px--;
	    if (g_px<0) g_px=0;
	    break;
	         
	default:
	    printf("\nUnknown key pressed.");
	    break;
    }  
    handleDrawLine();  
    handleDraw();
}


/*---------------------------------------------------------------------------*/

void handleKeyboard(unsigned char key, int x, int y)
{
    switch(key) {	
	case 27:
	    exit(1);
	    break;
	default:
	    printf("\nUnknown key pressed.");
    }
    
    handleDrawLine();
    return;
    
}

/*---------------------------------------------------------------------------*/

void handleMouse(int button, int state, int cx, int cy)
{    
    if ((button==0)&&(state==1))
    {	
	g_px=cx;
	g_py=cy;
	handleDrawLine();
	handleDraw();
    }
}

/*-------------------------------------------------------------------------- */

int main (int argc, char* argv[])
{

/* --------- initialisations ----------------------------------------------- */

g_bx=1;
g_by=1;

g_px=1;
g_py=1;


/* ---------- read in console ---------------------------------------------- */

if (argc==6) 
 {   
     strcpy(g_image1,argv[1]);
     strcpy(g_image2,argv[2]);  
     strcpy(g_matrix_t,argv[3]);
     strcpy(g_matrix_e,argv[4]);                
     g_g_pixel_ratio=(int)atoi(argv[5]);     
 }
if (argc!=6)
 {
     printf("\n\n   %s <image1.pgm> <image2.pgm> <matrix_truth.fm> <matrix_estimation.fm> <pixel_ratio>\n\n\n",argv[0]);
     exit(0);    
 }

/* ---------- allocate memory ---------------------------------------------- */

ALLOC_MATRIX_DOUBLE(2, 3+2, 3+2, &g_F_t, &g_F_e);

	
/* ---------- read in information of first image --------------------------- */

read_pgm_header(g_image1,&g_position,&g_nx,&g_ny);


/* ---------- memory allocation -------------------------------------------- */

ALLOC_MATRIX(2, g_nx+2*g_bx, g_ny+2*g_by, &g_f1,  &g_f2);
ALLOC_CUBIX(2, g_nx+2*g_bx, g_ny+2*g_by, 3, &g_f1_c,  &g_f2_c);
ALLOC_CUBIX(1, 2*g_nx+2*g_bx, g_ny+2*g_by, 3, &g_p6);
g_pixels = (GLubyte *) malloc (2*g_nx*g_ny*3*sizeof(GLubyte));   


/* ---------- read input images -------------------------------------------- */

/* first frame */
read_pgm_header(g_image1,&g_position,&g_nx,&g_ny);
read_pgm_data(g_image1,g_position,g_f1,g_nx,g_ny,g_bx,g_by);
 
/* second frame */
read_pgm_header(g_image2,&g_position,&g_nx,&g_ny);
read_pgm_data(g_image2,g_position,g_f2,g_nx,g_ny,g_bx,g_by);


/* ---------- read fundamental matrices  ----------------------------------- */

/* ground truth fundamental matrix */
read_fundamental_matrix(g_matrix_t,g_F_t);

/* estiamted fundamental matrix */
read_fundamental_matrix(g_matrix_e,g_F_e);

/* ---------- M A I N   L O O P -------------------------------------------- */

// OpenGL-Fenster einblenden
glutInit(&argc, argv);
glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
glutInitWindowSize((int)round(2*g_nx*g_g_pixel_ratio),
		   (int)round(g_ny*g_g_pixel_ratio));
glutCreateWindow("EPIPOLAR LINE VIEWER");

// Handler registrieren
glutDisplayFunc(handleDraw);
glutIdleFunc(handleDoNothing);
glutKeyboardFunc(handleKeyboard);
glutSpecialFunc(handleKeyboardspecial);
glutMouseFunc(handleMouse);
// Main 
glutMainLoop();


/* ---------- Free memory -------------------------------------------------- */

FREE_MATRIX_DOUBLE(2, 3+2, 3+2, g_F_t, g_F_e);
FREE_CUBIX (1, 2*g_nx+2*g_bx, g_ny+2*g_by, 3, g_p6);
free(g_pixels);
FREE_MATRIX (2, g_nx+2*g_bx, g_ny+2*g_by, g_f1, g_f2);
FREE_CUBIX (2, g_nx+2*g_bx, g_ny+2*g_by, 3, g_f1_c, g_f2_c);

printf("\n\n\n");

return(0);
}








