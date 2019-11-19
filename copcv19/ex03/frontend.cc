/*****************************************************************************/
/*                                                                           */
/*               Copyright 08/2006 by Dr. Andres Bruhn and,                  */
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
#include <sys/time.h>
#include <GL/glut.h>

/*---------------------------------------------------------------------------*/
/* include own libraries */

#include "alloc_mem_linear_mult.c"
#include "io_lib.c"
#include "bounds_lib.c"
#include "matrix_lib.c"
#include "conv_lib.c"
#include "funct_lib.c"
#include "of_lib.c"
#include "horn_schunck.c"

/*---------------------------------------------------------------------------*/

                         /****************************************************/
char   g_image1[80];     /* name of 1st image                                */
char   g_image2[80];     /* name of 2nd image                                */
char   g_ref[80];        /* name of problem ground truth file                */
                         /*                                                  */
float  **g_f1;           /* 1st image                                        */
float  **g_f2;           /* 2nd image                                        */
float  **g_f1_s;         /* smoothed 1st image                               */
float  **g_f2_s;         /* smoothed 2nd image                               */
float  **g_u;            /* displacement component in x-direction            */
float  **g_v;            /* displacement component in y-direction            */
float  **g_u_ref;        /* problem ground truth in x-direction              */
float  **g_v_ref;        /* problem ground truth in y-direction              */
                         /*                                                  */
float  ***g_p6;          /* colour array to draw image                       */
                         /*                                                  */
int    g_nx,g_ny;        /* size of both frames in x- and y-drection         */
int    g_bx,g_by;        /* size of image boundary                           */
float  g_hx,g_hy;        /* pixel size of grid                               */
                         /*                                                  */
                         /*                                                  */
float  g_e_sigma;        /* standard deviation of the Gaussian presmoothing  */
                         /*                                                  */
float  g_m_alpha;        /* weight of the smoothness term                    */
int    g_n_iter;         /* number of iterations                             */
float  g_g_pixel_ratio;  /* pixel zoom factor                                */
float  g_g_max_disp;     /* maxium displacement magnitude                    */
                         /*                                                  */
                         /*                                                  */
int    g_active_param;   /* active parameter                                 */
int    g_direct_compute; /* flag for direct computation                      */
                         /*                                                  */
                         /*                                                  */
float  g_aae;            /* average angualer error                           */
float  g_al2e;           /* average l2 norm error                            */
float  g_ref_density;    /* density of problem ground truth                  */
float  g_est_density;    /* density of estimation                            */
                         /*                                                  */
                         /*                                                  */
long   g_position;       /* position in stream for reading input             */
                         /*                                                  */
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



void showParams()
{
 print_console_header("Model");

 printf("\n Horn/Schunck Approach");
 
 print_console_header("Model Parameters");

 if (g_active_param==1)
     printf("\n (a) (Alpha)                     %4.6lf",(double)g_m_alpha);
 else
     printf("\n (a)  Alpha                      %4.6lf",(double)g_m_alpha);

                          
 print_console_header("Preprocessing");

 if (g_active_param==0)
     printf("\n (p) (Sigma)                     %4.3lf",(double)g_e_sigma);
 else
     printf("\n (p)  Sigma                      %4.3lf",(double)g_e_sigma);

 print_console_header("Numerical Parameters");

 if (g_active_param==2)
     printf("\n (i) (Iterations)                %d",g_n_iter);
 else
     printf("\n (i)  Iterations                 %d",g_n_iter);

 print_console_header("Visualisation Parameters");

 if (g_active_param==4)
     printf("\n (g) (Maximum Shown Disp)        %4.3lf",(double)g_g_max_disp);
 else
     printf("\n (g)  Maximum Shown Disp         %4.3lf",(double)g_g_max_disp);

 if (g_direct_compute==1)
     printf("\n (,)  Direct Computation         ON");
 else
     printf("\n (,)  Direct Computation         OFF");

 print_console_header("Quality Measures");

 if(g_u_ref!=NULL)
     printf("\n      Average Angular Error      %4.10lf",(double)g_aae);

 print_console_line();
 fflush(stdout);
 return;
}


/*---------------------------------------------------------------------------*/


void handleComputeDisplacements()
{           
            		
  /* presmooth both frames */
  presmooth_2d(g_f1,g_f1_s,g_nx,g_ny,g_bx,g_by,g_hx,g_hy,g_e_sigma,g_e_sigma);
  presmooth_2d(g_f2,g_f2_s,g_nx,g_ny,g_bx,g_by,g_hx,g_hy,g_e_sigma,g_e_sigma);
  
  
  /* initialise displacement vector field with zero */    
  set_matrix_2d(g_u,g_nx+2*g_bx,g_ny+2*g_by,0,0,(float)0.0);
  set_matrix_2d(g_v,g_nx+2*g_bx,g_ny+2*g_by,0,0,(float)0.0);
  
  
  /* set boundaries zero */
  set_bounds_2d(g_u,g_nx,g_ny,g_bx,g_by,(float)0.0);
  set_bounds_2d(g_v,g_nx,g_ny,g_bx,g_by,(float)0.0);
  set_bounds_2d(g_f1_s,g_nx,g_ny,g_bx,g_by,(float)0.0);
  set_bounds_2d(g_f2_s,g_nx,g_ny,g_bx,g_by,(float)0.0);     
  
  
  /* compute displacement vector field */    
  HORN_SCHUNCK(g_f1_s,g_f2_s,g_u,g_v,g_nx,g_ny,g_bx,g_by,
	       g_hx, g_hy, g_m_alpha, g_n_iter);
  
  
  /* if ground truth available compute average angular error */
  if(g_u_ref!=NULL)
    calculate_errors_2d(g_u_ref,g_v_ref,g_u,g_v,g_nx,g_ny,g_bx,g_by,
			&g_aae,&g_al2e,&g_ref_density,&g_est_density);
  
  return;
}


/*---------------------------------------------------------------------------*/

void handleDraw()
{
   
    /* create image */
    convert_displacements_to_image(g_u,g_v,g_f1_s,g_nx,g_ny,g_bx,g_by,
				   (float)0.0,g_g_max_disp,g_p6);      
    
    /* draw displacement field field in pixel array */
    drawGlutScene_from_image(g_p6,g_pixels,g_g_pixel_ratio,
			     2*g_nx,g_ny,g_bx,g_by);

    /* draw pixel array on screen */
    glutPostRedisplay();
}


/*---------------------------------------------------------------------------*/

void handleKeyboardspecial(int key, int x, int y)
{   
  /* keyboard handler */

  switch(key) {
  case GLUT_KEY_F1:
    /* call help menu */
    print_console_line();
    printf("\n\n F1     ........this help\n");
    printf(" F2     ........write out displacement field in colour code \n");
    printf(" ESC    ........program termination\n");
    printf(" p      ........select presmoothing parameter\n");	    
    printf(" a      ........select smoothness weight\n");
    printf(" i      ........select number of iterations\n");
    printf(" g      ........select max displacement magnitude\n"); 	      
    printf(" up     ........increase active parameter\n");
    printf(" down   ........decrease active parameter\n");
    printf(" ,      ........direct computation on/off \n");
    printf(" .      ........compute displacement field \n");
    print_console_line();
    fflush(stdout);
    break;
    
  case GLUT_KEY_F2:	
    /* write out displacement field */
    if (g_u_ref!=NULL)
      {
	int i,j;
	for (i=g_bx; i<g_nx+g_bx; i++)	
	  for (j=g_by; j<g_ny+g_by; j++)		
	    {
	      if((g_u_ref[i][j]==100.0)&&(g_v_ref[i][j]==100.0))
		{
		  g_u[i][j]=100.0;
		  g_v[i][j]=100.0;
		}
	    }
      }    
    write_ppm_blank_header("out.ppm",g_nx,g_ny);
    write_ppm_data("out.ppm",g_p6,g_nx,g_ny,g_bx+g_nx,g_by);
    break; 	   
 	 	   				    
  case GLUT_KEY_DOWN:    
    /* decrease sigma */
    if (g_active_param==0) 
      {
	g_e_sigma=g_e_sigma/1.1;
	if (g_e_sigma<0.3) 
	  g_e_sigma=0.0;
	if (g_direct_compute==1) handleComputeDisplacements();break;
      }	    
    /* decrease alpha */
    if (g_active_param==1)
      {
	g_m_alpha/=sqrt(sqrt(sqrt(sqrt(10))));
	if (g_m_alpha<0.01) g_m_alpha=0.01;
	if (g_direct_compute==1) handleComputeDisplacements();break;
      }     
    /* decrease number of iterations */
    if (g_active_param==2) 
      {
	  g_n_iter-=100;
	if (g_n_iter<0) 
	  g_n_iter=0;
	if (g_direct_compute==1) handleComputeDisplacements();break;
      }	       
    /* decrease displacement field magnitude (visualsiation) */
    if (g_active_param==4) 
      {
	if (g_g_max_disp>=0.01) 
	  g_g_max_disp=g_g_max_disp/1.1;
	break;
      }

  case GLUT_KEY_UP:    
    /* increase sigma */
    if (g_active_param==0)
      {
	if (g_e_sigma<0.30) 
	  g_e_sigma=0.3;
	else 
	  g_e_sigma=g_e_sigma*1.1;
	if (g_direct_compute==1) handleComputeDisplacements();break;
      }
    /* increase alpha */
    if (g_active_param==1)
      {
	g_m_alpha*=sqrt(sqrt(sqrt(sqrt(10))));     
	if (g_direct_compute==1) handleComputeDisplacements();break;
      }    
    /* increase number of iterations */
    if (g_active_param==2) 
      {
	  g_n_iter+=100;
	if (g_direct_compute==1) handleComputeDisplacements();break;
      }	  
    /* increase displacement field magnitude (visualisation) */
    if (g_active_param==4)
      {
	g_g_max_disp=g_g_max_disp*1.1;
	break;
      }	
  default:
    printf("\nUnknown key pressed (Code %d).",key);
   
  }    
  
  /* show new parameters and displacement field */
  showParams();    
  handleDraw();
  return;
}


/*---------------------------------------------------------------------------*/

void handleKeyboard(unsigned char key, int x, int y)
{
  /* keyboard handler */
  switch(key) {
  case 'p':
    g_active_param=0;break;
  case 'a':
    g_active_param=1;break; 
  case 'i':
    g_active_param=2;break; 
  case 'g':
    g_active_param=4;break;
  case 44:  //,
    if (g_direct_compute==1) {g_direct_compute=0;break;}
    if (g_direct_compute==0) {g_direct_compute=1;break;}
  case 46:  //.
    handleComputeDisplacements();break;
  case 27:
    exit(1);	
  default:
    printf("\nUnknown key pressed (Code %d).",key);
  }    

  /* show new parameters and displacement field */
  showParams();
  handleDraw();
  return;       
}

/*---------------------------------------------------------------------------*/

void handleComputeNothing()
{    
}

/*---------------------------------------------------------------------------*/

void handleMouse(int button, int state, int cx, int cy)
{
    printf("\nNo mouse handler yet.");
}


/*---------------------------------------------------------------------------*/

int main (int argc, char* argv[])
{
/* ---------- set boundary and grid size ----------------------------------- */

g_bx=1;
g_by=1;

g_hx=1;
g_hy=1;

/* ---------- read in arguments -------------------------------------------- */

 
if (argc==4) 
 {   
     strcpy(g_image1,argv[1]);
     strcpy(g_image2,argv[2]);                
     g_g_pixel_ratio=(int)atoi(argv[3]);
 }
if (argc==5) 
 {            
     strcpy(g_image1,argv[1]);
     strcpy(g_image2,argv[2]);               
     g_g_pixel_ratio=(int)atoi(argv[3]);
     strcpy(g_ref,argv[4]);    
 }
if ((argc<4)||(argc>5))
 {
   printf("\n\n%s image1.pgm image2.pgm zoom_ratio [ground_truth.F]\n\n",argv[0]);
     exit(0);    
 }


/* set default parameters */

g_g_max_disp=2.0;

g_direct_compute=0;

g_e_sigma=0.0;

g_m_alpha=100;

g_n_iter=100;

g_active_param=1000;

g_u_ref=NULL;
g_v_ref=NULL;


/* ---------- read in information of first image ------------------- */

read_pgm_header(g_image1,&g_position,&g_nx,&g_ny);


/* ---------- memory allocation ------------------------------------ */

ALLOC_MATRIX(2, g_nx+2*g_bx, g_ny+2*g_by, &g_f1,   &g_f2);
ALLOC_MATRIX(2, g_nx+2*g_bx, g_ny+2*g_by, &g_f1_s, &g_f2_s);
ALLOC_MATRIX(2, g_nx+2*g_bx, g_ny+2*g_by, &g_u,    &g_v);

if (argc==5)
ALLOC_MATRIX(2, g_nx+2*g_bx, g_ny+2*g_by, &g_u_ref,   &g_v_ref);

ALLOC_CUBIX(1, 2*g_nx+2*g_bx, g_ny+2*g_by, 3, &g_p6);
g_pixels = (GLubyte *) malloc (2*(g_nx+1)*g_ny*3*sizeof(GLubyte));   


/* ---------- read in image pair ------------------------------------------- */

/* read first frame */
read_pgm_header(g_image1,&g_position,&g_nx,&g_ny);
read_pgm_data(g_image1,g_position,g_f1,g_nx,g_ny,g_bx,g_by);

/* read second frame */ 
read_pgm_header(g_image2,&g_position,&g_nx,&g_ny);
read_pgm_data(g_image2,g_position,g_f2,g_nx,g_ny,g_bx,g_by);


/* ---------- read in ground truth displacement field ---------------------- */

if (argc==5)
  read_barron_data(g_ref,g_u_ref,g_v_ref,g_nx,g_ny,g_bx,g_by); 


/* ---------- M A I N   L O O P -------------------------------------------- */

// open OpenGL window */
glutInit(&argc, argv);
glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);

glutInitWindowSize((int)round(2*g_nx*g_g_pixel_ratio),
		   (int)round(g_ny*g_g_pixel_ratio));
glutCreateWindow("CORRESPONDENCE PROBLEMS - VISUALISATION FRONTEND");

// register handle routines
glutDisplayFunc(handleDraw);
glutIdleFunc(handleComputeNothing);
glutKeyboardFunc(handleKeyboard);
glutSpecialFunc(handleKeyboardspecial);
glutMouseFunc(handleMouse);

// main 
handleComputeDisplacements();
handleDraw();
showParams();
glutMainLoop();


/* ---------- free memory -------------------------------------------------- */

FREE_MATRIX (2, g_nx+2*g_bx, g_ny+2*g_by, g_f1,   g_f2);
FREE_MATRIX (2, g_nx+2*g_bx, g_ny+2*g_by, g_f1_s, g_f2_s);
FREE_MATRIX (2, g_nx+2*g_bx, g_ny+2*g_by, g_u,    g_v);

if (argc==5)
FREE_MATRIX(2,   g_nx+2*g_bx, g_ny+2*g_by, g_u_ref, g_v_ref);

FREE_CUBIX (1, 2*g_nx+2*g_bx, g_ny+2*g_by, 3, g_p6);
free(g_pixels);
	
printf("\n\n\n");


return(0);
}








