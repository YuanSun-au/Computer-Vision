#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "libraries.h"

/*************************/

/*--------------------------------------------------------------------------*/
inline float max(float a, float b) { return (a > b) ? a : b; }
/*--------------------------------------------------------------------------*/
inline float min(float a, float b) { return (a < b) ? a : b; }
/*--------------------------------------------------------------------------*/
inline float sq(float a) { return a * a; }
/*--------------------------------------------------------------------------*/

void brightness(
    float **u,
    int nx,
    int ny,
    float offset)
{
    for (int i = 1; i <= nx; i++)
        for (int j = 1; j <= ny; j++)
        {
            u[i][j] += offset;
        }
}

/*--------------------------------------------------------------------------*/

void contrast(
    float **u,
    int nx,
    int ny,
    float scaling)
{
    for (int i = 1; i <= nx; i++)
        for (int j = 1; j <= ny; j++)
        {
            u[i][j] -= 127.5;
            /* scale all pixel values by "scaling" */
            u[i][j] *= scaling;
            u[i][j] += 127.5;
        }
}

/*--------------------------------------------------------------------------*/

int main(int argc, char *argv[])
{
    char in[80];
    char out[80];
    char row[80];

    int nx, ny;

    float **u1, **u2, **u3;

    int i, j;

    unsigned char byte;

    long int selection;

    float offset, scaling;

    long pos;

    int color;

    FILE *ptr;

    /* get information: color image or grayscale? */
    printf("Color [0=No, 1=Yes]:                                    ");
    gets(row);
    sscanf(row, "%d", &color);

    /* get input image name */
    printf("input image:                                            ");
    gets(in);

    /* read image data */
    if (color == 0)
    {
        read_pgm_header(in, &pos, &nx, &ny);
        alloc_matrix(&u1, nx + 2, ny + 2);
        read_pgm_data(in, pos, u1, nx, ny, 1, 1);
    }
    else
    {
        read_ppm_header(in, &pos, &nx, &ny);
        alloc_matrix(&u1, nx + 2, ny + 2);
        alloc_matrix(&u2, nx + 2, ny + 2);
        alloc_matrix(&u3, nx + 2, ny + 2);
        read_ppm_data_channelwise(in, pos, u1, u2, u3, nx, ny, 1, 1);
    }

    /* get information: what operation should be performed? */
    printf("\n1 = Brightness, 2 = Contrast:                         ");
    gets(row);
    sscanf(row, "%ld", &selection);

    /* check if selection is valid */
    if (selection < 1 || selection > 2)
    {
        printf("Wrong selection\n");
        exit(0);
    }

    /* if selection*/
    if (selection == 1)
    {
        /* get brightness offset */
        printf("Brightness offset:                                      ");
        gets(row);
        sscanf(row, "%f", &offset);
    }
    else if (selection == 2)
    {
        /* get scaling */
        printf("Scaling:                                                ");
        gets(row);
        sscanf(row, "%f", &scaling);
    }

    /* get output image name */
    printf("output image:                                           ");
    gets(out);

    /* perform selected operation */
    switch (selection)
    {
    case 1:
        brightness(u1, nx, ny, offset);
        if (color)
            brightness(u2, nx, ny, offset);
        if (color)
            brightness(u3, nx, ny, offset);
        break;
    case 2:
        contrast(u1, nx, ny, scaling);
        if (color)
            contrast(u2, nx, ny, scaling);
        if (color)
            contrast(u3, nx, ny, scaling);
        break;
    default:
        printf("Wrong selection\n");
        exit(0);
        break;
    }

    /* open file and write header (incl. filter parameters) */
    ptr = fopen(out, "wb");
    fprintf(ptr, "P%d \n", 5 + color); //P5 for gray, P6 for color images
    fprintf(ptr, "%d %d \n255\n", nx, ny);

    fclose(ptr);

    /* write data */
    if (!color)
        write_pgm_data(out, u1, nx, ny, 1, 1);
    else
        write_ppm_data_channelwise(out, u1, u2, u3, nx, ny, 1, 1);

    /* free memory */
    disalloc_matrix(u1, nx + 2, ny + 2);
    /* u2 and u3 for color images only */
    if (color)
        disalloc_matrix(u2, nx + 2, ny + 2);
    if (color)
        disalloc_matrix(u3, nx + 2, ny + 2);

    return (0);
}