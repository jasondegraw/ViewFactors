/******************************************************************************/
/*                                                                            */
/*  Compute the view factors of a box-           _____________________        */
/*  shaped enclosure                            /  /  /  /  /  /  /  /|       */
/*                                             /__/__/__/__/__/__/__/ |       */
/*  The view factors are stored in a          /  /  /  /  /  /  /  /| |       */
/*  vector of vectors - which is a           /__/__/__/__/__/__/__/ |/|       */
/*  little inconvenient but will do         /  /  /  /  /  /  /  /| | |       */
/*  for now.                               /__/__/__/__/__/__/__/ |/|/|       */
/*                                         |  |  |  |  |  |  |  | | | |       */
/*                                         |__|__|__|__|__|__|__|/|/|/|       */
/*                                         |  |  |  |  |  |  |  | | | |       */
/*                                         |__|__|__|__|__|__|__|/|/|/        */
/*                                         |  |  |  |  |  |  |  | | /         */
/*                                         |__|__|__|__|__|__|__|/|/          */
/*                                         |  |  |  |  |  |  |  | /           */
/*                                         |__|__|__|__|__|__|__|/            */
/*                                                                            */
/*                                                                            */
/******************************************************************************/
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
    double **vf;
    int i,j,ij;
    int nx = 2;
    int ny = 2;
    int nz = 2;
    double dx = 1;
    double dy = 1;
    double dz = 1;
    double X = 2;
    double Y = 2;
    double Z = 2;
    int N = 2*nx*ny + 2*nx*nz + 2*ny*nz;
    double *base_x,*base_y,*base_z;
    vf = calloc(N,sizeof(double));
    base_x = calloc(3*N,sizeof(double));
    base_y = base_x + N;
    base_z = base_y + N;
    ij = 0;
    /* z = 0 plane */
    for(j=0;j<ny;j++)
    {
        for(i=0;i<nx;i++)
        {
            base_x[ij] = i*dx;
            base_y[ij] = j*dy;
            base_z[ij] = 0.0;
        }
    }
    /* z = Z plane */
    for(j=0;j<ny;j++)
    {
        for(i=0;i<nx;i++)
        {
            base_x[ij] = i*dx;
            base_y[ij] = j*dy;
            base_z[ij] = Z;
        }
    }
    /* y = 0 plane */
    for(j=0;j<nz;j++)
    {
        for(i=0;i<nx;i++)
        {
            base_x[ij] = i*dx;
            base_y[ij] = 0.0;
            base_z[ij] = j*dz;
        }
    }
    /* z = Z plane */
    for(j=0;j<ny;j++)
    {
        for(i=0;i<nx;i++)
        {
            base_x[ij] = i*dx;
            base_y[ij] = Y;
            base_z[ij] = j*dz;
        }
    }
    /* x = 0 plane */
    for(j=0;j<nz;j++)
    {
        for(i=0;i<ny;i++)
        {
            base_x[ij] = 0.0;
            base_y[ij] = i*dy;
            base_z[ij] = j*dz;
        }
    }
    /* x = X plane */
    for(j=0;j<nz;j++)
    {
        for(i=0;i<ny;i++)
        {
            base_x[ij] = X;
            base_y[ij] = i*dy;
            base_z[ij] = j*dz;
        }
    }
    return EXIT_SUCCESS;
}
