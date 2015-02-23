#include <stdio.h>
#include <stdlib.h>
#include "voz.h"

/* Open File */
/* Returns number of particles read */
long openfile(char *filename, FILE **f) {
    long np;

    *f = fopen(filename, "r");
    if (f == NULL) {
        printf("Unable to open position file %s\n\n", filename);
        exit(0);
    }

    /* Read number of particles */
    fread(&np, sizeof(long), 1, *f);

    return np;
}

/*
Read in a chunk of particles
Return number actually read in
np is the total number of particles.
nread is the offset.
*/
long
posread_chunk(FILE *f, realT **p, realT fact, long np, long nread)
{
    long i, j, d;
    long ntoread = N_CHUNK;
    realT *ptemp;

    /* fix read size if near the EOF. */
    if ((np - nread) < N_CHUNK) {
        ntoread = np - nread;
    }

    /*
    TODO
    This function is typically called from inside a loop.
    Maybe we shouldn't allocate a buffer every call...
    */
    ptemp = (realT *) malloc(3 * ntoread * sizeof(realT));

    /* seek to first particle */
    fseek(f, sizeof(long) + 3 * nread * sizeof(realT), SEEK_SET);
    /* read buffer size */
    fread(ptemp, sizeof(realT), 3 * ntoread, f);
    /* save to other format */
    j = 0;
    for (i = 0; i < ntoread; i++) {
        p[i][0] = ptemp[j++];
        p[i][1] = ptemp[j++];
        p[i][2] = ptemp[j++];
    }

    free(ptemp);

    for (i = 0; i < ntoread; i++) DL p[i][d] *= fact;

    return ntoread;
}

long posread_isol(char *posfile, float ***p, float fact, long *np, long *npreal) {

    FILE *pos;
    long dum,d,i;
    float xt,yt,zt;
    long npread, npreadreal;
    float xmin,xmax,ymin,ymax,zmin,zmax;

    pos = fopen(posfile, "r");
    if (pos == NULL) {
        printf("Unable to open position file %s\n\n",posfile);
        exit(0);
    }

    /* Read number of particles */
    if (fscanf(pos,"%d%d\n",&npread,&npreadreal) != 2) {
        printf("Problem reading no. particles\n");
        exit(0);
    }

    /* Allocate the arrays */
    (*p) = (float **)malloc(npread*sizeof(float *));
    /* Fill the arrays */
    for (i=0; (fscanf(pos,"%f%f%f\n",&xt,&yt,&zt) == 3); i++) {
      if (i >= npread) {
        printf("More particles in %s than first line claims!  Exiting.\n\n",posfile);
        exit(0);
      }
      (*p)[i] = (float *)malloc(3*sizeof(float));
      if ((*p)[i] == NULL) {
        printf("Unable to allocate particle array in readfiles!\n");
        fflush(stdout);
        exit(0);
      }
      (*p)[i][0]=xt*fact; (*p)[i][1]=yt*fact; (*p)[i][2] = zt*fact;
    }
    fclose(pos);
    printf("%d\n",i);
    if (npread != i) {
      printf("Read %d particles (not %d, as the first line claims)!  Exiting.\n\n", i, npread);
      exit(0);
    }
    /* Test range -- can comment out */
    xmin = BF; xmax = -BF; ymin = BF; ymax = -BF; zmin = BF; zmax = -BF;
    for (i=0; i<npread;i++) {
      if ((*p)[i][0]<xmin) xmin = (*p)[i][0]; if ((*p)[i][0]>xmax) xmax = (*p)[i][0];
      if ((*p)[i][1]<ymin) ymin = (*p)[i][1]; if ((*p)[i][1]>ymax) ymax = (*p)[i][1];
      if ((*p)[i][2]<zmin) zmin = (*p)[i][2]; if ((*p)[i][2]>zmax) zmax = (*p)[i][2];
    }
    printf("npread: %d, x: %f,%f; y: %f,%f; z: %f,%f\n",npread,xmin,xmax, ymin,ymax, zmin,zmax); fflush(stdout);

    *np = npread;
    *npreal = npreadreal;
    return(0);
}

/* Positions */
/* Returns number of particles read */
long posread(char *posfile, realT ***p, realT fact) {

    FILE *pos;
    long npr,dum,d,i;
    realT xmin,xmax,ymin,ymax,zmin,zmax;
    realT *ptemp;

    pos = fopen(posfile, "rb");
    if (pos == NULL) {
      printf("Unable to open position file %s\n\n",posfile);
      exit(0);
    }
    /* Fortran77 4-byte headers and footers */
    /* Delete "dum" statements if you don't need them */

    /* Read number of particles */
     fread(&npr,sizeof(long),1,pos);

    /* Allocate the arrays */
    (*p) = (realT **)malloc(npr*sizeof(realT *));
    ptemp = (realT *)malloc(npr*sizeof(realT));

    printf("np = %d\n",npr);

    /* Fill the arrays */

    for (i=0; i<npr; i++) {
      (*p)[i] = (realT *)malloc(3*sizeof(realT));
      if ((*p)[i] == NULL) {
        printf("Unable to allocate particle array in readfiles!\n");
        fflush(stdout);
        exit(0);
      }
    }

    fread(ptemp,sizeof(realT),npr,pos);
    for (i=0; i<npr; i++) (*p)[i][0] = ptemp[i];

    fread(ptemp,sizeof(realT),npr,pos);
    for (i=0; i<npr; i++) (*p)[i][1] = ptemp[i];

    fread(ptemp,sizeof(realT),npr,pos);
    for (i=0; i<npr; i++) (*p)[i][2] = ptemp[i];


    fclose(pos);
    free(ptemp);

    /* Get into physical units (Mpc/h) */

    for (i=0; i<npr; i++) DL (*p)[i][d] *= fact;


    /* Test range -- can comment out */
    xmin = BF; xmax = -BF; ymin = BF; ymax = -BF; zmin = BF; zmax = -BF;
    for (i=0; i<npr;i++) {
      if ((*p)[i][0]<xmin) xmin = (*p)[i][0]; if ((*p)[i][0]>xmax) xmax = (*p)[i][0];
      if ((*p)[i][1]<ymin) ymin = (*p)[i][1]; if ((*p)[i][1]>ymax) ymax = (*p)[i][1];
      if ((*p)[i][2]<zmin) zmin = (*p)[i][2]; if ((*p)[i][2]>zmax) zmax = (*p)[i][2];
    }
    printf("np: %d, x: %g,%g; y: %g,%g; z: %g,%g\n",npr,xmin,xmax, ymin,ymax, zmin,zmax); fflush(stdout);

    return(npr);
}

/* Velocities */
/* Returns number of particles read */
long velread(char *velfile, realT ***v, realT fact) {

    FILE *vel;
    long npr,dum,d,i;
    realT xmin,xmax,ymin,ymax,zmin,zmax;

    vel = fopen(velfile, "rb");
    if (vel == NULL) {
      printf("Unable to open velocity file %s\n\n",velfile);
      exit(0);
    }
    /* Fortran77 4-byte headers and footers */
    /* Delete "dum" statements if you don't need them */

    /* Read number of particles */
     fread(&npr,sizeof(long),1,vel);

    /* Allocate the arrays */
    (*v) = (realT **)malloc(3*sizeof(realT*));
    for (i=0;i<3;i++) (*v)[i] = (realT *)malloc(npr*sizeof(realT));

    /* Fill the arrays */
    fread((*v)[0],sizeof(realT),npr,vel);
    fread((*v)[1],sizeof(realT),npr,vel);
    fread((*v)[2],sizeof(realT),npr,vel);

    fclose(vel);

    /* Convert from code units into physical units (km/sec) */

    for (i=0; i<npr; i++) DL (*v)[d][i] *= fact;

    /* Test range -- can comment out */
    xmin = BF; xmax = -BF; ymin = BF; ymax = -BF; zmin = BF; zmax = -BF;
    for (i=0; i<npr;i++) {
      if ((*v)[0][i] < xmin) xmin = (*v)[0][i]; if ((*v)[0][i] > xmax) xmax = (*v)[0][i];
      if ((*v)[1][i] < ymin) ymin = (*v)[1][i]; if ((*v)[1][i] > ymax) ymax = (*v)[1][i];
      if ((*v)[2][i] < zmin) zmin = (*v)[2][i]; if ((*v)[2][i] > zmax) zmax = (*v)[2][i];
    }
    printf("vx: %g,%g; vy: %g,%g; vz: %g,%g\n",xmin,xmax, ymin,ymax, zmin,zmax);

    return(npr);
}
