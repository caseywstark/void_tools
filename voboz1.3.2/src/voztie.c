#include <stdio.h>
#include <stdlib.h>
#include "voz.h"

#define EQUALTHRESHOLD 1.52587890625e-5 /* 2^-16 */

int
main(int argc, char *argv[])
{

    FILE *part, *adj, *vol;
    char partfile[80], *suffix, adjfile[80], volfile[80];
    realT *vols, volstemp;

    PARTADJ *adjs;

    long numdiv, np, np2;
    int na;

    long i, j, k, p, nout;
    long nvp, npnotdone, nvpmax, nvpsum;
    int *orig;
    realT avgnadj, avgvol;

    if (argc != 3) {
        printf("Wrong number of arguments.\n");
        printf("arg1: number of divisions (default 2)\n");
        printf("arg2: suffix describing this run\n\n");
        exit(0);
    }

    numdiv = atoi(argv[1]);
    suffix = argv[2];

    if (numdiv < 2) {
        printf("Cannot have a number of divisions less than 2.  Resetting to 2:\n");
        numdiv = 2;
    }

    np = -1; nvpmax = -1; nvpsum = 0;

    for (i = 0; i < numdiv; i++) {
        for (j = 0; j < numdiv; j++) {
            for (k = 0; k < numdiv; k++) {
                sprintf(partfile, "part.%s.%02d.%02d.%02d", suffix, (int)i, (int)j, (int)k);
                part = fopen(partfile, "r");
                if (part == NULL) {
                    printf("Unable to open file %s.\n\n", partfile);
                    exit(0);
                }
                fread(&np2, 1, sizeof(long), part);
                fread(&nvp, 1, sizeof(long), part);
                if (np == -1)
                    np = np2;
                else if (np2 != np) {
                    printf("Incompatible total particle numbers: %ld, %ld\n\n", np, np2);
                    exit(0);
                }
                if (nvp > nvpmax) nvpmax = nvp;
                fclose(part);
            }
        }
    }

    printf("We have %ld particles to tie together.\n", np); fflush(stdout);
    printf("The maximum number of particles in a file is %ld.\n", nvpmax);

    adjs = (PARTADJ *)malloc(np * sizeof(PARTADJ));
    if (adjs == NULL) printf("Couldn't allocate adjs.\n");
    vols = (realT *)malloc(np * sizeof(realT));
    if (vols == NULL) printf("Couldn't allocate vols.\n");
    orig = (int *) malloc(nvpmax * sizeof(int));
    if (orig == NULL) printf("Couldn't allocate orig.\n");
    if ((vols == NULL) || (orig == NULL) || (adjs == NULL)) {
        printf("Not enough memory to allocate. Exiting.\n");
        exit(0);
    }
    for (p = 0; p < np; p++)
        vols[p] = -1.;

    for (i = 0; i < numdiv; i++) {
        for (j = 0; j < numdiv; j++) {
            for (k = 0; k < numdiv; k++) {
                sprintf(partfile, "part.%s.%02d.%02d.%02d", suffix, (int)i, (int)j, (int)k);
                part = fopen(partfile, "r");
                if (part == NULL) {
                    printf("Unable to open file %s.\n\n",partfile);
                    exit(0);
                }
                fread(&np2, sizeof(long), 1, part);
                fread(&nvp, sizeof(long), 1, part);
                nvpsum += nvp;

                fread(orig, sizeof(int), nvp, part);
                for (p = 0; p < nvp; p++) {
                    fread(&volstemp, sizeof(realT), 1, part);
                    if (vols[orig[p]] > -1.0)
                        if (fabs(vols[orig[p]] - volstemp)/volstemp < EQUALTHRESHOLD) {
                            printf("Warning: different vols measured for p.%d (%g,%g). Ignore if close enough.\n",
                                orig[p], vols[orig[p]], volstemp);
                            volstemp = 0.5 * (volstemp + vols[orig[p]]);
                        }
                    vols[orig[p]] = volstemp;
                }

                for (p = 0; p < nvp; p++) {
                    fread(&na, sizeof(int), 1, part);
                    if (na > 0) {
                        adjs[orig[p]].nadj = na;
                        adjs[orig[p]].adj = (int *)malloc(na * sizeof(int));
                        if (adjs[orig[p]].adj == NULL) {
                            printf("Couldn't allocate adjs[orig[%ld]].adj (size %i).\n", p, na);
                            exit(0);
                        }
                        fread(adjs[orig[p]].adj, sizeof(int), na, part);
                    } else {
                        printf("0"); fflush(stdout);
                    }
                }

                fclose(part);
                printf("%d ", k);
            }
        }
    }

    printf("\n");
    npnotdone = 0; avgnadj = 0.; avgvol = 0.;
    for (p = 0; p < np; p++) {
        if (vols[p] == -1.) npnotdone++;
        avgnadj += (realT)(adjs[p].nadj);
        avgvol += (realT)(vols[p]);
    }
    if (npnotdone > 0)
        printf("%ld particles not done!\n", npnotdone);

    printf("%ld particles done more than once.\n", nvpsum - np);
    avgnadj /= (realT)np;
    avgvol /= (realT)np;
    printf("Average # adjacencies = %g (%f for Poisson)\n", avgnadj,
        48.*3.141593*3.141593/35.+2.);
    printf("Average volume = %g\n", avgvol);

    /* Now the output! */

    sprintf(adjfile, "%s.adj", suffix);
    sprintf(volfile, "%s.vol", suffix);

    printf("Outputting to %s, %s\n\n", adjfile, volfile);

    adj = fopen(adjfile, "w");
    if (adj == NULL) {
        printf("Unable to open %s\n",adjfile);
        exit(0);
    }
    fwrite(&np, sizeof(long), 1, adj);
    /* Adjacencies: first the numbers of adjacencies,
       and the number we're actually going to write per particle */
    for (i = 0; i < np; i++)
        fwrite(&adjs[i].nadj, sizeof(int), 1, adj);

    /* Now the lists of adjacencies (without double counting) */
    for (i = 0; i < np; i++)
        if (adjs[i].nadj > 0) {
            nout = 0;
            for (j = 0; j < adjs[i].nadj; j++) if (adjs[i].adj[j] > i) nout++;
            fwrite(&nout, sizeof(int), 1, adj);
            for (j = 0; j < adjs[i].nadj; j++)
                if (adjs[i].adj[j] > i)
                    fwrite(&(adjs[i].adj[j]), sizeof(int), 1, adj);
        }

    fclose(adj);

    /* Volumes */
    vol = fopen(volfile,"w");
    if (vol == NULL) {
        printf("Unable to open %s\n",volfile);
        exit(0);
    }
    fwrite(&np, sizeof(long), 1, vol);
    fwrite(vols, sizeof(realT), np, vol);

    fclose(vol);

    return(0);
}
