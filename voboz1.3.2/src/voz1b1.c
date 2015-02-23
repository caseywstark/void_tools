#include "voz.h"

int delaunadj (coordT *points, int nvp, int nvpbuf, int nvpall, PARTADJ **adjs);
int vorvol (coordT *deladjs, coordT *points, pointT *intpoints, int numpoints, realT *vol);

long openfile(char *filename, FILE **f);
long posread_chunk(FILE *f, realT **p, realT fact, long np, long nread);

void voz1b1(char *posfile, realT border, realT boxsize,
    int numdiv, int b[], char *suffix){
    int exitcode;
    long i, j, np, np_current, np_tot;
    realT **r;
    coordT rtemp[3], rtemp_inner[3], *parts;
    coordT deladjs[3*MAXVERVER], points[3*MAXVERVER];
    pointT intpoints[3*MAXVERVER];
    FILE *pos, *out;
    char outfile[80];
    PARTADJ *adjs;
    realT *vols;
    realT predict, xmin,xmax,ymin,ymax,zmin,zmax;
    int *orig;
    int isitinbuf;
    int n_overlap;
    int isitinmain, d, i_x, i_y, i_z;
    long nvp, nvpall, nvpbuf;
    realT width, width2, totwidth, totwidth2, bf, s, g;
    realT c[3];
    realT totalvol;
    realT overlap_corner;

    /* Boxsize should be the range in r, yielding a range 0-1 */
    /* np = posread(posfile,&r,1./boxsize); */
    np = openfile(posfile, &pos);

    printf("%ld particles\n", np); fflush(stdout);

    /* chunked read doesn't allocate memory */
    r = (realT **)malloc(N_CHUNK * sizeof(realT *));
    for (i = 0; i < N_CHUNK; i++) {
        r[i] = (realT *)malloc(3 * sizeof(realT));
        if (r[i] == NULL) {
            printf("Unable to allocate particle array!\n");
            fflush(stdout);
            exit(0);
        }
    }

    xmin = BF; xmax = -BF; ymin = BF; ymax = -BF; zmin = BF; zmax = -BF;

    /*
    purely informational removing to save file reads
       for (i=0; i<np;i++) {
       if (r[i][0]<xmin) xmin = r[i][0]; if (r[i][0]>xmax) xmax = r[i][0];
       if (r[i][1]<ymin) ymin = r[i][1]; if (r[i][1]>ymax) ymax = r[i][1];
       if (r[i][2]<zmin) zmin = r[i][2]; if (r[i][2]>zmax) zmax = r[i][2];
       }
       printf("np: %d, x: %g,%g; y: %g,%g; z: %g,%g\n",np,xmin,xmax, ymin,ymax, zmin,zmax); fflush(stdout);
    */

    width = 1.0 / (realT)numdiv;
    width2 = 0.5 * width;
    if (border > 0.) bf = border;
    else bf = 0.1;

    /* In units of 0-1, the thickness of each subregion's buffer*/
    totwidth = width + 2.0 * bf;
    totwidth2 = width2 + bf;

    s = width / (realT)NGUARD;
    if ( bf * bf - 2.0 * s * s < 0.0) {
        printf("bf = %g, s = %g.\n", bf, s);
        printf("Not enough guard points for given border.\nIncrease guards to >= %g\n.",
            sqrt(2.)*width/bf);
        exit(0);
    }
    g = (bf / 2.) * (1. + sqrt(1 - 2.*s*s/(bf*bf)));
    printf("s = %g, bf = %g, g = %g.\n",s,bf,g);
    fflush(stdout);

    DL c[d] = ((realT)b[d] + 0.5) * width;

    /* inner edge of overlapping buffer region */
    overlap_corner = 1.0 - totwidth2;

    printf("c: %g, %g, %g\n", c[0], c[1], c[2]);
    /* Assign temporary array*/
    nvpbuf = 0; /* Number of particles to tesselate, including buffer */
    nvp = 0; /* Without the buffer */

    np_tot = 0;
    while (np_tot < np) {
        np_current = posread_chunk(pos, r, 1./boxsize, np, np_tot);

        /* loop over current buffer particles */
        for (i = 0; i < np_current; i++) {
            n_overlap = 0;
            isitinbuf = 1;
            isitinmain = 1;

            DL {
                rtemp[d] = (realT)r[i][d] - (realT)c[d];
                if (rtemp[d] > 0.5) rtemp[d] --;
                if (rtemp[d] < -0.5) rtemp[d] ++;
                isitinbuf = isitinbuf && (fabs(rtemp[d]) < totwidth2);
                isitinmain = isitinmain && (fabs(rtemp[d]) <= width2);
            }

            if (isitinmain){
                nvp++;
                nvpbuf++;
            }
            if (isitinbuf && !isitinmain) {
                /* handle the pathalogical case, when buffer zones overlap */
                for (i_x = 0; i_x < 2; i_x++) {
                    rtemp_inner[0] = rtemp[0];

                    if (i_x) {
                        if (rtemp_inner[0] > 0.0)  rtemp_inner[0] --;
                        else rtemp_inner[0] ++;
                    }
                    if (fabs(rtemp_inner[0]) < totwidth2){
                        for (i_y = 0; i_y < 2; i_y++) {
                            rtemp_inner[1] = rtemp[1];
                            if (i_y) {
                                if (rtemp_inner[1] > 0.0)  rtemp_inner[1] --;
                                else rtemp_inner[1] ++;
                            }
                            if (fabs(rtemp_inner[1]) < totwidth2) {
                                for (i_z = 0; i_z < 2; i_z++) {
                                    rtemp_inner[2] = rtemp[2];
                                    if (i_z) {
                                        if (rtemp_inner[2] > 0.0) rtemp_inner[2] --;
                                        else rtemp_inner[2] ++;
                                    }
                                    if (fabs(rtemp_inner[2]) < totwidth2) {
                                        nvpbuf++;
                                    }
                                }
                            }
                        }
                    }
                }
            }

        }

        /* update counter */
        np_tot += np_current;
        printf("%.4f\n", (double)np_tot / np); fflush(stdout);
    }

    printf("nvp = %ld, after count\n", nvp);
    printf("nvpbuf = %ld, after count\n", nvpbuf);

    nvpbuf += 6 * (NGUARD + 1) * (NGUARD + 1); /* number of guard points */

    parts = (coordT *)malloc(3 * nvpbuf * sizeof(coordT));
    orig = (int *)malloc(nvpbuf * sizeof(int));

    if (parts == NULL) {
        printf("Unable to allocate parts\n");
        fflush(stdout);
    }
    if (orig == NULL) {
        printf("Unable to allocate orig\n");
        fflush(stdout);
    }

    nvp = 0; nvpall = 0; /* nvp = number of particles without buffer */
    xmin = BF; xmax = -BF; ymin = BF; ymax = -BF; zmin = BF; zmax = -BF;

    np_tot = 0;
    while (np_tot < np) {
        np_current = posread_chunk(pos, r, 1./boxsize, np, np_tot);

        for (i = 0; i < np_current; i++) {

            isitinmain = 1;

            DL {
                rtemp[d] = r[i][d] - c[d];
                if (rtemp[d] > 0.5) rtemp[d] --;
                if (rtemp[d] < -0.5) rtemp[d] ++;
                isitinmain = isitinmain && (fabs(rtemp[d]) <= width2);
            }

            if (isitinmain) {
                parts[3*nvp] = rtemp[0];
                parts[3*nvp+1] = rtemp[1];
                parts[3*nvp+2] = rtemp[2];
                orig[nvp] = np_tot + i;
                nvp++;
                if (rtemp[0] < xmin) xmin = rtemp[0];
                if (rtemp[0] > xmax) xmax = rtemp[0];
                if (rtemp[1] < ymin) ymin = rtemp[1];
                if (rtemp[1] > ymax) ymax = rtemp[1];
                if (rtemp[2] < zmin) zmin = rtemp[2];
                if (rtemp[2] > zmax) zmax = rtemp[2];
            }
        }
        np_tot += np_current;
    }

    printf("nvp = %ld\n", nvp);
    printf("x: %g,%g; y: %g,%g; z:%g,%g\n",xmin,xmax,ymin,ymax,zmin,zmax);
    nvpbuf = nvp;

    np_tot = 0;
    while(np_tot < np){
        np_current = posread_chunk(pos, r, 1./boxsize, np, np_tot);
        for (i=0; i< np_current; i++) {

            isitinbuf = 1;
            DL {
                rtemp[d] = r[i][d] - c[d];
                if (rtemp[d] > 0.5) rtemp[d] --;
                if (rtemp[d] < -0.5) rtemp[d] ++;
                isitinbuf = isitinbuf && (fabs(rtemp[d])<totwidth2);
            }
            if ( (isitinbuf > 0) &&
                 ((fabs(rtemp[0])>width2)||(fabs(rtemp[1])>width2)
                  ||(fabs(rtemp[2])>width2))) {

                /* handle the pathalogical case, when buffer zones overlap */
                for (i_x = 0; i_x < 2; i_x++) {
                    rtemp_inner[0] = rtemp[0];
                    if (i_x) {
                        if (rtemp_inner[0] > 0.0) rtemp_inner[0] --;
                        else rtemp_inner[0] ++;
                    }
                    if(fabs(rtemp_inner[0]) < totwidth2){
                        for( i_y = 0; i_y < 2; i_y++){
                            rtemp_inner[1] = rtemp[1];
                            if(i_y){
                                if (rtemp_inner[1] > 0.0) rtemp_inner[1] --;
                                else rtemp_inner[1] ++;
                            }
                            if(fabs(rtemp_inner[1]) < totwidth2){
                                for( i_z = 0; i_z < 2; i_z++){
                                    rtemp_inner[2] = rtemp[2];
                                    if(i_z){
                                        if (rtemp_inner[2] > 0.0) rtemp_inner[2] --;
                                        else rtemp_inner[2] ++;
                                    }
                                    if(fabs(rtemp_inner[2]) < totwidth2){
                                        parts[3*nvpbuf] = rtemp_inner[0];
                                        parts[3*nvpbuf+1] = rtemp_inner[1];
                                        parts[3*nvpbuf+2] = rtemp_inner[2];
                                        orig[nvpbuf] = np_tot + i;

                                        nvpbuf++;
                                        if (rtemp_inner[0] < xmin) xmin = rtemp_inner[0];
                                        if (rtemp_inner[0] > xmax) xmax = rtemp_inner[0];
                                        if (rtemp_inner[1] < ymin) ymin = rtemp_inner[1];
                                        if (rtemp_inner[1] > ymax) ymax = rtemp_inner[1];
                                        if (rtemp_inner[2] < zmin) zmin = rtemp_inner[2];
                                        if (rtemp_inner[2] > zmax) zmax = rtemp_inner[2];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        np_tot += np_current;
    }

    printf("nvp = %ld, after assignment\n", nvp);
    printf("nvpbuf = %ld, after assignment\n", nvpbuf);
    printf("x: %g,%g; y: %g,%g; z:%g,%g\n",xmin,xmax,ymin,ymax,zmin,zmax);
    nvpall = nvpbuf;
    predict = pow(width + 2.0 * bf, 3) * np;
    printf("There should be ~ %g points; there are %ld\n", predict, nvpbuf);

    for (i = 0; i < N_CHUNK; i++) free(r[i]);
    free(r);
    fclose(pos);

    /* Add guard points */
    for (i=0; i<NGUARD+1; i++) {
        for (j=0; j<NGUARD+1; j++) {
            /* Bottom */
            parts[3*nvpall]   = -width2 + (realT)i * s;
            parts[3*nvpall+1] = -width2 + (realT)j * s;
            parts[3*nvpall+2] = -width2 - g;
            nvpall++;
            /* Top */
            parts[3*nvpall]   = -width2 + (realT)i * s;
            parts[3*nvpall+1] = -width2 + (realT)j * s;
            parts[3*nvpall+2] = width2 + g;
            nvpall++;
        }
    }
    for (i=0; i<NGUARD+1; i++) { /* Don't want to overdo the corners*/
        for (j=0; j<NGUARD+1; j++) {
            parts[3*nvpall]   = -width2 + (realT)i * s;
            parts[3*nvpall+1] = -width2 - g;
            parts[3*nvpall+2] = -width2 + (realT)j * s;
            nvpall++;

            parts[3*nvpall]   = -width2 + (realT)i * s;
            parts[3*nvpall+1] = width2 + g;
            parts[3*nvpall+2] = -width2 + (realT)j * s;
            nvpall++;
        }
    }
    for (i=0; i<NGUARD+1; i++) {
        for (j=0; j<NGUARD+1; j++) {
            parts[3*nvpall]   = -width2 - g;
            parts[3*nvpall+1] = -width2 + (realT)i * s;
            parts[3*nvpall+2] = -width2 + (realT)j * s;
            nvpall++;

            parts[3*nvpall]   = width2 + g;
            parts[3*nvpall+1] = -width2 + (realT)i * s;
            parts[3*nvpall+2] = -width2 + (realT)j * s;
            nvpall++;
        }
    }
    xmin = BF; xmax = -BF; ymin = BF; ymax = -BF; zmin = BF; zmax = -BF;
    for (i=nvpbuf;i<nvpall;i++) {
        if (parts[3*i] < xmin) xmin = parts[3*i];
        if (parts[3*i] > xmax) xmax = parts[3*i];
        if (parts[3*i+1] < ymin) ymin = parts[3*i+1];
        if (parts[3*i+1] > ymax) ymax = parts[3*i+1];
        if (parts[3*i+2] < zmin) zmin = parts[3*i+2];
        if (parts[3*i+2] > zmax) zmax = parts[3*i+2];
    }

    printf("Added guard points to total %ld points (should be %ld)\n", nvpall,
        nvpbuf + 6*(NGUARD+1)*(NGUARD+1));
    printf("x: %g,%g; y: %g,%g; z:%g,%g\n",xmin,xmax,ymin,ymax,zmin,zmax);

    /* allocate adjacencies using nvpall */
    adjs = (PARTADJ *)malloc(nvpall * sizeof(PARTADJ));
    if (adjs == NULL) {
        printf("Unable to allocate adjs\n");
        exit(0);
    }

    /* Do tesselation*/
    printf("File read.  Tessellating ...\n"); fflush(stdout);
    exitcode = delaunadj(parts, (int)nvp, (int)nvpbuf, (int)nvpall, &adjs);

    /* Now calculate volumes*/
    printf("Now finding volumes ...\n"); fflush(stdout);
    vols = (realT *)malloc(nvp * sizeof(realT));

    for (i = 0; i < nvp; i++) {
        /* Just the original particles
            Assign adjacency coordinate array*/
        /* Volumes */
        for (j = 0; j < adjs[i].nadj; j++)
            DL {
                deladjs[3*j + d] = parts[3*adjs[i].adj[j]+d] - parts[3*i+d];
                if (deladjs[3*j+d] < -0.5) deladjs[3*j+d]++;
                if (deladjs[3*j+d] > 0.5) deladjs[3*j+d]--;
            }

        exitcode = vorvol(deladjs, points, intpoints, adjs[i].nadj, &(vols[i]));
        vols[i] *= (realT)np;
        /*    if (i % 1000 == 0)
        printf("%d: %d, %g\n",i,adjs[i].nadj,vols[i]); */
    }

    /* Get the adjacencies back to their original values */
    for (i = 0; i < nvp; i++)
        for (j = 0; j < adjs[i].nadj; j++)
            adjs[i].adj[j] = orig[adjs[i].adj[j]];

    totalvol = 0.0;
    for (i = 0; i < nvp; i++) {
        totalvol += (realT)vols[i];
    }
    printf("Average volume = %g\n", totalvol / (realT)nvp);

    /* Now the output!
       First number of particles */
    sprintf(outfile, "part.%s.%02d.%02d.%02d", suffix, b[0], b[1], b[2]);

    printf("Output to %s\n\n",outfile);
    out = fopen(outfile, "w");
    if (out == NULL) {
        printf("Unable to open %s\n",outfile);
        exit(0);
    }
    fwrite(&np, sizeof(long), 1, out);
    fwrite(&nvp, sizeof(long), 1, out);
    printf("nvp = %ld\n", nvp);

    /* Tell us where the original particles were */
    fwrite(orig, sizeof(int), nvp, out);
    /* Volumes*/
    fwrite(vols, sizeof(realT), nvp, out);
    /* Adjacencies */
    for (i = 0; i < nvp; i++) {
        fwrite(&(adjs[i].nadj), sizeof(int), 1, out);
        if (adjs[i].nadj > 0)
            fwrite(adjs[i].adj, sizeof(int), adjs[i].nadj,out);
        else printf("0");
    }
    fclose(out);

    /* clean up memory */
    for (i = 0; i < nvp; i++) {
        if (adjs[i].nadj > 0)
            free(adjs[i].adj);
    }
    free(adjs);

    free(orig);
    free(parts);
    free(vols);
}
