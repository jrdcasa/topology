#ifndef __UNWRAP_H
#define __UNWRAP_H

#include <stdio.h>
#include <math.h>

#ifdef PARALLEL
  #include <omp.h>
  #define USED_OPENMP 1
#else
  #define USED_OPENMP 0
#endif

int iframe_nojump = 1;

static void c_unwrap(int nmols[], int neigh[], float coords[], float box[],
                     int nchains, int maxatomsch, int natoms, int maxneigh,
                     float coords_unwrap[])
{

    /*Unwrap molecules in a simulation box

    At the moment this only works with cubic boxes. Probalbly, although has not been checked, it works with prism boxes.

    nmols[nchains, maxatch] --> Int array. Index of the atoms for each chain. If the value is -1, the position
                               does not contain any atom. This is useful for trajectories with molecules (or chains)
                               of different number of atoms.
    neigh[natoms, 4]        -->  Int array. neigh[i*4+1] = [j, k, -1, -1] j and k are bonded to i.
    coords[natoms,3]        --> Float array. Wraped coordinates of the current step. Coordinates from the trajectory file.
    coords_0[natoms,3]      --> Float array. Unwraped coordinates from the previous step. This coordinates must not have
                                             the jumps corrected. Thus, these coordinates are the store in the trajectory file
                                             once the molecule has been unwrapped.
    n_image[natoms, 3]      --> Int array
    box[3]                  --> Float array. Dimensions of the box
    nchains                 --> Int. Number of chains (or molecules) in the system.
    maxatomsch              --> Int. The number of atoms of the molecule with maximum number of atoms. This is used for
                                     systems composed with molecules of different sizes.
    natoms                  --> Int. Total number of atoms in the system.
    maxneigh                --> Int. Maximum number of bonded (neighbours) atoms to a determinate atom
    coords_nojump           --> Float array. Unwraped coordinates with nojump correction applied. This is the return value

    */

    // Variables ==============================================================
    int i, ich, iat, jat, iatom, jatom, prev;
    int* isvisited = malloc(sizeof(int)*natoms);
    float d[3];

    // Initialize isvisited array
    for (int i=0; i<natoms; i++){
        isvisited[i] = 0;
    }

    // Unwrap the coordinates of the current step =============================
    #ifdef PARALLEL
    #pragma omp parallel for private(ich,iatom,jatom,prev, d, i, jat) \
    shared(coords,isvisited, nchains, maxatomsch, maxneigh, nmols, coords_unwrap)
    #endif
    for (ich=0; ich<nchains; ich++) {
        // First atom in the molecule acts as a pivot to unwrap the molecule
        iatom = nmols[ich*maxatomsch+0];
        for (i=0; i<3; i++){
            coords_unwrap[iatom*3+i] = coords[iatom*3+i];
        }
        isvisited[iatom] = 1;

        // Loop over the atoms of the current molecule
        for (iat=1; iat<maxatomsch; iat++) {
            iatom = nmols[ich*maxatomsch+iat];
            prev = -1;
            if (iatom == -1) continue;
            for (jat=0; jat<maxneigh; jat++) {
                jatom = neigh[iatom*maxneigh+jat];
                if (jatom == -1) continue;
                if (isvisited[jatom] == 1) prev = jatom;
            }
            isvisited[iatom] = 1;

            for (i=0; i<3; i++){
                d[i] = coords[iatom*3+i] - coords_unwrap[prev*3+i];
                d[i] = d[i] - box[i] * round(d[i]/box[i]);
                coords_unwrap[iatom*3+i] = coords_unwrap[prev*3+i] + d[i];
            }

        }
    }
    free(isvisited);
}

static void c_unwrap_nojump(int nmols[], int neigh[], float coords[], float coords_0[], int n_image[],
                     float box[], int nchains, int maxatomsch, int natoms, int maxneigh,
                     float coords_nojump[])
{

    /*

    At the moment this only works with cubic boxes. Probalbly, although has not been checked, it works with prism boxes.

    nmols[nchains, maxatch] --> Int array. Index of the atoms for each chain. If the value is -1, the position
                               does not contain any atom. This is useful for trajectories with molecules (or chains)
                               of different number of atoms.
    neigh[natoms, 4]        -->  Int array. neigh[i*4+1] = [j, k, -1, -1] j and k are bonded to i.
    coords[natoms,3]        --> Float array. Wraped coordinates of the current step. Coordinates from the trajectory file.
    coords_0[natoms,3]      --> Float array. Unwraped coordinates from the previous step. This coordinates must not have
                                             the jumps corrected. Thus, these coordinates are the store in the trajectory file
                                             once the molecule has been unwrapped.
    n_image[natoms, 3]      --> Int array
    box[3]                  --> Float array. Dimensions of the box
    nchains                 --> Int. Number of chains (or molecules) in the system.
    maxatomsch              --> Int. The number of atoms of the molecule with maximum number of atoms. This is used for
                                     systems composed with molecules of different sizes.
    natoms                  --> Int. Total number of atoms in the system.
    maxneigh                --> Int. Maximum number of bonded (neighbours) atoms to a determinate atom
    coords_nojump           --> Float array. Unwraped coordinates with nojump correction applied. This is the return value

    */


    // Variables ==============================================================
    int i, ich, iat, jat, iatom, jatom, prev, natch;
    int* isvisited = malloc(sizeof(int)*natoms);
    float* cog = malloc(sizeof(float)*nchains*3);
    float* cog_0 = malloc(sizeof(float)*nchains*3);
    float d[3], d_ich[3];

   // Initialize isvisited array ==============================================
    for (int i=0; i<natoms; i++){
        isvisited[i] = 0;
    }
    for (ich=0; ich<nchains; ich++) {
        for (i=0; i<3; i++) {
            cog[ich*3+i] = 0.0;
            cog_0[ich*3+i] = 0.0;
        }
    }

    // Center of geometry of the previous step ====================================
    // This value is used to determinate if there is a jump between steps.
    // The coordinates used must be unwrapped without the jumps correction applied.
    // In brief, these are the coordinates of the previous step from the trajectory file with the unwrapped
    // correction applied (aka apply PBC conditions)
    for (ich=0; ich<nchains; ich++) {

        iatom = nmols[ich*maxatomsch+0];

        for (i=0; i<3; i++){
            cog_0[ich*3+i] = coords_0[iatom*3+i];
        }
        isvisited[iatom] = 1;
        natch = 1;

        // Loop over the atoms of the current molecule
        for (iat=1; iat<maxatomsch; iat++) {
            iatom = nmols[ich*maxatomsch+iat];
            prev = -1;
            if (iatom == -1) continue;
            natch += 1;
            for (jat=0; jat<maxneigh; jat++) {
                jatom = neigh[iatom*maxneigh+jat];
                if (jatom == -1) continue;
                if (isvisited[jatom] == 1) prev = jatom;
            }
            isvisited[iatom] = 1;

            for (i=0; i<3; i++){
                cog_0[ich*3+i] += coords_0[iatom*3+i];
            }
        }
        for (i=0; i<3; i++){
            cog_0[ich*3+i] = cog_0[ich*3+i]/natch;
        }
    }

    // Initialize isvisited array =============================================
    for (int i=0; i<natoms; i++){
        isvisited[i] = 0;
    }

    // Unwrap the coordinates of the current step =============================
    #ifdef PARALLEL
    #pragma omp parallel for private(ich,iatom,jatom,prev, d, i, jat) \
    shared(coords,isvisited, nchains, maxatomsch, maxneigh, nmols, coords_nojump)
    #endif
    for (ich=0; ich<nchains; ich++) {
        // First atom in the molecule acts as a pivot to unwrap the molecule
        iatom = nmols[ich*maxatomsch+0];
        for (i=0; i<3; i++){
            coords_nojump[iatom*3+i] = coords[iatom*3+i];
        }
        isvisited[iatom] = 1;
        natch = 1;

        // Loop over the atoms of the current molecule
        for (iat=1; iat<maxatomsch; iat++) {
            iatom = nmols[ich*maxatomsch+iat];
            prev = -1;
            if (iatom == -1) continue;
            natch += 1;
            for (jat=0; jat<maxneigh; jat++) {
                jatom = neigh[iatom*maxneigh+jat];
                if (jatom == -1) continue;
                if (isvisited[jatom] == 1) prev = jatom;
            }
            isvisited[iatom] = 1;

            for (i=0; i<3; i++){
                d[i] = coords[iatom*3+i] - coords_nojump[prev*3+i];
                d[i] = d[i] - box[i] * round(d[i]/box[i]);
                coords_nojump[iatom*3+i] = coords_nojump[prev*3+i] + d[i];
                cog[ich*3+i] += coords_nojump[iatom*3+i];
            }
        }
        for (i=0; i<3; i++){
            cog[ich*3+i] = cog[ich*3+i]/natch;
        }

    }

    // Apply nojump correction ================================================
    // The trajectory needs to be sampled in such a way that the maximum move of a molecule in the system
    // sholud not be more that the half of the dimensions box.
    for (ich=0; ich<nchains; ich++) {
        for (i=0; i<3; i++){
             // Here cog is an unwrapped chain
             //      cog_0 is an unwrapped chain without jumps corrected
            d_ich[i] = cog[ich*3+i] - cog_0[ich*3+i];

            for (iat=0; iat<maxatomsch; iat++) {
                iatom = nmols[ich*maxatomsch+iat];
                if (iatom == -1) continue;
                  if (d_ich[i] >box[i]/2.0) {
                    n_image[iatom*3+i] -= 1;
                  } else if (d_ich[i] < -box[i]/2.0)  {
                    n_image[iatom*3+i] += 1;
                  }
                coords_nojump[iatom*3+i] = coords_nojump[iatom*3+i] + box[i]*n_image[iatom*3+i];
            }
        }
    }

    iframe_nojump += 1;
    free(isvisited);
    free(cog);
    free(cog_0);

}

static void c_center_of_geom(int natoms, float coord_current[], float cog[])
{

    /*

    Calculate the center of geometry of a set of coordinates. The coordinates must be unwrapped.

    natoms                  --> Int. Total number of atoms in the system.
    coords_current[natoms,3]--> Float array. Unwrapped coordinates of the current step. Coordinates from the trajectory file.
    cog[3]                  --> Coordinates if the center of geometry

    */

    int ndim = 3;
    int iat;
    int idim;


    for (idim=0; idim<ndim; idim++) cog[idim] = 0.0;

    // Sum coordinates over all atoms
    for (iat = 0; iat<natoms; iat++) {
        cog[0] += coord_current[iat*ndim+0];
        cog[1] += coord_current[iat*ndim+1];
        cog[2] += coord_current[iat*ndim+2];
    }

    // Divide by the number of atoms
    for (idim=0; idim<ndim; idim++) {
        cog[idim] = cog[idim] / natoms;
    }

}

static void c_center_of_mass(int natoms, float coord_current[], float mass[], float cog[])
{

    /*

    Calculate the center of mass of a set of coordinates. The coordinates must be unwrapped.

    natoms                  --> Int. Total number of atoms in the system.
    coords_current[natoms,3]--> Float array. Unwrapped coordinates of the current step. Coordinates from the trajectory file.
    mass[natoms]            --> Float array. Masses of the atoms
    cog[3]                  --> Coordinates if the center of geometry

    */

    int ndim = 3;
    int iat;
    int idim;
    float mass_total = 0.;


    for (idim=0; idim<ndim; idim++) cog[idim] = 0.0;

    // Sum coordinates over all atoms
    for (iat = 0; iat<natoms; iat++) {
        cog[0] += coord_current[iat*ndim+0]*mass[iat];
        cog[1] += coord_current[iat*ndim+1]*mass[iat];
        cog[2] += coord_current[iat*ndim+2]*mass[iat];
        mass_total += mass[iat];
    }

    // Divide by the number of atoms
    for (idim=0; idim<ndim; idim++) {
        cog[idim] = cog[idim] / mass_total;
    }

}

#endif