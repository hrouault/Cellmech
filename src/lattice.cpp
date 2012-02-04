/*
 * Copyright (C) 2004-2012 Hervé Rouault <herve.rouault@pasteur.fr>
 *
 * This file is part of Cellmech.
 * 
 * Cellmech is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Cellmech is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with Cellmech. If not, see <http://www.gnu.org/licenses/>.
 */

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_randist.h>

#include "opti.hpp"
#include "lattice.hpp"
#include "main.hpp"

unsigned int length_first_lattice=4;

Lattice::Lattice()
{

}

/* Generation of the first lattice (some pointers to cells remains void and shrink_datafirst() must be applied)
*/
Lattice::Lattice(unsigned int size)
{
   Vertex *pvertex;
   Cell *pcell;

   unsigned int last_first_lattice=size*(2*size-1);

   for (unsigned int i=0;i<last_first_lattice;i++){
      cells.push_back(NULL);
      verts.push_back(NULL);
   }

   /* Allocation of memory for the vertices
    * The code is not optimized at all but it is only used for initialisation and it is more readable
    */

   for (unsigned int row=2 ; row<2*size-1; row++){

      if (row%2){
         unsigned int column_lim=size + 1;
         for (unsigned int column=1 ; column<column_lim; column++){

            int vertex=vertex_from_coord(row, column);
            verts[vertex]=new Vertex();
         }
      }

      if ((row+1)%2){
         unsigned int column_lim=size - 1;
         for (unsigned int column=1 ; column<column_lim; column++){

            unsigned int vertex=vertex_from_coord(row, column);
            verts[vertex]=new Vertex();
         }

         for (unsigned int column=column_lim ; column < column_lim + 2; column++){
            unsigned int vertex=vertex_from_coord(row, column);
            verts[vertex]=NULL;
         }
      }
   }

   verts[0]=NULL;
   verts[1]=NULL;
   verts[size]=NULL;

   for (unsigned int column=2 ; column < size; column++){
      unsigned int vertex=column;
      verts[vertex]=new Vertex;
   }

   verts[size*(2*size-2)+1]=NULL;

   for (unsigned int column=2 ; column < size; column++){
      unsigned int vertex=size*(2*size-2)+column;
      verts[vertex]=new Vertex;
   }

   /* Allocation of memory for the cells */

   for (unsigned int row=1 ; row<2*size-1; row++){

      if (row%2){
         unsigned int column=1;
         unsigned int vertex=vertex_from_coord(row, column);
         cells[vertex]=NULL;

         unsigned int column_lim=size + 1;
         cells[size*(row-1)+1]=NULL;
         for (column=2 ; column<column_lim; column++){
            vertex=vertex_from_coord(row, column);
            if (column%2){
               cells[vertex]=new Cell;
            }
            else {
               cells[vertex]=NULL;
            }
         }
      }

      if ((row+1)%2){
         unsigned int column_lim=size+1;
         for (unsigned int column=1 ; column<column_lim; column++){
            unsigned int vertex=vertex_from_coord(row, column);

            if (column%2){
               cells[vertex]=new Cell;
            }
            else {
               cells[vertex]=NULL;
            }
         }
      }
   }
   unsigned int row=2*size-1;
   unsigned int column=1;
   unsigned int vertex=vertex_from_coord(row, column);
   cells[vertex]=NULL;

   unsigned int column_lim=size;
   for (unsigned int column=2 ; column<column_lim; column++){
      vertex=vertex_from_coord(row, column);
      if (column%2){
         cells[vertex]=new Cell;
      }
      else {
         cells[vertex]=NULL;
      }
   }


   cells[0]=new Cell;

   /* Effective construction of the lattice */

   for (row=2 ; row < 2*size -1 ; row++){

      if (row%2){

         column_lim=size;

         for (column=2 ; column<column_lim; column++){
            vertex=vertex_from_coord(row, column);
            pvertex=verts[vertex];
            if (pvertex){

               pvertex -> y=(row-1)*sqrt(3)/2;

               if (column%2){
                  pvertex -> x=-1.5+1.5*column;

                  pvertex -> pneighb[0]=verts[vertex+1];
                  pvertex -> pneighb[1]=verts[vertex+size-1];
                  pvertex -> pneighb[2]=verts[vertex-size-1];

                  pvertex -> pncell[0]=cells[vertex+size];	
                  pvertex -> pncell[1]=cells[vertex];
                  pvertex -> pncell[2]=cells[vertex-size];
               }

               if ((column+1)%2){
                  pvertex -> x=-2+1.5*column;

                  pvertex -> pneighb[0]=verts[vertex-1];
                  pvertex -> pneighb[1]=verts[vertex-size-1];
                  pvertex -> pneighb[2]=verts[vertex+size-1];

                  pvertex -> pncell[0]=cells[vertex-size-1];	
                  pvertex -> pncell[1]=cells[vertex+1];
                  pvertex -> pncell[2]=cells[vertex+size-1];
               }

            }
            else {
               locerror("coord_gen", "The vertex does not exist");
            }
         }
      }


      if ((row+1)%2){

         column_lim=size-1;

         for (column=1 ; column<column_lim; column++){
            vertex=vertex_from_coord(row, column);
            pvertex=verts[vertex];
            if (pvertex){

               pvertex -> y=(row-1)*sqrt(3)/2;

               if (column%2){
                  pvertex -> x=1.5*column;

                  pvertex -> pneighb[0]=verts[vertex+1];
                  pvertex -> pneighb[1]=verts[vertex+size+1];
                  pvertex -> pneighb[2]=verts[vertex-size+1];

                  pvertex -> pncell[0]=cells[vertex+size+2];	
                  pvertex -> pncell[1]=cells[vertex];
                  pvertex -> pncell[2]=cells[vertex-size+2];
               }

               if ((column+1)%2){
                  pvertex -> x=-0.5+1.5*column;

                  pvertex -> pneighb[0]=verts[vertex-1];
                  pvertex -> pneighb[1]=verts[vertex-size+1];
                  pvertex -> pneighb[2]=verts[vertex+size+1];

                  pvertex -> pncell[0]=cells[vertex-size+1];	
                  pvertex -> pncell[1]=cells[vertex+1];
                  pvertex -> pncell[2]=cells[vertex+size+1];
               }

            }
            else {
               locerror("coord_gen", "The vertex does not exist");
            }

         }
      }
   }

   column_lim=size-1;
   for (column=3 ; column<column_lim; column++){
      if (verts[column] && verts[size*(2*size-2)+column]){
         if (column%2){
            row=1;
            vertex=vertex_from_coord(row, column);
            pvertex=verts[vertex];
            pvertex -> x = -1.5 + 1.5*column;
            pvertex -> y = 0;

            pvertex -> pneighb[0]=verts[vertex+1];
            pvertex -> pneighb[1]=verts[vertex+size-1];
            pvertex -> pneighb[2]=verts[vertex-1];

            pvertex -> pncell[0]=cells[vertex+size];	
            pvertex -> pncell[1]=cells[vertex];
            pvertex -> pncell[2]=cells[0];


            row=2*size-1;
            vertex=vertex_from_coord(row, column);
            pvertex=verts[vertex];

            pvertex -> x = -1.5 + 1.5*column;
            pvertex -> y = (row-1)*sqrt(3)/2;

            pvertex -> pneighb[0]=verts[vertex+1];
            pvertex -> pneighb[1]=verts[vertex-1];
            pvertex -> pneighb[2]=verts[vertex-size-1];

            pvertex -> pncell[0]=cells[0];	
            pvertex -> pncell[1]=cells[vertex];
            pvertex -> pncell[2]=cells[vertex-size];
         }

         if ((column+1)%2){
            row=1;
            vertex=vertex_from_coord(row, column);
            pvertex=verts[vertex];

            pvertex -> x = -2 + 1.5*column;
            pvertex -> y = 0;

            pvertex -> pneighb[0]=verts[vertex-1];
            pvertex -> pneighb[1]=verts[vertex+1];
            pvertex -> pneighb[2]=verts[vertex+size-1];

            pvertex -> pncell[0]=cells[0];	
            pvertex -> pncell[1]=cells[vertex+1];
            pvertex -> pncell[2]=cells[vertex + size-1];


            row=2*size-1;
            vertex=vertex_from_coord(row, column);
            pvertex=verts[vertex];


            pvertex -> x = -2 + 1.5*column;
            pvertex -> y = (row-1)*sqrt(3)/2;

            pvertex -> pneighb[0]=verts[vertex-1];
            pvertex -> pneighb[1]=verts[vertex-size-1];
            pvertex -> pneighb[2]=verts[vertex+1];

            pvertex -> pncell[0]=cells[vertex-size-1];
            pvertex -> pncell[1]=cells[vertex+1];
            pvertex -> pncell[2]=cells[0];
         }
      }
      else {
         locerror("coord_gen", "The vertex does not exist");
      }
   }

   for (row=5; row< 2 * size - 4; row++){

      if (row%2){
         column=1;
         vertex=vertex_from_coord(row, column);
         pvertex=verts[vertex];

         if (pvertex){

            pvertex -> x = 0;
            pvertex -> y = (row-1)*sqrt(3)/2;

            pvertex -> pneighb[0]=verts[vertex+1];
            pvertex -> pneighb[1]=verts[vertex+2*size];
            pvertex -> pneighb[2]=verts[vertex-2*size];

            pvertex -> pncell[0]=cells[vertex+size];
            pvertex -> pncell[1]=cells[0];
            pvertex -> pncell[2]=cells[vertex-size];
         }
         else {
            locerror("coord_gen", "The vertex does not exist");
         }

         column=size;
         vertex=vertex_from_coord(row, column);
         pvertex=verts[vertex];

         if (pvertex){

            pvertex -> x = -2 + 1.5*size;
            pvertex -> y = (row-1)*sqrt(3)/2;

            pvertex -> pneighb[0]=verts[vertex-1];
            pvertex -> pneighb[1]=verts[vertex-2*size];
            pvertex -> pneighb[2]=verts[vertex+2*size];

            pvertex -> pncell[0]=cells[vertex-size-1];
            pvertex -> pncell[1]=cells[0];
            pvertex -> pncell[2]=cells[vertex+size-1];
         }
         else {
            locerror("coord_gen", "The vertex does not exist");
         }
      }
   }

   /* Corner generation */
   vertex=2;

   pvertex=verts[vertex];

   pvertex -> x = 1;
   pvertex -> y = 0;

   pvertex -> pneighb[0]=verts[2*size+1];
   pvertex -> pneighb[1]=verts[3];
   pvertex -> pneighb[2]=verts[size+1];

   pvertex -> pncell[0]=cells[0];
   pvertex -> pncell[1]=cells[3];
   pvertex -> pncell[2]=cells[size + 1];


   vertex=2*size+1;
   pvertex=verts[vertex];

   pvertex -> x = 0;
   pvertex -> y = sqrt(3);

   pvertex -> pneighb[0]=verts[2*size+2];
   pvertex -> pneighb[1]=verts[4*size+1];
   pvertex -> pneighb[2]=verts[2];

   pvertex -> pncell[0]=cells[3*size+1];
   pvertex -> pncell[1]=cells[0];
   pvertex -> pncell[2]=cells[size+1];


   vertex=size-1;
   pvertex=verts[vertex];

   pvertex -> x = -3 + 1.5*size;
   pvertex -> y = 0;

   pvertex -> pneighb[0]=verts[3*size];
   pvertex -> pneighb[1]=verts[2*size-2];
   pvertex -> pneighb[2]=verts[size-2];

   pvertex -> pncell[0]=cells[2*size-1];
   pvertex -> pncell[1]=cells[size-1];
   pvertex -> pncell[2]=cells[0];


   vertex=3*size;
   pvertex=verts[vertex];

   pvertex -> x = -2 + 1.5*size;
   pvertex -> y = sqrt(3);

   pvertex -> pneighb[0]=verts[3*size-1];
   pvertex -> pneighb[1]=verts[size-1];
   pvertex -> pneighb[2]=verts[5*size];

   pvertex -> pncell[0]=cells[2*size-1];
   pvertex -> pncell[1]=cells[0];
   pvertex -> pncell[2]=cells[4*size-1];


   vertex=size*(2*size-4)+1;
   pvertex=verts[vertex];

   pvertex -> x = 0;
   pvertex -> y = (2*size-4)*sqrt(3)/2;

   pvertex -> pneighb[0]=verts[size*(2*size-4)+2];
   pvertex -> pneighb[1]=verts[size*(2*size-2)+2];
   pvertex -> pneighb[2]=verts[size*(2*size-6)+1];

   pvertex -> pncell[0]=cells[size*(2*size-3)+1];
   pvertex -> pncell[1]=cells[0];
   pvertex -> pncell[2]=cells[size*(2*size-5)+1];


   vertex=size*(2*size-2)+2;
   pvertex=verts[vertex];

   pvertex -> x = 1;
   pvertex -> y = (2*size-2)*sqrt(3)/2;

   pvertex -> pneighb[0]=verts[size*(2*size-4)+1];
   pvertex -> pneighb[1]=verts[size*(2*size-3)+1];
   pvertex -> pneighb[2]=verts[size*(2*size-2)+3];

   pvertex -> pncell[0]=cells[size*(2*size-3)+1];
   pvertex -> pncell[1]=cells[size*(2*size-2)+3];
   pvertex -> pncell[2]=cells[0];


   vertex=size*(2*size-3);
   pvertex=verts[vertex];

   pvertex -> x = -2 + 1.5*size;
   pvertex -> y = (2*size-4)*sqrt(3)/2;

   pvertex -> pneighb[0]=verts[size*(2*size-3)-1];
   pvertex -> pneighb[1]=verts[size*(2*size-5)];
   pvertex -> pneighb[2]=verts[size*(2*size-1)-1];

   pvertex -> pncell[0]=cells[size*(2*size-4)-1];
   pvertex -> pncell[1]=cells[0];
   pvertex -> pncell[2]=cells[size*(2*size-2)-1];


   vertex=size*(2*size-1)-1;
   pvertex=verts[vertex];

   pvertex -> x = -3 + 1.5*size;
   pvertex -> y = (2*size-2)*sqrt(3)/2;

   pvertex -> pneighb[0]=verts[size*(2*size-3)];
   pvertex -> pneighb[1]=verts[size*(2*size-1)-2];
   pvertex -> pneighb[2]=verts[size*(2*size-2)-2];

   pvertex -> pncell[0]=cells[0];
   pvertex -> pncell[1]=cells[size*(2*size-1)-1];
   pvertex -> pncell[2]=cells[size*(2*size-2)-1];

   /* List of vertices in cells */

   /* Hexagons */

   for (row=2 ; row<2*size-1; row++){

      if (row%2){
         column_lim=size;
         for (column=3 ; column<column_lim; column++){
            if (column%2){

               int cell=vertex_from_coord(row, column);
               pcell=cells[cell];
               if (pcell){
                  pcell -> verts.push_back(verts[cell]);
                  pcell -> verts.push_back(verts[cell+size-1]);
                  pcell -> verts.push_back(verts[cell+size-2]);
                  pcell -> verts.push_back(verts[cell-1]);
                  pcell -> verts.push_back(verts[cell-size-2]);
                  pcell -> verts.push_back(verts[cell-size-1]);
               }
               else {
                  locerror("coord_gen", "The cell does not exist");
               }
            }

         }
      }

      if ((row+1)%2){
         column_lim=size - 2;
         for (column=3 ; column<column_lim; column++){

            if (column%2){

               int cell=vertex_from_coord(row, column);
               pcell=cells[cell];
               if (pcell){
                  pcell -> verts.push_back(verts[cell]);
                  pcell -> verts.push_back(verts[cell+size+1]);
                  pcell -> verts.push_back(verts[cell+size]);
                  pcell -> verts.push_back(verts[cell-1]);
                  pcell -> verts.push_back(verts[cell-size]);
                  pcell -> verts.push_back(verts[cell-size+1]);
               }
               else {
                  locerror("coord_gen", "The cell does not exist");
               }
            }
         }
      }
   }

   for (row=4 ; row<2*size-3; row++){

      if ((row+1)%2){
         column=1;
         int cell=vertex_from_coord(row, column);
         pcell=cells[cell];

         if (pcell){
            pcell -> verts.push_back(verts[cell]);
            pcell -> verts.push_back(verts[cell+size+1]);
            pcell -> verts.push_back(verts[cell+size]);
            pcell -> verts.push_back(verts[cell-size]);
            pcell -> verts.push_back(verts[cell-size+1]);
         }
         else {
            locerror("coord_gen", "The cell does not exist");
         }


         column=size-2;
         cell=vertex_from_coord(row, column);
         pcell=cells[cell+1];

         if (pcell){
            pcell -> verts.push_back(verts[cell]);
            pcell -> verts.push_back(verts[cell-size+1]);
            pcell -> verts.push_back(verts[cell-size+2]);
            pcell -> verts.push_back(verts[cell+size+2]);
            pcell -> verts.push_back(verts[cell+size+1]);
         }
         else {
            locerror("coord_gen", "The cell does not exist");
         }

      }
   }

   column_lim=size-1;
   for (column=2 ; column<column_lim; column++){
      if ((column+1)%2){
         row=1;
         int cell=vertex_from_coord(row, column);
         pcell=cells[cell+1];
         if (pcell){
            pcell -> verts.push_back(verts[cell]);
            pcell -> verts.push_back(verts[cell+1]);
            pcell -> verts.push_back(verts[cell+size]);
            pcell -> verts.push_back(verts[cell+size-1]);
         }
         else {
            locerror("coord_gen", "The cell+1 does not exist");
         }


         row=2*size-1;
         cell=vertex_from_coord(row, column);
         pcell=cells[cell+1];
         if (pcell){
            pcell -> verts.push_back(verts[cell]);
            pcell -> verts.push_back(verts[cell-size-1]);
            pcell -> verts.push_back(verts[cell-size]);
            pcell -> verts.push_back(verts[cell+1]);
         }
         else {
            locerror("coord_gen", "The cell+1 does not exist");
         }
      }
   }


   pcell=cells[size+1];
   if (pcell){
      pcell -> verts.push_back(verts[size+1]);
      pcell -> verts.push_back(verts[2*size+2]);
      pcell -> verts.push_back(verts[2*size+1]);
      pcell -> verts.push_back(verts[2]);
   }
   else {
      locerror("coord_gen", "The size+1 does not exist");
   }

   pcell=cells[2*size-1];
   if (pcell){
      pcell -> verts.push_back(verts[3*size]);
      pcell -> verts.push_back(verts[3*size-1]);
      pcell -> verts.push_back(verts[2*size-2]);
      pcell -> verts.push_back(verts[size-1]);
   }
   else {
      locerror("coord_gen", "The size+1 does not exist");
   }

   int cell=size*(2*size-3)+1;
   pcell=cells[cell];
   if (pcell){
      pcell -> verts.push_back(verts[cell]);
      pcell -> verts.push_back(verts[size*(2*size-2)+2]);
      pcell -> verts.push_back(verts[size*(2*size-4)+1]);
      pcell -> verts.push_back(verts[size*(2*size-4)+2]);
   }
   else {
      locerror("coord_gen", "The size+1 does not exist");
   }

   cell=size*(2*size-2)-1;
   pcell=cells[cell];
   if (pcell){
      pcell -> verts.push_back(verts[size*(2*size-1)-1]);
      pcell -> verts.push_back(verts[size*(2*size-2)-2]);
      pcell -> verts.push_back(verts[size*(2*size-3)-1]);
      pcell -> verts.push_back(verts[size*(2*size-3)]);
   }
   else {
      locerror("coord_gen", "The size+1 does not exist");
   }


   /* Initialisation of the outside cell */

   pcell=cells[0];
   pcell -> verts = allocverts(4*(size-2));


   for (row=3 ; row<2*size-2; row++){

      if (row%2){
         column=1;
         vertex=vertex_from_coord(row, column);
         pvertex=verts[vertex];
         pcell -> verts[(row-1)/2-1]=pvertex;

         column=size;
         vertex=vertex_from_coord(row, column);
         pvertex=verts[vertex];
         pcell -> verts[3*size-(row-1)/2-6]=pvertex;
      }
   }


   column_lim=size;
   for (column=2 ; column<column_lim; column++){
      row=2*size-1;
      vertex=vertex_from_coord(row, column);
      pvertex=verts[vertex];
      pcell -> verts[size+column-4]=pvertex;

      row=1;
      vertex=vertex_from_coord(row, column);
      pvertex=verts[vertex];
      pcell -> verts[4*size-column-7]=pvertex;
   }

   Cell *pcentercell=cells[vertex_from_coord(size, size/2)];
   if (!pcentercell){
      pcentercell=cells[vertex_from_coord(size, size/2) - 1];
      if (!pcentercell){
         locerror(" coord_gen ", " The central cell does not exist");
      }
   }


   // Shrink data
   vpvert vertstemp;
   vpcell cellstemp;
   /* Node list shrinking */

   unsigned int nb_vertex_tot=0;

   for (unsigned int ivertex=0 ; ivertex < last_first_lattice ; ivertex++){
      if (verts[ivertex]){
         nb_vertex_tot++;
      }
   }
   vertstemp=allocverts(nb_vertex_tot);

   int idefvert=0;

   for (unsigned int ivertex=0 ; ivertex < last_first_lattice ; ivertex++){
      if (verts[ivertex]){
         vertstemp[idefvert]=verts[ivertex];
         idefvert++;
      }
   }
   verts=vertstemp;

   /* Cell list shrinking */    
   unsigned int nb_cell_tot=0;

   for (unsigned int icell=0 ; icell<last_first_lattice ; icell++){
      if (cells[icell]){
         nb_cell_tot++;
      }
   }
   cellstemp=alloccells(nb_cell_tot);

   //    printf("Total number of cells : %i\n", nb_cell_tot);

   int idefcell=0;

   for (unsigned int icell=0 ; icell<last_first_lattice ; icell++){
      if (cells[icell]){
         cellstemp[idefcell]=cells[icell];
         idefcell++;
      }
   }
   cells=cellstemp;

   time=0;
   for (icell ic=cells.begin()+1;ic!=cells.end();ic++){
      Cell & c=**ic;
      c.area_intr=1.0;
      c.vol_intr=1.0;
      c.thick_intr=1.0;
      c.compart=0;
      c.clone_index=0;
      c.cyc_shape=args_info.cycle_arg;
      c.cyc_scale=1.0;
      c.t_next_div=gsl_ran_gamma(rng,c.cyc_shape,c.cyc_scale/c.cyc_shape);
      c.just_div=false;
      c.thick=c.thick_intr;
   }
   (*cells.begin())->clone_index=-1;
   (*cells.begin())->vol_intr=0;
   (*cells.begin())->thick=0;

   nbt1=0;
   nbt1closetodiv=0;
   dirborderdir=-1;

   exploded=false;
}


   void
Simulation::relax()
{
   int icheck=1;

   do {
      if (args_info.third_dim_flag){
         optimize3d();
      } else {
         optimize();
      }
      icheck=lattice.check_t1();

   } while(icheck);
}

void
Simulation::optimize3d()
{
   gsl_vector *startpos;
   gsl_multimin_function_fdf hamilton;
   const gsl_multimin_fdfminimizer_type *minitype;
   gsl_multimin_fdfminimizer *minimize;

   lattice.update_pa();
   lattice.update_centers();
   lattice.update_texture();

   int N=2*lattice.verts.size()+lattice.cells.size()-1;

   hamilton.f=&energy3d;
   hamilton.df=&denergy3d;
   hamilton.fdf=&enerdener3d;
   hamilton.n=N;
   hamilton.params=this;

   startpos = gsl_vector_alloc(N);
   unsigned int i=0;
   for (civert iv=lattice.verts.begin();iv!=lattice.verts.end();++iv){
      gsl_vector_set(startpos,2*i,(*iv)->x);
      gsl_vector_set(startpos,2*i+1,(*iv)->y);
      i++;
   }
   i=0;
   for (cicell iv=lattice.cells.begin()+1;iv!=lattice.cells.end();++iv){
      gsl_vector_set(startpos,2*lattice.verts.size()+i,(*iv)->thick);
      i++;
   }

   minitype=gsl_multimin_fdfminimizer_conjugate_fr;
   minimize=gsl_multimin_fdfminimizer_alloc(minitype, N);

   gsl_multimin_fdfminimizer_set(minimize, &hamilton, startpos, 0.001, 0.001);

   unsigned int iter=0;
   int status;
   do {
      iter++;
      status = gsl_multimin_fdfminimizer_iterate(minimize);

//      if (status){
//         cout << status << endl;
//         break;
//      }

      status = gsl_multimin_test_gradient(minimize->gradient, 1e-3);
   } while (status==GSL_CONTINUE && iter<500);

   i=0;
   for (ivert iv=lattice.verts.begin();iv!=lattice.verts.end();++iv){
      (*iv)->x=gsl_vector_get(minimize->x,2*i);
      (*iv)->y=gsl_vector_get(minimize->x,2*i+1);
      i++;
   }
   i=0;
   for (icell iv=lattice.cells.begin()+1;iv!=lattice.cells.end();++iv){
      (*iv)->thick=gsl_vector_get(minimize->x,2*lattice.verts.size()+i);
      i++;
   }

   gsl_multimin_fdfminimizer_free(minimize);
   gsl_vector_free(startpos);
}

void
Simulation::optimize()
{
   gsl_vector *startpos;
   gsl_multimin_function_fdf hamilton;
   const gsl_multimin_fdfminimizer_type *minitype;
   gsl_multimin_fdfminimizer *minimize;

   int N=2*lattice.verts.size();


   hamilton.f=&energy;
   hamilton.df=&denergy;
   hamilton.fdf=&enerdener;
   hamilton.n=N;
   hamilton.params=this;

   startpos = gsl_vector_alloc(N);
   unsigned int i=0;
   for (civert iv=lattice.verts.begin();iv!=lattice.verts.end();++iv){
      gsl_vector_set(startpos,2*i,(*iv)->x);
      gsl_vector_set(startpos,2*i+1,(*iv)->y);
      i++;
   }

   minitype=gsl_multimin_fdfminimizer_conjugate_fr;
   minimize=gsl_multimin_fdfminimizer_alloc(minitype, N);

   gsl_multimin_fdfminimizer_set(minimize, &hamilton, startpos, 0.001, 0.001);

   unsigned int iter=0;
   int status;
   do {
      iter++;
      status = gsl_multimin_fdfminimizer_iterate(minimize);

//      if (status){
//         break;
//      }

      status = gsl_multimin_test_gradient(minimize->gradient, 1e-3);
   } while (status==GSL_CONTINUE && iter<500);

   i=0;
   for (ivert iv=lattice.verts.begin();iv!=lattice.verts.end();++iv){
      (*iv)->x=gsl_vector_get(minimize->x,2*i);
      (*iv)->y=gsl_vector_get(minimize->x,2*i+1);
      i++;
   }

   gsl_multimin_fdfminimizer_free(minimize);
   gsl_vector_free(startpos);
}

int
Lattice::check_t1()
{
   update_pa();
   double area=0;
   for (icell ic=cells.begin()+1;ic!=cells.end();ic++){
      area+=(*ic)->area;
   }
   area/=(double)(cells.size()-1);
   for (ivert iv=verts.begin();iv!=verts.end();iv++){
      for (unsigned j=0;j<3;j++){
         Vertex * pvertn=(*iv)->pneighb[j];
         if (dist(*iv,pvertn)<0.01*sqrt(area)){
            t1transform(*iv,j);
            nbt1++;
            return 1;
         }
      }
   }

   return 0;
}

/* Implementation of the t1 process
 * It has to be noticed that the choice of the vertex c is not arbitrary 
 * and can induce a rotation on the lattice. This effect should be negligible though (it is in practice)
 */
void
Lattice::t1transform(Vertex * pvert, unsigned int neighbor)
{
   Vertex * pvertexb=pvert->pneighb[neighbor];

   Cell * pcellalpha=pvert->pncell[neighbor];
   Cell * pcellbeta=pvert->pncell[(neighbor+2)%3];
   unsigned int indexbneighbc=pvertexb->indexcell(pcellalpha);
   Vertex * pvertexc=pvertexb->pneighb[indexbneighbc];
   Vertex * pvertexd=pvert->pneighb[(neighbor+2)%3];

   Cell * pcellgamma=pvert->pncell[(neighbor+1)%3];
   Cell * pcelldelta=pvertexb->pncell[(indexbneighbc+2)%3];

   if (pcellgamma->just_div || pcelldelta->just_div){
      nbt1closetodiv++;
   }

   unsigned int indexcneighbb=pvertexc->indexcell(pcelldelta);
   unsigned int indexdneighba=pvertexd->indexcell(pcellgamma);

   Vertex * pvertexe=pvert->pneighb[(neighbor+1)%3];
   Vertex * pvertexf=pvertexb->pneighb[(indexbneighbc+2)%3];


   /* Changes */

   pvertexd->pneighb[indexdneighba]=pvertexb;
   pvertexc->pneighb[indexcneighbb]=pvert;

   pcellgamma->insvertex(pvertexb, pvert);
   pcelldelta->insvertex(pvert, pvertexb);

   delvertex(*pcellalpha, pvertexb);
   delvertex(*pcellbeta, pvert);

   pvert->pneighb[0]=pvertexc;
   pvert->pneighb[1]=pvertexe;
   pvert->pneighb[2]=pvertexb;
   pvert->pncell[0]=pcellalpha;
   pvert->pncell[1]=pcellgamma;
   pvert->pncell[2]=pcelldelta;

   pvertexb->pneighb[0]=pvertexd;
   pvertexb->pneighb[1]=pvertexf;
   pvertexb->pneighb[2]=pvert;
   pvertexb->pncell[0]=pcellbeta;
   pvertexb->pncell[1]=pcelldelta;
   pvertexb->pncell[2]=pcellgamma;
}

   void
Lattice::delvertex(Cell & cell, Vertex * pvert)
{
   unsigned int ivert;
   if ((ivert=cell.indexvertex(pvert)) < 0){
      printf("In delvertex ivert < 0\n");
      exit(1);
   };

   unsigned int nbvertices=cell.verts.size();
   vpvert vertstmp;

   for (unsigned int n=0;n<nbvertices-1;n++){
      vertstmp.push_back(cell.verts[(ivert+n+1)%nbvertices]);
   }

   cell.verts = vertstmp;
}

/* Gives the indice of a vertex from the cartesian coordinates on the first lattice
 */
   int
vertex_from_coord(unsigned int row, unsigned int column)
{
   if (row>2*length_first_lattice-1 || column>length_first_lattice){
      locerror("vertex_from_coord", "coordinates went over the limits");
   }
   return length_first_lattice*(row-1)+column;
}

void
Lattice::division(Cell & c, unsigned int iv1, unsigned int iv2)
{
   c.comp_area();
   vd ccm=c.center_mass();
   vd clat0=center();
   Cell * pc0=closest_center_onborder();
   pc0->comp_area();
   vd cfcob0=pc0->center_mass();
   double dfcob0=dist(clat0,cfcob0);

   double dc0=dist(clat0,ccm);

   fprintf(fdivpos, "%f %i %f %f\n", time, cells.size()-1, dc0, dfcob0);
   fflush(fdivpos);
   
      

   Vertex & v=*c.verts[iv1];

   Vertex * pnvert1=new Vertex();
   Vertex * pnvert2=new Vertex();
   verts.push_back(pnvert1);
   verts.push_back(pnvert2);

   Cell * pnc=new Cell();
   cells.push_back(pnc);
   Cell & nc=*pnc;
   
   nc.thick=c.thick;

   unsigned int nbvertices=c.verts.size();

   Vertex * pvertex1=c.verts[(iv1+1)%nbvertices];
   Vertex * pvertex2=c.verts[iv2];
   Vertex * pvertex3=c.verts[(iv2+1)%nbvertices];

   unsigned int icell0=v.indexcell(&c);
   unsigned int icell1=pvertex1->indexcell(&c);
   unsigned int icell2=pvertex2->indexcell(&c);
   unsigned int icell3=pvertex3->indexcell(&c);

   Cell * pncell0m1=v.pncell[(icell0+2)%3];
   Cell * pncell2m1=pvertex2->pncell[(icell2+2)%3];

   pnvert1->x=(v.x + pvertex1->x)/2;
   pnvert1->y=(v.y + pvertex1->y)/2;
   pnvert1->pneighb[0]=&v;
   pnvert1->pneighb[1]=pvertex1;
   pnvert1->pneighb[2]=pnvert2;
   pnvert1->pncell[0]=pncell0m1;
   pnvert1->pncell[1]=pnc;
   pnvert1->pncell[2]=&c;

   pnvert2->x=(pvertex2->x + pvertex3->x)/2;
   pnvert2->y=(pvertex2->y + pvertex3->y)/2;
   pnvert2->pneighb[0]=pvertex2;
   pnvert2->pneighb[1]=pvertex3;
   pnvert2->pneighb[2]=pnvert1;
   pnvert2->pncell[0]=pncell2m1;
   pnvert2->pncell[1]=&c;
   pnvert2->pncell[2]=pnc;


   v.pneighb[icell0]=pnvert1;

   pvertex1->pneighb[(icell1+1)%3]=pnvert1;

   pvertex2->pneighb[icell2]=pnvert2;

   pvertex3->pneighb[(icell3+1)%3]=pnvert2;

   for (unsigned int n=iv1+1 ; n < iv2 +1 ; n++){
      Vertex * pvert=c.verts[n%nbvertices];
      unsigned int icelltemp=pvert->indexcell(&c);
      pvert->pncell[icelltemp]=pnc;
   }

   nc.verts.push_back(pnvert1);
   for (unsigned int n=iv1+1 ; n < iv2 +1 ; n++){
      nc.verts.push_back(c.verts[n%nbvertices]);
   }
   nc.verts.push_back(pnvert2);

   vpvert pvertexlisttemp;

   pvertexlisttemp.push_back(pnvert2);
   for (unsigned int n=iv2+1 ; n < iv1+nbvertices+1 ; n++){
      pvertexlisttemp.push_back(c.verts[n%nbvertices]);
   }
   pvertexlisttemp.push_back(pnvert1);

   c.verts=pvertexlisttemp;

   pncell0m1->insvertex(pnvert1 ,&v);

   pncell2m1->insvertex(pnvert2 ,pvertex2);

   nc.area_intr=c.area_intr=1.0;
   nc.vol_intr=c.vol_intr=1.0;
   nc.thick_intr=c.thick_intr=1.0;
   nc.clone_index=c.clone_index;
   nc.compart=c.compart;
   nc.cyc_scale=c.cyc_scale;
   nc.cyc_shape=c.cyc_shape;
   nc.thick=nc.thick_intr;
   int nbtot=cells.size();
   c.comp_area();
   nc.comp_area();
   vd c1=c.center_mass();
   vd c2=nc.center_mass();
   vd clat=center();
   Cell * pc=closest_center_onborder();
   pc->comp_area();
//   vd cfcob=pc->center_mass();
//   double dfcob=dist(clat,cfcob);
//
//   double dc1=dist(clat,c1);
////      double gamc=1/args_info.clone_cycle_arg;
////      double gam=gamc-(gamc-1)*dc1*dc1/dfcob/dfco 
//      double sigmagam=dfcob*dfcob/20;
////      double gam=exp(-dc1*dc1/2/sigmagam)1-9*dc1*dc1/dfcob/dfcob;
//
////      double gam=exp(-dc1*dc1/2/sigmagam);
//      double gam=1-3.0*dc1/dfcob;
//      if (gam < 0) gam=0;
//      c.cyc_scale=1/gam;
//   double dc2=dist(clat,c2);
////      double gamc=1/args_info.clone_cycle_arg;
////      double gam=gamc-(gamc-1)*dc2*dc2/dfcob/dfcob;
//
////      gam=exp(-dc2*dc2/2/sigmagam);
//
//      gam=1-3.0*dc2/dfcob;
//      if (gam < 0) gam=0;
//      nc.cyc_scale=1/gam;
//
//      cout << c.cell_center[0] << " " << c.cell_center[1] << " " << dfcob << endl;
   
   c.cyc_scale=1.0;
   nc.cyc_scale=1.0;
   c.t_next_div=time+gsl_ran_gamma(rng,c.cyc_shape,c.cyc_scale/c.cyc_shape);
   nc.t_next_div=time+gsl_ran_gamma(rng,nc.cyc_shape,nc.cyc_scale/nc.cyc_shape);

   for (icell ic=cells.begin();ic!=cells.end();++ic){
      (*ic)->just_div=false;
   }
   nc.just_div=true;
   c.just_div=true;
}

void
Lattice::division(Cell & c, double theta)
{
   vint listv;
   vd ct=c.center();
   vd directv;
   directv.push_back(cos(theta));
   directv.push_back(sin(theta));
   for (unsigned int i=0;i<c.verts.size();i++){
      Vertex & v1=*c.verts[i];
      Vertex & v2=*c.verts[(i+1)%c.verts.size()];
      if (((v1.x-ct[0])*directv[1]-(v1.y-ct[1])*directv[0])*((v2.x-ct[0])*directv[1]-(v2.y-ct[1])*directv[0])<0){
         listv.push_back(i);
      }
   }
   division(c,listv[0],listv[1]);
}

   void
Lattice::division_alongbound(Cell & c)
{

   if (dirborderdir>-1){
      for (ivert iv=c.verts.begin();iv!=c.verts.end();++iv){
         Vertex & v=**iv;
         unsigned int ind=v.indexcell(&c);
         Cell & c2=*v.pncell[(ind+2)%3];
         Vertex & v2=*v.pneighb[ind];
         if (c2.compart!=-1 && c2.compart!=c.compart){
            double x=v.x-v2.x;
            double y=v.y-v2.y;
            double theta=atan(y/x);
            if (dirborderdir==0){
               division(c,theta+M_PI/2);
            } else {
               division(c,theta);
            }
            return;
         }
      }
   }

   division(c);
}

void
Lattice::division_polarclone(Cell & c)
{
   division(c);
}

void
Lattice::division(Cell & c)
{

   vd ccent=c.center_mass();
   vd clat=center();
   Cell * pc=closest_center_onborder();
   pc->comp_area();
   vd cfcob=pc->center_mass();
   double dfcob=dist(clat,cfcob);
   double dc1=dist(clat,ccent);


   double sig2 = dfcob*dfcob/20;
   double k = 0.4*exp(-dc1*dc1/2/sig2);

   double disp=von_mises_dist(k)/2.0;

   double angle=0;
   if (fabs(ccent[0]-clat[0])<0.1){
      angle=M_PI/2.0;
   } else {
      angle=atan((ccent[1]-clat[1])/(ccent[0]-clat[0]));
   }
   angle+=M_PI/2.0+disp;

   division(c,angle);
}

void
Lattice::division_randdir(Cell & c)
{
   double angle=2*M_PI*gsl_rng_uniform(rng);

   division(c,angle);
}

Cell *
Lattice::closest_center()
{
   Cell * cur_cell=cells[1];

   vd clat=center();
   double cx=clat[0];
   double cy=clat[1];
   double dist=1e8;
   for (icell ic=cells.begin()+1;ic!=cells.end();ic++){
      Cell & c=**ic;

      vd centerc=c.center();
      double distance=(centerc[0]-cx)*(centerc[0]-cx)+(centerc[1]-cy)*(centerc[1]-cy);
      if (distance<dist){
         dist=distance;
         cur_cell=&c;
      }
   }
   if (dist>1e7){
      printf("dist > 1e7\n");
      exploded=true;
   }
   return cur_cell;
}

Cell *
Lattice::closest_center_onborder()
{
   Vertex * cur_vert=verts[0];

   vd clat=center();
   double cx=clat[0];
   double cy=clat[1];
   double dist=1e8;
   for (ivert iv=verts.begin();iv!=verts.end();iv++){
      Vertex & v=**iv;

      double distance=(v.x-cx)*(v.x-cx)+(v.y-cy)*(v.y-cy);
      if (belongborder(v) && distance<dist){
         dist=distance;
         cur_vert=&v;
      }
   }
   if (dist>1e7){
      printf("dist > 1e7\n");
      exploded=true;
   }
   if (cur_vert->pncell[0]!=cells[0]) return cur_vert->pncell[0];
   return cur_vert->pncell[1];
}

/* Differentiate compartments
 */
void
Lattice::diffcompart()
{
   vd clat=center();
   (*cells.begin())->compart=-1;
   for (icell ic=cells.begin()+1;ic!=cells.end();ic++){
      vd cc=(*ic)->center();
      if (cc[0]<clat[0]){
         (*ic)->compart=0;
      } else {
         (*ic)->compart=1;
      }
   }
}

/* Length of the boundary between compartments
 */
double
Lattice::boundarylength()
{
   double length=0;
   for (ivert iv=verts.begin();iv!=verts.end();++iv){
      for (unsigned int i=0;i<3;i++){
         int comp1=(*iv)->pncell[i]->compart;
         int comp2=(*iv)->pncell[(i+2)%3]->compart;
         if (comp1!=-1 && comp2!=-1 && comp1!=comp2){
            length+=dist(*iv,(*iv)->pneighb[i]);
         }
      }
   }
   length/=2;
   return length;
}

bool
Lattice::belongborder(Vertex & v)
{
   for (icell ic=v.pncell.begin();ic!=v.pncell.end();ic++){
      if (*ic==cells[0]) return true;
   }
   return false;
}

   /* Select a clone with nb_cells cells
 * The markers size, div_rate and clone_index are set to 1
 */
void
Lattice::diffclone(unsigned int nb_cells)
{
   double radius=0.05; /* Increments the radius until it gets nb_cells around lattice_center */

   /* Center cell method */
   vd latc=center();

   unsigned int nb_clone_rad=0;
   while (nb_clone_rad<nb_cells){
      nb_clone_rad=0;
      for (icell ic=cells.begin()+1;ic!=cells.end();ic++){
         vd centerc=(*ic)->center();

         double dx=latc[0]-centerc[0];
         double dy=latc[1]-centerc[1];
         if (dx*dx+dy*dy<radius){
            (*ic)->clone_index=1;
            (*ic)->cyc_scale*=args_info.clone_cycle_arg;
            nb_clone_rad++;
         }
      }
      radius += 0.05;
   }
}

/* Select a clone composed of one cell in the center of the lattice
 *    
 */  
void
Lattice::diffclone()
{
   Cell * c=closest_center();

   c->clone_index=1;
   c->cyc_scale*=args_info.clone_cycle_arg;
}


/* Computes the center of the lattice
 */
vd
Lattice::center()
{
   Cell & c=*cells[1];

   vd centerc=c.center();
   double left,right,bottom,top;
   left=right=centerc[0];
   bottom=top=centerc[1];

   for (icell ic=cells.begin()+2;ic!=cells.end();ic++){
      centerc=(*ic)->center();
      if (centerc[0] < left) left=centerc[0];
      if (centerc[0] > right) right=centerc[0];
      if (centerc[1] < bottom) bottom=centerc[1];
      if (centerc[1] > top) top=centerc[1];
   }

   vd clat;
   clat.push_back(0.5*(left+right));
   clat.push_back(0.5*(bottom+top));

   return clat;
}


vpvert
allocverts(int nbverts)
{
   vpvert res;
   for (int i=0;i<nbverts;i++){
      res.push_back(new Vertex);
   }
   return res;
}

vpcell
alloccells(int nbcells)
{
   vpcell res;
   for (int i=0;i<nbcells;i++){
      res.push_back(new Cell);
   }
   return res;
}

   void
locerror(string function, string message)
{
   cerr << "Fatal error : function " << function << ", " << message << endl;
   exit(1);
}

/* Svg file generation : takes a snapshop of the lattice
 */
   void
Lattice::latticeprint(string filename, unsigned int type)
{
   /* File initialization */
   ofstream file(filename.c_str());
   string ext=filename.substr(filename.length()-3,3);

   update_pa();
   update_centers();

   if (ext=="svg"){
      file << "<?xml version=\"1.0\" encoding=\"utf-8\" standalone=\"no\"?>\n\
         <!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"\n\
         \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n\
         <svg\n   xmlns=\"http://www.w3.org/2000/svg\"\n\
         version=\"1.1\"\n\
         x=\"0mm\"\n\
         y=\"0mm\"\n\
         width=\"250mm\"\n\
         height=\"250mm\">\n\
         <desc>Growing lattice\n\
         </desc>\n\
         <defs>\n\
         </defs>\n\
         <g \n\
         style=\"fill:white;stroke:black;stroke-width:0;stroke-linecap:round;stroke-linejoin:round\">\n";
   } else if (ext=="x3d"){
      file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
      file << "<!DOCTYPE X3D PUBLIC \"ISO//Web3D//DTD X3D 3.0//EN\" \"http://www.web3d.org/specifications/x3d-3.0.dtd\">\n";
      file << "\n";
      file << "<X3D profile='Interchange' version='3.0' xmlns:xsd='http://www.w3.org/2001/XMLSchema-instance' xsd:noNamespaceSchemaLocation=' http://www.web3d.org/specifications/x3d-3.0.xsd '>\n";
      file << "   <head>\n";
      file << "   <meta name='title' content='PointSet.x3d'/>\n";
      file << "   <meta name='description' content='X3D specification example: PointSet demonstration.'/>\n";
      file << "   <meta name='editors' content='Hervé Rouault'/>\n";
      file << "   <meta name='created' content='April 2010'/>\n";
      file << "   </head>\n";
      file << "   <Scene>\n";
   } else {
      cerr << "EROR : filetype not recognized : " << ext << endl;
   }

   double scale=0;
   for (icell ic=cells.begin()+1;ic!=cells.end();ic++){
      double pres=0;
      if (args_info.third_dim_flag){
         pres=fabs((*ic)->thick-(*ic)->thick_intr);
      } else {
         pres=fabs((*ic)->pressure());
      }
      if (scale<pres) scale=pres;
   }

   unsigned int i=1;
   for (icell ic=cells.begin()+1;ic!=cells.end();ic++){
      stringstream sticker;
      sticker << i;
      (*ic)->latticeprint(file,type,sticker.str(),scale,ext);
      i++;
   }

   double sca=20;
   double transx=300;
   double transy=300;
   vd clat=center();
   file << "         <circle cx=\"" << clat[0]*sca+transx << "\" cy=\"" << clat[1]*sca+transy << "\" r=\"5.0\" fill=\"red\" />\n";
   Cell * pcbound=closest_center_onborder();
   vd cbound=pcbound->center();
   file << "         <circle cx=\"" << cbound[0]*sca+transx << "\" cy=\"" << cbound[1]*sca+transy << "\" r=\"5.0\" fill=\"red\" />\n";

   printcompartline(file);
//   printgradient(file);
   if (ext=="svg"){
      file << "   </g>\n</svg>";
   } else if (ext=="x3d"){
      file << "   </Scene>\n";
      file << "   </X3D>\n";
   }
   file.close();
}

   void
Lattice::printcompartline(ofstream & file)
{
   double scale=20;
   double transx=300;
   double transy=300;
   for (ivert iv=verts.begin();iv!=verts.end();++iv){
      for (unsigned int i=0;i<3;i++){
         int comp1=(*iv)->pncell[i]->compart;
         int comp2=(*iv)->pncell[(i+2)%3]->compart;
         if (comp1!=-1 && comp2!=-1 && comp1!=comp2){
            file << "         <line x1=\"" << (*iv)->x*scale+transx << "\" y1=\"" << (*iv)->y*scale+transy << "\" x2=\"" << (*iv)->pneighb[i]->x*scale+transx << "\" y2=\"" << (*iv)->pneighb[i]->y*scale+transy <<"\" stroke-width=\"3\" stroke=\"blue\" />\n";
         }
      }
   }
}

void
Lattice::update_pa()
{
   for (icell ic=cells.begin();ic!=cells.end();ic++){
      (*ic)->comp_area();
      (*ic)->comp_peri();
   }
}

void
Lattice::update_centers()
{
   for (icell ic=cells.begin();ic!=cells.end();ic++){
      (*ic)->cell_center=(*ic)->center_mass();
   }
}

void
Lattice::update_texture()
{
   for (icell ic=cells.begin()+1;ic!=cells.end();ic++){
      vd text=(*ic)->texture();
      (*ic)->ellipticity=text[3];
      (*ic)->angle=text[4];
   }
}

void
Lattice::printneighbors(ofstream & file)
{
   unsigned int i=1;
   for (icell ic=cells.begin()+1;ic!=cells.end();++ic){
      Cell & c=**ic;
      c.index=i;
      i++;
   }
   for (icell ic=cells.begin()+1;ic!=cells.end();++ic){
      Cell & c=**ic;
      if (c.clone_index>0){
         int test_lone=1;
         vpcell cs=c.neighbs();
         for (icell icc=cs.begin();icc!=cs.end();++icc){
            if ((*icc)->clone_index>0){
               file << c.index << "   " << (*icc)->index << endl;
               test_lone=0;
            }
         }
         if (test_lone){
            file << c.index << "   " << c.index << endl;
         }
      }
   }
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  reinittime
 *  Description:  After initial growth, update the phase of the cell cycle to meet 
 *  growth rate pattern
 * =====================================================================================
 */
void
Lattice::reinittime()
{
   // The phase is drawn randomly (see explanatory text for further detail)
   time=0;
   for (icell ic=cells.begin()+1;ic!=cells.end();ic++){
      Cell & c=**ic;

      c.comp_area();
      vd cc=c.center_mass();
      vd clat=center();
      Cell * pc=closest_center_onborder();
      pc->comp_area();
//      vd cfcob=pc->center_mass();
//      double dfcob=dist(clat,cfcob);
////
////      double dc1=dist(clat,cc);
////      if (dc1<dfcob){
////         c.cyc_scale=(1-args_info.clone_cycle_arg)*dc1/dfcob+args_info.clone_cycle_arg;
////      } else {
////         c.cyc_scale=1.0;
////      }
//   double dc1=dist(clat,cc);
////      double gamc=1/args_info.clone_cycle_arg;
////      double gam=gamc-(gamc-1)*dc1*dc1/dfcob/dfco 
//      double sigmagam=dfcob*dfcob/20;
////      double gam=exp(-dc1*dc1/2/sigmagam)1-9*dc1*dc1/dfcob/dfcob;
//
//
////      double gam=exp(-dc1*dc1/2/sigmagam);
//      double gam=1-3.0*dc1/dfcob;
//      if (gam < 0) gam=0;
//      
//      c.cyc_scale=1/gam;
      c.cyc_scale=1.0;
      double rnb=gsl_rng_uniform(rng);
      c.t_next_div=c.cyc_scale*log(rnb+1)/log(2);
   }
}
         
void
Simulation::div_next()
{
   double next_time=1e8;
   Cell * nc=NULL;
   for (icell ic=lattice.cells.begin()+1;ic!=lattice.cells.end();++ic){
      if (next_time>(*ic)->t_next_div){
         next_time=(*ic)->t_next_div;
         nc=*ic;
      }
   }
   double old_time=lattice.time;
   lattice.time=next_time;

   unsigned int nbtot=lattice.cells.size();
   for (icell ic=lattice.cells.begin()+1;ic!=lattice.cells.end();++ic){
      (*ic)->area_intr+=(lattice.time-old_time)/(*ic)->cyc_scale;
      (*ic)->vol_intr+=(lattice.time-old_time)/(*ic)->cyc_scale;
      if (lattice.time>(*ic)->t_next_div-0.0625 && old_time<(*ic)->t_next_div-0.0625){
         (*ic)->area_intr+=1.0;
         (*ic)->vol_intr*=1.36;
         relax();
         (*ic)->area_intr+=1.0;
         (*ic)->vol_intr*=1.36;
         relax();
         (*ic)->area_intr+=1.0;
         (*ic)->vol_intr*=1.36;
         relax();
      }
   }

   lattice.division_polarclone(*nc);
}

/* 
 * void
 * Lattice::clone_outline(ofstream & file)
 * {
 *    Dual *pc0=web_dual[0];
 *    Node *nvert;
 *    Node *pvert;
 *    double tot_len=perimeter(pc0);
 *    double s=0;
 *    double d,cost;
 *    int iv;
 * 
 *    Node **vertl=pc0->vertexlist;
 *    for (iv=1;iv<pc0->nb_vertices;iv++){
 *       pvert=vertl[iv-1];
 *       nvert=vertl[iv];
 *       d=dist(pvert,nvert);
 *       cost=(nvert->x-pvert->x)/d;
 *       printf("%f %f\n",s/tot_len,cost);
 *       s+=d;
 *    }
 *    printf("\n");
 * }
 */

unsigned int
Lattice::clone_size()
{
   unsigned int nb=0;
   for (icell ic=cells.begin()+1;ic!=cells.end();++ic){
      if ((*ic)->clone_index>0) nb++;
   }
   return nb;
}

vd
Lattice::centermass_clone()
{
   vd cmass;
   cmass.push_back(0);
   cmass.push_back(0);
   unsigned int i=0;
   for (icell ic=cells.begin()+1;ic!=cells.end();ic++){
      if ((*ic)->clone_index>0){
         vd centerc=(*ic)->center();
         cmass[0]+=centerc[0];
         cmass[1]+=centerc[1];
         i++;
      }
   }
   cmass[0]/=i;
   cmass[1]/=i;

   return cmass;
}

double
Lattice::anisotropy_clone()
{
   double mxx=0;
   double mxy=0;
   double mass=0;
   vd cmass=centermass_clone();
   for (icell ic=cells.begin()+1;ic!=cells.end();ic++){
      if ((*ic)->clone_index>0){
         vd centerc=(*ic)->center();
         double rx=centerc[0]-cmass[0];
         double ry=centerc[1]-cmass[1];
         mass+=rx*rx+ry*ry;
         mxx+=rx*rx-ry*ry;
         mxy+=2*rx*ry;
      }
   }
   return sqrt(mxx*mxx+mxy*mxy)/mass;
}

double
Lattice::meanpressure_clone()
{
   double areaclone=0;
   double pressure=0;

   for (icell ic=cells.begin()+1;ic!=cells.end();ic++){
      if ((*ic)->clone_index>0){
         (*ic)->comp_area();
         double area=(*ic)->area;
         areaclone+=area;
         pressure+=area*(*ic)->pressure();
      }
   }
   pressure/=areaclone;
   return pressure;
}

   void
Lattice::adddivdir()
{
   vd divcenter;
   for (icell ic=cells.begin()+1;ic!=cells.end();++ic){
      if ((*ic)->just_div && (*ic)->clone_index>0){
         vd cent=(*ic)->center();
         divcenter.push_back(cent[0]);
         divcenter.push_back(cent[1]);
      }
   }
   if (divcenter.size()){
      division_dirs.push_back(divcenter);
   }
}

vvd
Lattice::angles_aroundclone()
{
   vd cclone=centermass_clone();
   vvd res;
   for (icell ic=cells.begin()+1;ic!=cells.end();++ic){
      vd line;
      vd centerc=(*ic)->center();
      double dist2=(centerc[0]-cclone[0])*(centerc[0]-cclone[0])+(centerc[1]-cclone[1])*(centerc[1]-cclone[1]);
      double dist=sqrt(dist2);
      line.push_back(dist);
      double ux=(centerc[0]-cclone[0])/dist;
      double uy=(centerc[1]-cclone[1])/dist;
//      vd cshape=(*ic)->shape();
//      double angle=cshape[1];
      vd cshape=(*ic)->texture();
      double angle=cshape[4];
      double vx=cos(angle);
      double vy=sin(angle);
      line.push_back(atan((ux*vy-uy*vx)/(ux*vx+uy*vy)));
      res.push_back(line);
   }
   return res;
}

void
Lattice::addcellshapes()
{
   vd cclone=centermass_clone();
   for (icell ic=cells.begin()+1;ic!=cells.end();++ic){
      vd line;
      line.push_back(time);
      vd centerc=(*ic)->center();
      double dist2=(centerc[0]-cclone[0])*(centerc[0]-cclone[0])+(centerc[1]-cclone[1])*(centerc[1]-cclone[1]);
      double dist=sqrt(dist2);
      line.push_back(dist);
      double ux=(centerc[0]-cclone[0])/dist;
      double uy=(centerc[1]-cclone[1])/dist;
//      vd cshape=(*ic)->shape();
//      double angle=cshape[1];
      vd cshape=(*ic)->texture();
      double angle=cshape[4];
      double vx=cos(angle);
      double vy=sin(angle);
      line.push_back(atan((ux*vy-uy*vx)/(ux*vx+uy*vy)));
      line.push_back((*ic)->clone_index);
      cell_shapes.push_back(line);
   }
}

vd
Cell::shape()
{
   vd res;
   double mxx=0;
   double mxy=0;
   double mass=0;
   vd cmass=center();
   for (ivert iv=verts.begin();iv!=verts.end();iv++){
      double rx=(*iv)->x-cmass[0];
      double ry=(*iv)->y-cmass[1];
      mass+=rx*rx+ry*ry;
      mxx+=rx*rx-ry*ry;
      mxy+=2*rx*ry;
   }
   res.push_back(sqrt(mxx*mxx+mxy*mxy)/mass);
   if (fabs(mxy)<1e-4){
      if (mxx<0){
         res.push_back(M_PI/2);
      } else {
         res.push_back(0);
      }
   } else {
      if (mxy>0){
         res.push_back(atan(sqrt(mxx*mxx/(mxy*mxy)+1)-mxx/mxy));
      } else {
         res.push_back(atan(-sqrt(mxx*mxx/(mxy*mxy)+1)-mxx/mxy));
      }
   }
   return res;
}

vd
Cell::texture()
{
   double mxx,mxy,myy,mass;
   mxx=mxy=myy=mass=0;

   vd res;
   
   for (ivert iv=verts.begin();iv!=verts.end();++iv){
      unsigned int icelln=(*iv)->indexcell(this);
      Cell & cn=*((*iv)->pncell[(icelln+1)%3]);
      if (cn.clone_index!=-1){
         double dx=cn.cell_center[0]-cell_center[0];
         double dy=cn.cell_center[1]-cell_center[1];
         mxx+=dx*dx;
         mxy+=dx*dy;
         myy+=dy*dy;
         mass+=dx*dx+dy*dy;
      }
   }

   res.push_back(mxx);
   res.push_back(mxy);
   res.push_back(myy);

   double trace=mxx+myy;
   mxx-=trace/2.0;
   myy-=trace/2.0;
   res.push_back(2*sqrt(mxx*mxx+mxy*mxy)/mass);
   if (fabs(mxy)<1e-4){
      if (mxx<0){
         res.push_back(M_PI/2);
      } else {
         res.push_back(0);
      }
   } else {
      if (mxy>0){
         res.push_back(atan(sqrt(mxx*mxx/(mxy*mxy)+1)-mxx/mxy));
      } else {
         res.push_back(atan(-sqrt(mxx*mxx/(mxy*mxy)+1)-mxx/mxy));
      }
   }
   return res;
}

bool
Lattice::cbelongborder(Cell & c)
{
   for (ivert iv=c.verts.begin();iv!=c.verts.end();++iv){
      Vertex & v=**iv;
      if (belongborder(v)){
         return true;
      }
   }
   return false;
}


vpcell
Cell::neighbs()
{
   vpcell res;
   for (ivert iv=verts.begin();iv!=verts.end();iv++){
      Vertex & v=**iv;
      for (unsigned int i=0;i<3;i++){
         if (v.pncell[i]==this){
            res.push_back(v.pncell[(i+1)%3]);
         }
      }
   }
   return res;
}

   int
Cell::indexvertex(Vertex *pvertex)
{
   unsigned int i=0;
   for (ivert iv=verts.begin();iv!=verts.end();iv++){
      if (*iv==pvertex){
         return i;
      }
      i++;
   }
   //locerror("indexvertex", "The vertex does not belong to the cell");
   printf("The vertex does not belong to the cell\n");
   return -1;
}


/* Insert a new vertex in vertexlist */
   void
Cell::insvertex(Vertex *pnewvert, Vertex *pnextvert) /* The new vertex is at the first position */
{
   unsigned int ivert;
   if ((ivert=indexvertex(pnextvert)) < 0){
      printf("In insvertex ivert < 0\n");
      exit(1);
   };

   unsigned int nbvertices=verts.size();
   vpvert pvertexlisttemp;

   pvertexlisttemp.push_back(pnewvert);
   for (unsigned int n=0 ; n < nbvertices ; n++){
      pvertexlisttemp.push_back(verts[(ivert+n)%nbvertices]);
   }

   verts=pvertexlisttemp;
}


void
Cell::comp_area()
{
   Vertex *pvprev=*(verts.end()-1);
   Vertex *pv=*verts.begin();
   area=pvprev->x*pv->y-pvprev->y*pv->x;

   for (ivert iv=verts.begin()+1; iv!=verts.end(); iv++){
      pvprev=*(iv-1);
      pv=*iv;
      area += pvprev->x*pv->y-pvprev->y*pv->x;
   }
   area*=0.5;

   if (args_info.model_arg==1){
      sqrtarea=sqrt(fabs(area));
   }
}

void
Cell::comp_peri()
{
   peri=0;

   ivert iv;
   for (iv=verts.begin();iv!=verts.end()-1;iv++){
      peri+=dist(*iv,*(iv+1));
   }
   peri+=dist(verts[0], *iv);
}

   void
Cell::latticeprint(ofstream & file, unsigned int type, string sticker, double scalepres, string extension)
{
   double scale=20;
   double transx=300;
   double transy=300;

   if (extension=="svg"){
      file << "         <polygon\n";
      file << "            id=\"cell" << sticker << "\"\n";
   } else if (extension=="x3d"){
      file << "         <Shape>\n";
      file << "            <Extrusion crossSection='";
      for (civert civ=verts.begin();civ!=verts.end();++civ){
         file << (*civ)->x << " " << (*civ)->y << " ";
      }
      file << (*(verts.begin()))->x << " " << (*(verts.begin()))->y << "' scale='1 1 1 1' spine='0 0 0 0 0 " << thick << "'/>\n";
      file << "            <Appearance>\n";
   } else {
      cerr << "EROR : filetype not recognized : " << extension << endl;
   }

   switch (type){
      case 0:
         //         printdefects(pcell, pfile);
         break;
      case 1:
         printpressure(file,scalepres,extension);
         break;
      case 2:
         printthick(file, scalepres,extension);
         break;
      case 3:
         //         printgradcell(pcell, pfile);
         break;
      case 4:
         //         printampli(pcell, pfile);
         break;
      case 6:
         //         printgradcentr(pcell,pfile);
         break;
      case 5:
         break;
   }
   if (extension=="svg"){
      Vertex & vert=*(verts[0]);
      file << setprecision(4);
      file << "            points=\"" << vert.x*scale+transx << "," << vert.y*scale+transy;
      for (unsigned int ivert=1; ivert< verts.size(); ivert++){
         Vertex & vert2=*(verts[ivert]);
         file << " " << vert2.x*scale+transx << "," << vert2.y*scale+transy;
      }
      file << "\" />\n";

      // Variable boundary thickness
      vd text=texture();
      double angle=text[4];
      for (unsigned int ivert=1; ivert< verts.size(); ivert++){
         file << "         <line stroke-width=\"";
         Vertex & vert1=*(verts[ivert-1]);
         Vertex & vert2=*(verts[ivert]);
         double dx=vert2.x-vert1.x;
         double dy=vert2.y-vert1.y;

         double scal=cos(angle)*dx+sin(angle)*dy;
         double thickl=4*text[3]*scal*scal/(dx*dx+dy*dy);
         file << thickl << "\" x1=\"";
         file << vert1.x*scale+transx << "\" y1=\"" << vert1.y*scale+transy;
         file << "\" x2=\"";
         file << vert2.x*scale+transx << "\" y2=\"" << vert2.y*scale+transy;
         file << "\" />\n";
      }
      file << "         <line stroke-width=\"";
      Vertex & vert1=*(verts[verts.size()-1]);
      Vertex & vert2=*(verts[0]);
      double dx=vert2.x-vert1.x;
      double dy=vert2.y-vert1.y;

      double scal=cos(angle)*dx+sin(angle)*dy;
      double thickl=4*text[3]*scal*scal/(dx*dx+dy*dy);
      file << thickl << "\" x1=\"";
      file << vert1.x*scale+transx << "\" y1=\"" << vert1.y*scale+transy;
      file << "\" x2=\"";
      file << vert2.x*scale+transx << "\" y2=\"" << vert2.y*scale+transy;
      file << "\" />\n";

      if (clone_index>0){
         vd c=center();
         file << "         <circle cx=\"" << c[0]*scale+transx << "\" cy=\"" << c[1]*scale+transy << "\" r=\"5.0\" fill=\"black\" />\n";
      }
   } else if (extension=="x3d"){
      file << "            </Appearance>\n";
      file << "         </Shape>\n";
   }
}

/* Computes the centroid of one cell (considering vertices)
*/
   vd
Cell::center()
{
   vd center;

   center.push_back(0);
   center.push_back(0);
   for (ivert iv=verts.begin();iv!=verts.end();iv++){
      Vertex & v=*(*iv);
      center[0]+=v.x;
      center[1]+=v.y;
   }

   center[0]/=verts.size();
   center[1]/=verts.size();

   return center;
}

/* Computes the centroid of one cell (true center of mass)
*/
   vd
Cell::center_mass()
{
   vd center;

   center.push_back(0);
   center.push_back(0);
   for (ivert iv=verts.begin();iv!=verts.end()-1;iv++){
      Vertex & v=*(*iv);
      Vertex & vn=*(*(iv+1));
      center[0]+=(v.x+vn.x)*(v.x*vn.y-vn.x*v.y);
      center[1]+=(v.y+vn.y)*(v.x*vn.y-vn.x*v.y);
   }
   Vertex & v=*(*(verts.end()-1));
   Vertex & vn=*(*verts.begin());
   center[0]+=(v.x+vn.x)*(v.x*vn.y-vn.x*v.y);
   center[1]+=(v.y+vn.y)*(v.x*vn.y-vn.x*v.y);
   center[0]/=6*area;
   center[1]/=6*area;

   return center;
}

   void
Cell::printpressure(ofstream & file,double scale,string extension)
{
   double pres;

   pres=0;
   pres=pressure();

   if (pres<0){
      int color=255*(1+pres/scale);
      file << "            fill=\"rgb(255," << color << "," << color << ")\"\n";
   } else{
      int color=255*(1-pres/scale);
      file << "            fill=\"rgb(" << color << "," << color << ",255)\"\n";
   }
}

   void
Cell::printthick(ofstream & file,double scale,string extension)
{
   double pres=thick-thick_intr;
   if (extension=="svg"){
      if (pres>0){
         int color=255*(1-pres/scale);
         file << "            fill=\"rgb(255," << color << "," << color << ")\"\n";
      } else{
         int color=255*(1+pres/scale);
         file << "            fill=\"rgb(" << color << "," << color << ",255)\"\n";
      }
   } else if (extension=="x3d"){
      file << "               <Material diffuseColor='";
      if (pres>0){
         double color=1-pres/scale;
         file << "1.0 " << color << " " << color << "'/>\n";
      } else{
         double color=1+pres/scale;
         file << color << " " << color << " 1.0'/>\n";
      }
   }
}

   double
Cell::area_fct()
{
   double fct;

   fct =area-area_intr;
   if (area>0){
      return 0.5*fct*fct/area;
   } else {
      return -0.5*fct*fct/area;
   }
}

/* energy function of the area part
*/
   double
Cell::area_fctp()
{
   double diff=area-area_intr;
   if (area>0){
      return diff/area*(1-0.5*diff/area);
   } else {
      return diff/area*(-1+0.5*diff/area);
   }
}

   double
Cell::vol_fct()
{
   double vol=thick*area;
   double diff=vol-vol_intr;
//   return 0.5*diff*diff;
   if (area>0){
      return 0.5*diff*diff/vol;
   } else {
      return -0.5*diff*diff/vol;
   }
}

   double
Cell::myopol_en()
{
   double cosangle=1-angle*angle/2;
   double sinangle=angle-angle*angle*angle/6;
   double res=0;
   for (ivert iv=verts.begin()+1;iv!=verts.end();++iv){
      double dx=(*iv)->x-(*(iv-1))->x;
      double dy=(*iv)->y-(*(iv-1))->y;
      double scal=cosangle*dx+sinangle*dy;
      res+=scal*scal;
   }
   double dx1=(*verts.begin())->x-(*(verts.end()-1))->x;
   double dy1=(*verts.begin())->y-(*(verts.end()-1))->y;
   double scal1=cosangle*dx1+sinangle*dy1;
   res+=scal1*scal1;

   res*=0.5*ellipticity;
   return res;
}

/* energy function of the area part
*/
   double
Cell::vol_fctpa()
{
   double vol=area*thick;
   double diff=vol-vol_intr;
//   return (thick*area-vol_intr)*thick;
   if (area>0){
      return diff/vol*(thick-0.5*diff/area);
   } else {
      return -diff/vol*(thick-0.5*diff/area);
   }
}

/* energy function of the area part
*/
   double
Cell::vol_fctpe()
{
   double vol=area*thick;
   double diff=vol-vol_intr;
//   return (thick*area-vol_intr)*area;
   if (area>0){
      return diff/vol*(area-0.5*diff/thick);
   } else {
      return -diff/vol*(area-0.5*diff/thick);
   }
}

   double
Cell::thick_fct()
{
   double fct=0;
   for (ivert iv=verts.begin();iv!=verts.end();++iv){
      unsigned int icelln=(*iv)->indexcell(this);
      Cell & cn=*((*iv)->pncell[(icelln+1)%3]);
      if (cn.clone_index!=-1){
         double diff=(thick-cn.thick)/(thick+cn.thick);
         fct+=0.5*diff*diff;
      }
   }
   return fct;
}

   double
Cell::thick_fctp()
{
   double dfct=0;
   for (ivert iv=verts.begin();iv!=verts.end();++iv){
      unsigned int icelln=(*iv)->indexcell(this);
      Cell & cn=*((*iv)->pncell[(icelln+1)%3]);
      if (cn.clone_index!=-1){
         double diff=thick-cn.thick;
         double sum=thick+cn.thick;
         dfct+=diff/(sum*sum)*(1-diff/sum);
      }
   }
   return dfct;
}

   double
Cell::thick_1_fct()
{
   double diff;

   diff=thick-thick_intr;
   return 0.5*diff*diff;
}

   double
Cell::thick_1_fctp()
{
   return thick-thick_intr;
}

/* energy function of the perimeter part
*/
   double
Cell::rho_fct()
{
   return peri/sqrtarea;
}

/* derivative  of the perimeter part of the Hamiltonian
*/
   double
Cell::rho_fctpp()
{
   return 1/sqrtarea;
}

/* derivative  of the perimeter part of the Hamiltonian
*/
   double
Cell::rho_fctpa()
{
   return -0.5*peri/(area*sqrtarea);
}

/* energy function of the perimeter part
*/
   double
Cell::rho_fct2()
{
   if (area>0){
      return 0.5*peri*peri/area;
   } else {
      return -0.5*peri*peri/area;
   }
   return 0;
}

/* derivative  of the perimeter part of the Hamiltonian
*/
   double
Cell::rho_fctpp2()
{
   if (area>0){
      return peri/area;
   } else {
      return -peri/area;
   }
   return 0;
}

/* derivative  of the perimeter part of the Hamiltonian
*/
   double
Cell::rho_fctpa2()
{
   double ratio=peri/area;
   if (area>0){
      return -0.5*ratio*ratio;
   } else {
      return 0.5*ratio*ratio;
   }
   return 0;
}

   double
Cell::pressure()
{
   double areap=area_fctp();
   unsigned int i=1;
   double pressx=0;
   double pressy=0;
   double x_0=(*verts.begin())->x;
   double y_0=(*verts.begin())->y;
   for (ivert iv=verts.begin();iv!=verts.end();iv++){
      ivert ivnext;
      ivert ivprev;
      if (i==1){
         ivnext=iv+1;
         ivprev=verts.end()-1;
      } else if (i==verts.size()){
         ivnext=verts.begin();
         ivprev=iv-1;
      } else {
         ivnext=iv+1;
         ivprev=iv-1;
      }
      double vertx=(*iv)->x;
      double verty=(*iv)->y;
      double dxa,dya;
      dxa=dya=0;

      dxa += areap*0.5*((*ivnext) -> y - (*ivprev) -> y); 
      dya += areap*0.5*((*ivprev) -> x - (*ivnext) -> x);

      pressx+=(vertx-x_0)*dxa;
      pressy+=(verty-y_0)*dya;
      i++;
   }
   comp_area();
   return (pressx+pressy)*0.5/area;
}

Vertex::Vertex()
{
   for (unsigned int i=0;i<3;i++){
      pneighb.push_back(NULL);
      pncell.push_back(NULL);
   }
}

/* indexcell return the index of cell for a given vertex and cell */ 
   int
Vertex::indexcell(Cell * pcell)
{
   for (unsigned int i=0; i<3 ; i++ ){
      if ((pncell[i])==pcell) return i;
   }
   //locerror("indexcell", "the cell does not belong to the vertex");
   return -1;
}


double
dist(Vertex * v1,Vertex * v2)
{
   double dx=v1->x-v2->x;
   double dy=v1->y-v2->y;
   return sqrt(dx*dx+dy*dy);
}

double
dist(vd & v1,vd & v2)
{
   double dx=v1[0]-v2[0];
   double dy=v1[1]-v2[1];
   return sqrt(dx*dx+dy*dy);
}

double von_mises_dist(double k)
{
  if (k<1e-6) return 2*M_PI*gsl_rng_uniform(rng); 

  double result = 0.0;

  double a = 1.0 + sqrt(1 + 4.0 * (k * k));
  double b = (a - sqrt(2.0 * a))/(2.0 * k);
  double r = (1.0 + b * b)/(2.0 * b);

  double f=0;
  while (1)
  {
    double U1 = gsl_rng_uniform(rng);
    double z = cos(M_PI * U1);
    f = (1.0 + r * z)/(r + z);
    double c = k * (r - f);
    double U2 = gsl_rng_uniform(rng);

    if (c * (2.0 - c) - U2 > 0.0 || U2 <= c*exp(1.0-c)) break;
  }

  double theta=0;
  double u3=gsl_rng_uniform(rng);
  if (u3 > 0.5) theta=acos(f);
  else theta=-acos(f);

  return theta;
}
