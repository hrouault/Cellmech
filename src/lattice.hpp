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

#ifndef LATTICE_HPP
#define LATTICE_HPP 1


#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include "vectortypes.hpp"

using namespace std;


class Lattice;
class Vertex;
class Cell;
class Parameters;

typedef vector<Vertex *> vpvert;
typedef vpvert::iterator ivert;
typedef vpvert::const_iterator civert;

typedef vector<Cell *> vpcell;
typedef vpcell::iterator icell;
typedef vpcell::const_iterator cicell;

class Parameters
{// Modèle et constantes
   public:
      // Mechanical parameters
      int model;
      double alpha;
      // Parameters for 3D model
      double beta;
      double gamma;
      double pol_force;
};

class Lattice
{// Un réseau de cellules
   public:
      double time;

      vpvert verts;
      vpcell cells;

      unsigned int nbt1;
      unsigned int nbt1closetodiv;

      int dirborderdir;

      vvd division_dirs;

      vvd cell_shapes;

      bool exploded;

      // division postitions
      FILE * fdivpos;

      Lattice();
      Lattice(unsigned int size);
      void shrink_datafirst();

      void latticeprint(string filename, unsigned int type);

      int check_t1();
      void t1transform(Vertex * pvert, unsigned int neighbor);

      void update_pa();
      void update_centers();
      void update_texture();
      
      void division(Cell & c);
      void division(Cell & c, unsigned int iv1, unsigned int iv2);
      void division(Cell & c, double theta);
      void division_alongbound(Cell & c);
      void division_polarclone(Cell & c);
      void division_randdir(Cell & c);
      void delvertex(Cell & cell, Vertex * pvert);

      vd center();

      void diffclone(unsigned int nb_cells);
      void diffclone();
      void diffcompart();

      Cell * closest_center();
      Cell * closest_center_onborder();
      bool belongborder(Vertex & v);
      bool cbelongborder(Cell & c);

      void printneighbors(ofstream & file);
      void printcompartline(ofstream & file);

      void reinittime();

      double boundarylength();

      unsigned int clone_size();
      double anisotropy_clone();
      double meanpressure_clone();
      vd centermass_clone();

      void adddivdir();

      void addcellshapes();

      vvd angles_aroundclone();
};


class Vertex
{
   public:
      double x;
      double y;

      vpvert pneighb;
      vpcell pncell;

      Vertex();
      int indexcell(Cell * pcell);
};

class Cell
{
   public:
      unsigned int index;

      vpvert verts;

      double thick;
      double thick_intr;

      double area;
      double area_intr;
      double vol_intr;

      double sqrtarea;
      double peri;

      // Clone properties
      int clone_index;

      int compart;

      // Division properties
      double cyc_shape; // Shape parameter of the gamma distribution
      double cyc_scale; // Time scale of the cell cycle
      double t_next_div;

      bool just_div;

      vd cell_center;

      double ellipticity;
      double angle;

      void insvertex(Vertex *pnewvert, Vertex *pnextvert);
      int indexvertex(Vertex *pvertex);

      vd center();
      vd center_mass();

      void latticeprint(ofstream & file, unsigned int type, string sticker, double scalepres, string extension);
      void printpressure(ofstream & file, double scale,string extension);
      void printthick(ofstream & file,double scale,string extension);

      void comp_area();
      void comp_peri();

      double area_fct();
      double area_fctp();
      double vol_fct();
      double vol_fctpa();
      double vol_fctpe();
      double rho_fct();
      double rho_fctpp();
      double rho_fctpa();
      
      double rho_fct2();
      double rho_fctpp2();
      double rho_fctpa2();

      double thick_fct();
      double thick_fctp();
      double thick_1_fct();
      double thick_1_fctp();

      double myopol_en();

      double pressure();

      vpcell neighbs();

      vd shape();
      vd texture();
};

class Simulation
{// Lattice et paramatres
   public:
      unsigned int index;

      Lattice lattice;
      Parameters params;

      void relax();
      void optimize3d();
      void optimize();
      void div_next();
};

typedef vector<Simulation> vsim;
typedef vsim::iterator ivsim;

int vertex_from_coord(unsigned int row, unsigned int column);

vpvert allocverts(int nbverts);
vpcell alloccells(int nbcells);

void locerror(string function, string message);

double dist(Vertex * v1,Vertex * v2);
double dist(vd & v1,vd & v2);

double von_mises_dist(double k);

#endif /* LATTICE_HPP */
