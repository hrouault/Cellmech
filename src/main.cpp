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

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>
#include <pthread.h>

#include "cmdline.h"
#include "main.hpp"
#include "lattice.hpp"


struct gengetopt_args_info args_info;
gsl_rng * rng;

void
cleansimus(vsim & vs)
{
   vint toerase;
   int i=0;
   for (ivsim is=vs.begin();is!=vs.end();++is){
      if (is->lattice.exploded){
         toerase.push_back(i);
      }
      i++;
   }
   for (ivint iv=toerase.end()-1;iv!=toerase.begin()-1;++iv){
      vs.erase(vs.begin()+*iv);
   }
}

void *
evolsimu_thr(void * simu)
{
   Simulation & sim=*(Simulation *)simu;
   Lattice & lat=sim.lattice;
   lat.nbt1=0;
   lat.nbt1closetodiv=0;

   // Division position recording
   stringstream fndivpos;
   fndivpos << "divpos-" << sim.index << ".dat";
   lat.fdivpos = fopen(fndivpos.str().c_str(),"w");

   int nbcelltotinit=lat.cells.size();
   for (unsigned int j=nbcelltotinit-1;j<1000;j++){
      if (!lat.exploded){
         Cell * pc;
         if (gsl_rng_uniform(rng)<0.2 && j<200){
            pc=lat.closest_center_onborder();
         } else {
            unsigned int icell=gsl_rng_uniform(rng)*(lat.cells.size()-1)+1;
            pc=lat.cells[icell];
         }

         lat.division_randdir(*pc);
         sim.relax();
         if (j%100==0){
            lat.update_centers();
            stringstream filename;
            filename << "simu-" << sim.index << "-";
            filename.fill('0');
            filename.width(5);
            filename << j << "-" << sim.lattice.dirborderdir;
            filename << ".svg";
            if (args_info.third_dim_flag){
               sim.lattice.latticeprint(filename.str(),2);
            } else {
               sim.lattice.latticeprint(filename.str(),1);
            }
         }
         /*
         cout << "Start of records" << endl;
         for (icell ic=lat.cells.begin()+1;ic!=lat.cells.end();++ic){
            cout << (*ic)->area << " " << (*ic)->peri << " " << (*ic)->thick;
            cout << " " << (*ic)->vol_fct() << " " << (*ic)->rho_fct2()+(*ic)->thick_1_fct() << " " << (*ic)->thick_1_fctp();
            cout << endl;
         }
         cout << "End of records" << endl;
         */
      }
   }

   if (!lat.exploded){
      lat.reinittime();
//      lat.diffclone();
//      lat.diffcompart();
   }

   for (unsigned int j=1000;j<3000;j++){
      if (!lat.exploded){
         //      Cell * pc;
         //      unsigned int icell=gsl_rng_uniform(rng)*(lat.cells.size()-1)+1;
         //      pc=lat.cells[icell];

         //      lat.division(*pc);
         sim.div_next();
         sim.relax();
         lat.adddivdir();

         /*       if (lat.clone_size()>10){
          * 
          *          stringstream filename2;
          *          filename2 << "simu-" << sim.index << ".dat";
          *          ofstream file(filename2.str().c_str());
          *          sim.lattice.printneighbors(file);
          * 
          *          stringstream filename;
          *          filename << "simu-" << sim.index << ".svg";
          *          lat.svgprint(filename.str(),1);
          * 
          *          break;
          *       }
          */

         if (j%100==0){
            lat.update_centers();
            stringstream filename;
            filename << "simu-" << sim.index << "-";
            filename.fill('0');
            filename.width(4);
            filename << j << "-" << sim.lattice.dirborderdir;
            if (args_info.third_dim_flag){
               filename << ".svg";
            } else {
               filename << ".svg";
            }
            if (args_info.third_dim_flag){
               sim.lattice.latticeprint(filename.str(),2);
            } else {
               sim.lattice.latticeprint(filename.str(),1);
            }
            // Thickness
            stringstream filename2;
            filename2 << "simu-" << sim.index << "-";
            filename2.fill('0');
            filename2.width(4);
            filename2 << j << "-pro.dat";
            ofstream filecellpro(filename2.str().c_str());
            vd clat=lat.center();
            for (icell ic=lat.cells.begin()+1;ic!=lat.cells.end();++ic){
               vd cc=(*ic)->cell_center;
               filecellpro << dist(clat,cc) << " " << (*ic)->thick << " " << (*ic)->area << endl;
            }
            filecellpro.close();
            // Angle relative to radius
            stringstream fnangles;
            fnangles << "simu-" << sim.index << "-";
            fnangles.fill('0');
            fnangles.width(4);
            fnangles << j << "-angle.dat";
            ofstream fcangles(fnangles.str().c_str());
            for (icell ic=lat.cells.begin()+1;ic!=lat.cells.end();++ic){
               vd cc=(*ic)->cell_center;
               double dist2=(cc[0]-clat[0])*(cc[0]-clat[0])+(cc[1]-clat[1])*(cc[1]-clat[1]);
               double distfromc=sqrt(dist2);
               double ux=(cc[0]-clat[0])/distfromc;
               double uy=(cc[1]-clat[1])/distfromc;
               vd cshape=(*ic)->texture();
               double angle=cshape[4];
               double vx=cos(angle);
               double vy=sin(angle);
               double relangle=atan((ux*vy-uy*vx)/(ux*vx+uy*vy));
               fcangles << dist(clat,cc) << " " << relangle << " " << cshape[3] << endl;
            }
            fcangles.close();
         }


         if (j%200==0){
            lat.addcellshapes();
         }

         /*       stringstream filename;
          *       filename << "simu-" << sim.index << "-";
          *       filename.fill('0');
          *       filename.width(4);
          *       filename << j << ".svg";
          *       lat.svgprint(filename.str(),1);
          */
      }
   }
   for (unsigned int k=0;k<30;k++){
      sim.optimize3d(); 
   }
   lat.update_centers();
   stringstream filename;
   filename << "simu-before-feedback-" << sim.index << "-";
   filename << sim.lattice.dirborderdir;
   if (args_info.third_dim_flag){
      filename << ".svg";
   } else {
      filename << ".svg";
   }
   if (args_info.third_dim_flag){
      sim.lattice.latticeprint(filename.str(),2);
   } else {
      sim.lattice.latticeprint(filename.str(),1);
   }
   // Thickness
   stringstream filename2;
   filename2 << "simu-before-feedback-" << sim.index << "-";
   filename2 << "pro.dat";
   ofstream filecellpro(filename2.str().c_str());
   vd clat=lat.center();
   for (icell ic=lat.cells.begin()+1;ic!=lat.cells.end();++ic){
      vd cc=(*ic)->cell_center;
      filecellpro << dist(clat,cc) << " " << (*ic)->thick << endl;
   }
   filecellpro.close();
   // Angle relative to radius
   stringstream fnangles;
   fnangles << "simu-before-feedback-" << sim.index << "-";
   fnangles << "angle.dat";
   ofstream fcangles(fnangles.str().c_str());
   for (icell ic=lat.cells.begin()+1;ic!=lat.cells.end();++ic){
      vd cc=(*ic)->cell_center;
      double dist2=(cc[0]-clat[0])*(cc[0]-clat[0])+(cc[1]-clat[1])*(cc[1]-clat[1]);
      double distfromc=sqrt(dist2);
      double ux=(cc[0]-clat[0])/distfromc;
      double uy=(cc[1]-clat[1])/distfromc;
      vd cshape=(*ic)->texture();
      double angle=cshape[4];
      double vx=cos(angle);
      double vy=sin(angle);
      double relangle=atan((ux*vy-uy*vx)/(ux*vx+uy*vy));
      fcangles << dist(clat,cc) << " " << relangle << " " << cshape[3] << endl;
   }
   fcangles.close();
   sim.params.pol_force=0.2;
   for (unsigned int k=0;k<30;k++){
      sim.optimize3d(); 
   }
   lat.update_centers();
   stringstream filename_2;
   filename_2 << "simu-after-feedback-" << sim.index << "-";
   filename_2 << sim.lattice.dirborderdir;
   if (args_info.third_dim_flag){
      filename_2 << ".svg";
   } else {
      filename_2 << ".svg";
   }
   if (args_info.third_dim_flag){
      sim.lattice.latticeprint(filename_2.str(),2);
   } else {
      sim.lattice.latticeprint(filename_2.str(),1);
   }
   // Thickness
   stringstream filename2_2;
   filename2_2 << "simu-after-feedback-" << sim.index << "-";
   filename2_2 << "-pro.dat";
   ofstream filecellpro_2(filename2_2.str().c_str());
   clat=lat.center();
   for (icell ic=lat.cells.begin()+1;ic!=lat.cells.end();++ic){
      vd cc=(*ic)->cell_center;
      filecellpro_2 << dist(clat,cc) << " " << (*ic)->thick << endl;
   }
   filecellpro_2.close();
   // Angle relative to radius
   stringstream fnangles_2;
   fnangles_2 << "simu-after-feedback-" << sim.index << "-";
   fnangles_2 << "-angle.dat";
   ofstream fcangles_2(fnangles_2.str().c_str());
   for (icell ic=lat.cells.begin()+1;ic!=lat.cells.end();++ic){
      vd cc=(*ic)->cell_center;
      double dist2=(cc[0]-clat[0])*(cc[0]-clat[0])+(cc[1]-clat[1])*(cc[1]-clat[1]);
      double distfromc=sqrt(dist2);
      double ux=(cc[0]-clat[0])/distfromc;
      double uy=(cc[1]-clat[1])/distfromc;
      vd cshape=(*ic)->texture();
      double angle=cshape[4];
      double vx=cos(angle);
      double vy=sin(angle);
      double relangle=atan((ux*vy-uy*vx)/(ux*vx+uy*vy));
      fcangles_2 << dist(clat,cc) << " " << relangle << " " << cshape[3] << endl;
   }
   fcangles_2.close();
   for (unsigned int k=0;k<30;k++){
      sim.relax(); 
   }
   lat.update_centers();
   stringstream filename_3;
   filename_3 << "simu-after-T1-" << sim.index << "-";
   filename_3 << sim.lattice.dirborderdir;
   if (args_info.third_dim_flag){
      filename_3 << ".svg";
   } else {
      filename_3 << ".svg";
   }
   if (args_info.third_dim_flag){
      sim.lattice.latticeprint(filename_3.str(),2);
   } else {
      sim.lattice.latticeprint(filename_3.str(),1);
   }
   // Thickness
   stringstream filename2_3;
   filename2_3 << "simu-after-T1-" << sim.index << "-";
   filename2_3 << "-pro.dat";
   ofstream filecellpro_3(filename2_3.str().c_str());
   clat=lat.center();
   for (icell ic=lat.cells.begin()+1;ic!=lat.cells.end();++ic){
      vd cc=(*ic)->cell_center;
      filecellpro_3 << dist(clat,cc) << " " << (*ic)->thick << endl;
   }
   filecellpro_3.close();
   // Angle relative to radius
   stringstream fnangles_3;
   fnangles_3 << "simu-after-T1-" << sim.index << "-";
   fnangles_3 << "-angle.dat";
   ofstream fcangles3_3(fnangles_3.str().c_str());
   for (icell ic=lat.cells.begin()+1;ic!=lat.cells.end();++ic){
      vd cc=(*ic)->cell_center;
      double dist2=(cc[0]-clat[0])*(cc[0]-clat[0])+(cc[1]-clat[1])*(cc[1]-clat[1]);
      double distfromc=sqrt(dist2);
      double ux=(cc[0]-clat[0])/distfromc;
      double uy=(cc[1]-clat[1])/distfromc;
      vd cshape=(*ic)->texture();
      double angle=cshape[4];
      double vx=cos(angle);
      double vy=sin(angle);
      double relangle=atan((ux*vy-uy*vx)/(ux*vx+uy*vy));
      fcangles3_3 << dist(clat,cc) << " " << relangle << " " << cshape[3] << endl;
   }
   fcangles3_3.close();

   fclose(lat.fdivpos);

   return 0;
}

/* Main function
 * Only call the content functions
 */
   int
main(int argc,char** argv)
{

   if ( cmdline_parser(argc,argv, & args_info)!=0) 
      exit(1);

   unsigned long int seed;
   FILE * devrandom;

   rng=gsl_rng_alloc(gsl_rng_default);

   if (args_info.rand_seed_flag){
      devrandom=fopen("/dev/urandom","r");
      fread(&seed,sizeof(seed),1,devrandom);
      fclose(devrandom);
      gsl_rng_set(rng,seed);
   }

   vsim vs;
   for (unsigned int i=0;i<args_info.nblat_arg;i++){
      Simulation simu;
      simu.lattice=Lattice(4);
      simu.index=i;

      Parameters & par=simu.params;
      par.alpha=args_info.alpha_arg;
      par.beta=args_info.beta_arg;
      par.gamma=args_info.gamma_arg;
      par.pol_force=args_info.pol_force_arg;
      par.model=args_info.model_arg;

      simu.relax();

      simu.lattice.update_pa();

      vs.push_back(simu);
   }

   vthds vt;
   for (ivsim is=vs.begin();is!=vs.end();++is){
      pthread_t idth;
      Simulation & simu=*is;
      pthread_create(&idth,NULL,evolsimu_thr,&simu);
      vt.push_back(idth);
   }
   for (ivthds it=vt.begin();it!=vt.end();it++){
      pthread_join(*it,NULL);
   }

   cleansimus(vs);

   ofstream filestat("stats.dat");
   filestat << "Clone_anisotropy\tClone_pressure\tClone_connectedcomponents" << endl;
   for (ivsim is=vs.begin();is!=vs.end();++is){
      Simulation & simu=*is;
      filestat << simu.lattice.anisotropy_clone() << "\t" << -simu.lattice.meanpressure_clone();
      filestat << endl;
      stringstream filename;
      filename << "simu-" << simu.index << "-" << simu.lattice.dirborderdir << "-";
      filename << ".svg";
      if (args_info.third_dim_flag){
         simu.lattice.latticeprint(filename.str(),2);
      } else {
         simu.lattice.latticeprint(filename.str(),1);
      }
   }
   filestat.close();

   ofstream filecellshapes("cellstats.dat");
   for (ivsim is=vs.begin();is!=vs.end();++is){
      Simulation & simu=*is;
      for (icell ic=simu.lattice.cells.begin()+1;ic!=simu.lattice.cells.end();++ic){
         if ((*ic)->clone_index){
//            vd shape=(*ic)->shape();
//            filecellshapes << shape[0] << " " << shape[1] << endl;
            vd shape=(*ic)->texture();
            filecellshapes << shape[3] << " " << shape[4] << endl;
         }
      }
      filecellshapes << endl;
   }
   filecellshapes.close();

   ofstream filecellthicks("cellthicks.dat");
   for (ivsim is=vs.begin();is!=vs.end();++is){
      Simulation & simu=*is;
      for (icell ic=simu.lattice.cells.begin()+1;ic!=simu.lattice.cells.end();++ic){
            filecellthicks << (*ic)->thick << endl;
      }
      filecellthicks << endl;
   }
   filecellthicks.close();

   ofstream filecellareas("cellstats-area.dat");
   for (ivsim is=vs.begin();is!=vs.end();++is){
      Simulation & simu=*is;
      vd cent=simu.lattice.center();
      for (icell ic=simu.lattice.cells.begin()+1;ic!=simu.lattice.cells.end();++ic){
         if (simu.lattice.cbelongborder(**ic)==false){
            (*ic)->comp_area();
            vd cellc=(*ic)->cell_center;
            filecellareas << dist(cent,cellc) << " " <<  (*ic)->area << " " << (*ic)->verts.size() << endl;
         }
      }
      filecellareas << endl;
   }
   filecellareas.close();

   ofstream filecelledges("cellstats-edges.dat");
   for (ivsim is=vs.begin();is!=vs.end();++is){
      Simulation & simu=*is;
      for (icell ic=simu.lattice.cells.begin()+1;ic!=simu.lattice.cells.end();++ic){
         Cell *pcell=*ic;
         filecelledges << dist(pcell->verts[0],*(pcell->verts.end()-1)) << endl;
         for (ivert iv=pcell->verts.begin();iv!=pcell->verts.end()-1;++iv){
            filecelledges << dist(*iv,*(iv+1)) << endl;
         }
      }
      filecelledges << endl;
   }
   filecelledges.close();

   ofstream filearoundcl("cellstats-aroundclone.dat");
   for (ivsim is=vs.begin();is!=vs.end();++is){
      Simulation & simu=*is;
      for (ivvd ian=simu.lattice.cell_shapes.begin();ian!=simu.lattice.cell_shapes.end();++ian){
         vd & line=*ian;
         for (ivd iline=line.begin();iline!=line.end();++iline){
            filearoundcl << *iline << " " ;
         }
         filearoundcl << endl;
      }
      filearoundcl << endl;
   }
   filearoundcl.close();

   ofstream filedivs("divstats.dat");
   for (ivsim is=vs.begin();is!=vs.end();++is){
      Simulation & simu=*is;
      for (ivvd id=simu.lattice.division_dirs.begin();id!=simu.lattice.division_dirs.end();++id){
         vd & coords=*id;
         double dx=coords[2]-coords[0];
         double dy=coords[3]-coords[1];
         double norm=dx*dx+dy*dy;
         dx/=norm;
         dy/=norm;
         if (fabs(dx)<1e-3){
            filedivs << M_PI/2 << endl;
         } else {
            filedivs << atan(dy/dx) << endl;
         }
      }
      filedivs << endl;
   }
   filedivs.close();

   return 0;
}
