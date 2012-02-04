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



#include <iostream>
#include <gsl/gsl_multimin.h>

#include "lattice.hpp"
#include "main.hpp"


   double
energy(const gsl_vector *x,void *params)
{   /* Compute the energy */
   
   double h=0;
   Simulation & simu=*(Simulation *)params;
   Lattice & lat=simu.lattice;
   Parameters & par=simu.params;

   unsigned int i=0;
   for (civert iv=lat.verts.begin();iv!=lat.verts.end();++iv){
      (*iv)->x=gsl_vector_get(x,2*i);
      (*iv)->y=gsl_vector_get(x,2*i+1);
      i++;
   }
   lat.update_pa();
   for (cicell ic=lat.cells.begin()+1;ic!=lat.cells.end();++ic){
      Cell & c=**ic;
      if (args_info.model_arg==1){
         h+=c.area_fct()+par.alpha*c.rho_fct();
      } else if (args_info.model_arg==2){
         h+=c.area_fct()+par.alpha*c.rho_fct2();
      }
//      else if (args_info.model_arg==2){
//         h+=area_fct2(pcell)+((double *)params)[0]*rho_fct2(pcell);
//      } else if (args_info.model_arg==3){
//         h+=area_fct2(pcell)+((double *)params)[0]*rho_fct(pcell);
//      }
//      if (args_info.rand_force_flag){
//         h+=drift_fct(pcell);
//      }
   }

   return h;
}

   void
denergy(const gsl_vector *x, void *params, gsl_vector *df)
{
   /* Compute the gradient of the energy */
   Simulation & simu=*(Simulation *)params;
   Lattice & lat=simu.lattice;
   Parameters & par=simu.params;

   unsigned int i=0;
   for (civert iv=lat.verts.begin();iv!=lat.verts.end();iv++){
      (*iv)->x=gsl_vector_get(x,2*i);
      (*iv)->y=gsl_vector_get(x,2*i+1);
      i++;
   }
   lat.update_pa();
   Cell & c0=**lat.cells.begin();
   i=0;
   for (civert iv=lat.verts.begin();iv!=lat.verts.end();iv++){
      Vertex & v=**iv;
      double vertx=v.x;
      double verty=v.y;
      double dx,dy,dxa,dya,dxd,dyd;
      dx=dy=dxa=dya=dxd=dyd=0;
      for (unsigned int j=0; j<3; j++){
         Cell & c    = *v.pncell[j];
         Vertex & vnext = *v.pneighb[j];
         Vertex & vprev = *v.pneighb[(j+1)%3];

         double prefact=0;

         /* Displacement along x and y due to perimeter extension */
         if (&c!=&c0){
            double locdist=dist(&vnext,&v);
            if (args_info.model_arg==1){
               prefact = par.alpha*c.rho_fctpp()/locdist;
            } else if (args_info.model_arg==2){
               prefact = par.alpha*c.rho_fctpp2()/locdist;
            }
            //            else if (args_info.model_arg==2){
            //               prefact = par.alpha*rho_fct2pp(pcell)/locdist;
            //            } else if (args_info.model_arg==3){
            //               prefact = par.alpha*rho_fctpp(pcell)/locdist;
            //            }
            dx += (vertx - vnext.x)*prefact;
            dy += (verty - vnext.y)*prefact;

            locdist=dist(&vprev,&v);
            if (args_info.model_arg==1){
               prefact = par.alpha*c.rho_fctpp()/locdist;
            } else if (args_info.model_arg==2){
               prefact = par.alpha*c.rho_fctpp2()/locdist;
            }
            //            } else if (args_info.model_arg==3){
            //               prefact = ((double *)params)[0]*rho_fctpp(pcell)/locdist;
            //            }
            dx += (vertx - vprev.x)*prefact;
            dy += (verty - vprev.y)*prefact;

            /* Displacement along x and y due to area extension */
            if (args_info.model_arg==1){
               prefact = c.area_fctp()+par.alpha*c.rho_fctpa();
            } else if (args_info.model_arg==2){
               prefact = c.area_fctp()+par.alpha*c.rho_fctpa2();
            }
            //            else if (args_info.model_arg==2){
            //               prefact = area_fct2p(pcell)+((double *)params)[0]*rho_fct2pa(pcell);
            //            } else if (args_info.model_arg==3){
            //               prefact = area_fct2p(pcell)+((double *)params)[0]*rho_fctpa(pcell);
            //            }
            dxa += prefact*0.5*(vnext.y - vprev.y); 
            dya += prefact*0.5*(vprev.x - vnext.x);

            /* Random force */
            //            if (args_info.rand_force_flag){
            //               dxd+=-pcell->randfx/pcell->nb_vertices;
            //               dyd+=-pcell->randfy/pcell->nb_vertices;
            //            }

         }
      }

      /*      if (args_info.rand_force_flag){
              gsl_vector_set(df, 2*i, dx+dxa+dxd);
              gsl_vector_set(df, 2*i+1, dy+dya+dyd);
              } else {*/
      gsl_vector_set(df, 2*i, dx+dxa);
      gsl_vector_set(df, 2*i+1, dy+dya);
      /*        }*/
      i++;
   }
}

   void
enerdener(const gsl_vector *x, void *params, double *f, gsl_vector *df)
{
   /* Compute both the gradient and the energy */
   /* Compute the gradient of the energy */
   Simulation & simu=*(Simulation *)params;
   Lattice & lat=simu.lattice;
   Parameters & par=simu.params;

   unsigned int i=0;
   for (civert iv=lat.verts.begin();iv!=lat.verts.end();iv++){
      (*iv)->x=gsl_vector_get(x,2*i);
      (*iv)->y=gsl_vector_get(x,2*i+1);
      i++;
   }
   lat.update_pa();

   Cell & c0=**lat.cells.begin();
   i=0;
   for (ivert iv=lat.verts.begin();iv!=lat.verts.end();iv++){
      Vertex & v=**iv;
      double vertx=v.x;
      double verty=v.y;
      double dx,dy,dxa,dya,dxd,dyd;
      dx=dy=dxa=dya=dxd=dyd=0;
      for (unsigned int j=0; j<3; j++){
         Cell & c = *v.pncell[j];
         Vertex & vnext = *v.pneighb[j];
         Vertex & vprev = *v.pneighb[(j+1)%3];

         double prefact=0;
         double locdist=0;

         /* Displacement along x and y due to perimeter extension */
         if (&c!=&c0){
            locdist=dist(&vnext,&v);
            if (args_info.model_arg==1){
               prefact = par.alpha*c.rho_fctpp()/locdist;
            } else if (args_info.model_arg==2){
               prefact = par.alpha*c.rho_fctpp2()/locdist;
            }
//            else if (args_info.model_arg==2){
//               prefact = ((double *)params)[0]*rho_fct2pp(pcell)/locdist;
//            } else if (args_info.model_arg==3){
//               prefact = ((double *)params)[0]*rho_fctpp(pcell)/locdist;
//            }
            dx += (vertx - vnext.x)*prefact;
            dy += (verty - vnext.y)*prefact;

            locdist=dist(&vprev,&v);
            if (args_info.model_arg==1){
               prefact = par.alpha*c.rho_fctpp()/locdist;
            } else if (args_info.model_arg==2){
               prefact = par.alpha*c.rho_fctpp2()/locdist;
            }
//            else if (args_info.model_arg==2){
//               prefact = ((double *)params)[0]*rho_fct2pp(pcell)/locdist;
//            } else if (args_info.model_arg==3){
//               prefact = ((double *)params)[0]*rho_fctpp(pcell)/locdist;
//            }
            dx += (vertx - vprev.x)*prefact;
            dy += (verty - vprev.y)*prefact;

            /* Displacement along x and y due to area extension */
            if (args_info.model_arg==1){
               prefact = c.area_fctp()+par.alpha*c.rho_fctpa();
            } else if (args_info.model_arg==2){
               prefact = c.area_fctp()+par.alpha*c.rho_fctpa2();
            }
//            else if (args_info.model_arg==2){
//               prefact = area_fct2p(pcell)+((double *)params)[0]*rho_fct2pa(pcell);
//            } else if (args_info.model_arg==3){
//               prefact = area_fct2p(pcell)+((double *)params)[0]*rho_fctpa(pcell);
//            }
            dxa += prefact*0.5*(vnext.y - vprev.y); 
            dya += prefact*0.5*(vprev.x - vnext.x);

//            if (args_info.rand_force_flag){
//               dxd+=-pcell->randfx/pcell->nb_vertices;
//               dyd+=-pcell->randfy/pcell->nb_vertices;
//            }

         }
      }
//      if (args_info.rand_force_flag){
//         gsl_vector_set(df, 2*i, dx+dxa+dxd);
//         gsl_vector_set(df, 2*i+1, dy+dya+dyd);
//      } else {
         gsl_vector_set(df, 2*i, dx+dxa);
         gsl_vector_set(df, 2*i+1, dy+dya);
//      }
      i++;
   }
   /*
      gsl_vector_set(df, 2*nb_vertex_tot, 0.0);
      for (i=1;i<nb_cell_tot;i++){
      pcell=web_dual[i];
      gsl_vector_set(df, 2*nb_vertex_tot+i, ((double *)params)[1]*thick_fctp(pcell));
      }*/
   double h=0;
   for (cicell ic=lat.cells.begin()+1;ic!=lat.cells.end();++ic){
      Cell & c=**ic;
      if (args_info.model_arg==1){
         h+=c.area_fct()+par.alpha*c.rho_fct();
      } else if (args_info.model_arg==2){
         h+=c.area_fct()+par.alpha*c.rho_fct2();
      }
//         h+=((double *)params)[0]*rho_fct2(pcell)+area_fct2(pcell);
//      } else if (args_info.model_arg==3){
//         h+=((double *)params)[0]*rho_fct(pcell)+area_fct2(pcell);
//      }
//      if (args_info.rand_force_flag){
//         h+=drift_fct(pcell);
//      }
   }
   *f=h;
}

// Functions for 3D simulations
   double
energy3d(const gsl_vector *x,void *params)
{   /* Compute the energy */
   double h=0;

   Simulation & simu=*(Simulation *)params;
   Lattice & lat=simu.lattice;
   Parameters & par=simu.params;

   unsigned int i=0;
   for (civert iv=lat.verts.begin();iv!=lat.verts.end();++iv){
      (*iv)->x=gsl_vector_get(x,2*i);
      (*iv)->y=gsl_vector_get(x,2*i+1);
      i++;
   }
   i=0;
   for (cicell iv=lat.cells.begin()+1;iv!=lat.cells.end();++iv){
      (*iv)->thick=gsl_vector_get(x,2*lat.verts.size()+i);
      i++;
   }

   lat.update_pa();
   for (cicell ic=lat.cells.begin()+1;ic!=lat.cells.end();++ic){
      Cell & c=**ic;
      if (args_info.model_arg==1){
         h+=c.vol_fct()+par.alpha*c.rho_fct();
      } else if (args_info.model_arg==2){
         h+=c.vol_fct()+par.alpha*c.rho_fct2();
      }
//      else if (args_info.model_arg==2){
//         h+=area_fct2(pcell)+((double *)params)[0]*rho_fct2(pcell);
//      } else if (args_info.model_arg==3){
//         h+=area_fct2(pcell)+((double *)params)[0]*rho_fct(pcell);
//      }
//      if (args_info.rand_force_flag){
//         h+=drift_fct(pcell);
//      }
      h+=par.beta*c.thick_1_fct();
      h+=par.gamma*c.thick_fct();

      // Effect of MyoII polarization
      h+=par.pol_force*c.myopol_en();
   }

   return h;
}

   void
denergy3d(const gsl_vector *x, void *params, gsl_vector *df)
{
   /* Compute the gradient of the energy */
   Simulation & simu=*(Simulation *)params;
   Lattice & lat=simu.lattice;
   Parameters & par=simu.params;

   unsigned int i=0;
   for (civert iv=lat.verts.begin();iv!=lat.verts.end();iv++){
      (*iv)->x=gsl_vector_get(x,2*i);
      (*iv)->y=gsl_vector_get(x,2*i+1);
      i++;
   }
   i=0;
   for (cicell iv=lat.cells.begin()+1;iv!=lat.cells.end();++iv){
      (*iv)->thick=gsl_vector_get(x,2*lat.verts.size()+i);
      i++;
   }
   lat.update_pa();

   Cell & c0=**lat.cells.begin();
   i=0;
   for (civert iv=lat.verts.begin();iv!=lat.verts.end();iv++){
      Vertex & v=**iv;
      double vertx=v.x;
      double verty=v.y;
      double dx,dy,dxa,dya,dxd,dyd,dxm,dym;
      dx=dy=dxa=dya=dxd=dyd=dxm=dym=0;
      for (unsigned int j=0; j<3; j++){
         Cell & c    = *v.pncell[j];
         Cell & cellp    = *v.pncell[(j+2)%3];
         Vertex & vnext = *v.pneighb[j];
         Vertex & vprev = *v.pneighb[(j+1)%3];

         double prefact=0;

         /* Displacement along x and y due to perimeter extension */
         if (&c!=&c0){
            double locdist=dist(&vnext,&v);
            if (args_info.model_arg==1){
               prefact = par.alpha*c.rho_fctpp()/locdist;
            } else if (args_info.model_arg==2){
               prefact = par.alpha*c.rho_fctpp2()/locdist;
            }
            //            else if (args_info.model_arg==2){
            //               prefact = par.alpha*rho_fct2pp(pcell)/locdist;
            //            } else if (args_info.model_arg==3){
            //               prefact = par.alpha*rho_fctpp(pcell)/locdist;
            //            }
            dx += (vertx - vnext.x)*prefact;
            dy += (verty - vnext.y)*prefact;

            locdist=dist(&vprev,&v);
            if (args_info.model_arg==1){
               prefact = par.alpha*c.rho_fctpp()/locdist;
            } else if (args_info.model_arg==2){
               prefact = par.alpha*c.rho_fctpp2()/locdist;
            }
            //            else if (args_info.model_arg==2){
            //               prefact = ((double *)params)[0]*rho_fct2pp(pcell)/locdist;
            //            } else if (args_info.model_arg==3){
            //               prefact = ((double *)params)[0]*rho_fctpp(pcell)/locdist;
            //            }
            dx += (vertx - vprev.x)*prefact;
            dy += (verty - vprev.y)*prefact;

            /* Displacement along x and y due to area extension */
            if (args_info.model_arg==1){
               prefact = c.vol_fctpa()+par.alpha*c.rho_fctpa();
            } else if (args_info.model_arg==2){
               prefact = c.vol_fctpa()+par.alpha*c.rho_fctpa2();
            }
            //            else if (args_info.model_arg==2){
            //               prefact = area_fct2p(pcell)+((double *)params)[0]*rho_fct2pa(pcell);
            //            } else if (args_info.model_arg==3){
            //               prefact = area_fct2p(pcell)+((double *)params)[0]*rho_fctpa(pcell);
            //            }
            dxa += prefact*0.5*(vnext.y - vprev.y); 
            dya += prefact*0.5*(vprev.x - vnext.x);


            /* Random force */
            //            if (args_info.rand_force_flag){
            //               dxd+=-pcell->randfx/pcell->nb_vertices;
            //               dyd+=-pcell->randfy/pcell->nb_vertices;
            //            }

         }
         if (&c!=&c0){
            double el1=c.ellipticity;
            double angle1=c.angle;
//            double ux1=cos(angle1);
            double ux1=1.0-angle1*angle1/2;
            double uy1=angle1-angle1*angle1*angle1/6.0;
//            double uy1=sin(angle1);
            double scal1=ux1*(vertx-vnext.x)+uy1*(verty-vnext.y);
            dxm+=par.pol_force*el1*scal1*ux1;
            dym+=par.pol_force*el1*scal1*uy1;
         }
         if (&cellp!=&c0){
            double el2=cellp.ellipticity;
            double angle2=cellp.angle;
//            double ux2=cos(angle2);
            double ux2=1.0-angle2*angle2/2;
            double uy2=angle2-angle2*angle2*angle2/6.0;
//            double uy2=sin(angle2);
            double scal2=ux2*(vertx-vnext.x)+uy2*(verty-vnext.y);
            dxm+=par.pol_force*el2*scal2*ux2;
            dym+=par.pol_force*el2*scal2*uy2;
         }
      }

      /*      if (args_info.rand_force_flag){
              gsl_vector_set(df, 2*i, dx+dxa+dxd);
              gsl_vector_set(df, 2*i+1, dy+dya+dyd);
              } else {*/
      gsl_vector_set(df, 2*i, dx+dxa+dxm);
      gsl_vector_set(df, 2*i+1, dy+dya+dym);
      /*        }*/
      i++;
   }
   gsl_vector_set(df,2*lat.verts.size(),0.0);
   i=0;
   for (cicell ic=lat.cells.begin()+1;ic!=lat.cells.end();++ic){
      Cell & c=**ic;
//      if (args_info.model_arg==1){
         gsl_vector_set(df, 2*lat.verts.size()+i,
               c.vol_fctpe()+par.beta*c.thick_1_fctp()+2.0*par.gamma*c.thick_fctp()
               );
//      }
      i++;
   }
}

   void
enerdener3d(const gsl_vector *x, void *params, double *f, gsl_vector *df)
{
   /* Compute both the gradient and the energy */
   /* Compute the gradient of the energy */
   Simulation & simu=*(Simulation *)params;
   Lattice & lat=simu.lattice;
   Parameters & par=simu.params;

   unsigned int i=0;
   for (civert iv=lat.verts.begin();iv!=lat.verts.end();iv++){
      (*iv)->x=gsl_vector_get(x,2*i);
      (*iv)->y=gsl_vector_get(x,2*i+1);
      i++;
   }
   i=0;
   for (cicell iv=lat.cells.begin()+1;iv!=lat.cells.end();++iv){
      (*iv)->thick=gsl_vector_get(x,2*lat.verts.size()+i);
      i++;
   }
   lat.update_pa();

   Cell & c0=**lat.cells.begin();
   i=0;
   for (civert iv=lat.verts.begin();iv!=lat.verts.end();iv++){
      Vertex & v=**iv;
      double vertx=v.x;
      double verty=v.y;
      double dx,dy,dxa,dya,dxd,dyd,dxm,dym;
      dx=dy=dxa=dya=dxd=dyd=dxm=dym=0;
      for (unsigned int j=0; j<3; j++){
         Cell & c    = *v.pncell[j];
         Cell & cellp    = *v.pncell[(j+2)%3];
         Vertex & vnext = *v.pneighb[j];
         Vertex & vprev = *v.pneighb[(j+1)%3];

         double prefact=0;

         /* Displacement along x and y due to perimeter extension */
         if (&c!=&c0){
            double locdist=dist(&vnext,&v);
            if (args_info.model_arg==1){
               prefact = par.alpha*c.rho_fctpp()/locdist;
            } else if (args_info.model_arg==2){
               prefact = par.alpha*c.rho_fctpp2()/locdist;
            }
            //            else if (args_info.model_arg==2){
            //               prefact = par.alpha*rho_fct2pp(pcell)/locdist;
            //            } else if (args_info.model_arg==3){
            //               prefact = par.alpha*rho_fctpp(pcell)/locdist;
            //            }
            dx += (vertx - vnext.x)*prefact;
            dy += (verty - vnext.y)*prefact;

            locdist=dist(&vprev,&v);
            if (args_info.model_arg==1){
               prefact = par.alpha*c.rho_fctpp()/locdist;
            } else if (args_info.model_arg==2){
               prefact = par.alpha*c.rho_fctpp2()/locdist;
            }
            //            else if (args_info.model_arg==2){
            //               prefact = ((double *)params)[0]*rho_fct2pp(pcell)/locdist;
            //            } else if (args_info.model_arg==3){
            //               prefact = ((double *)params)[0]*rho_fctpp(pcell)/locdist;
            //            }
            dx += (vertx - vprev.x)*prefact;
            dy += (verty - vprev.y)*prefact;

            /* Displacement along x and y due to area extension */
            if (args_info.model_arg==1){
               prefact = c.vol_fctpa()+par.alpha*c.rho_fctpa();
            } else if (args_info.model_arg==2){
               prefact = c.vol_fctpa()+par.alpha*c.rho_fctpa2();
            }
            //            else if (args_info.model_arg==2){
            //               prefact = area_fct2p(pcell)+((double *)params)[0]*rho_fct2pa(pcell);
            //            } else if (args_info.model_arg==3){
            //               prefact = area_fct2p(pcell)+((double *)params)[0]*rho_fctpa(pcell);
            //            }
            dxa += prefact*0.5*(vnext.y - vprev.y); 
            dya += prefact*0.5*(vprev.x - vnext.x);

            /* Random force */
            //            if (args_info.rand_force_flag){
            //               dxd+=-pcell->randfx/pcell->nb_vertices;
            //               dyd+=-pcell->randfy/pcell->nb_vertices;
            //            }

         }
         if (&c!=&c0){
            double el1=c.ellipticity;
            double angle1=c.angle;
//            double ux1=cos(angle1);
            double ux1=1.0-angle1*angle1/2;
            double uy1=angle1-angle1*angle1*angle1/6.0;
//            double uy1=sin(angle1);
            double scal1=ux1*(vertx-vnext.x)+uy1*(verty-vnext.y);
            dxm+=par.pol_force*el1*scal1*ux1;
            dym+=par.pol_force*el1*scal1*uy1;
         }
         if (&cellp!=&c0){
            double el2=cellp.ellipticity;
            double angle2=cellp.angle;
//            double ux2=cos(angle2);
            double ux2=1.0-angle2*angle2/2;
            double uy2=angle2-angle2*angle2*angle2/6.0;
//            double uy2=sin(angle2);
            double scal2=ux2*(vertx-vnext.x)+uy2*(verty-vnext.y);
            dxm+=par.pol_force*el2*scal2*ux2;
            dym+=par.pol_force*el2*scal2*uy2;
         }
      }

      /*      if (args_info.rand_force_flag){
              gsl_vector_set(df, 2*i, dx+dxa+dxd);
              gsl_vector_set(df, 2*i+1, dy+dya+dyd);
              } else {*/
      gsl_vector_set(df, 2*i, dx+dxa+dxm);
      gsl_vector_set(df, 2*i+1, dy+dya+dym);
      /*        }*/
      i++;
   }
   gsl_vector_set(df,2*lat.verts.size(),0.0);
   i=0;
   for (cicell ic=lat.cells.begin()+1;ic!=lat.cells.end();++ic){
      Cell & c=**ic;
      if (args_info.model_arg==1){
         gsl_vector_set(df, 2*lat.verts.size()+i,
               c.vol_fctpe()+par.beta*c.thick_1_fctp()+2.0*par.gamma*c.thick_fctp()
               );
      } else if (args_info.model_arg==2){
         gsl_vector_set(df, 2*lat.verts.size()+i,
               c.vol_fctpe()+par.beta*c.thick_1_fctp()+2.0*par.gamma*c.thick_fctp()
               );
      }
      i++;
   }
   /*
      gsl_vector_set(df, 2*nb_vertex_tot, 0.0);
      for (i=1;i<nb_cell_tot;i++){
      pcell=web_dual[i];
      gsl_vector_set(df, 2*nb_vertex_tot+i, ((double *)params)[1]*thick_fctp(pcell));
      }*/
   double h=0;
   for (cicell ic=lat.cells.begin()+1;ic!=lat.cells.end();++ic){
      Cell & c=**ic;
      if (args_info.model_arg==1){
         h+=c.vol_fct()+par.alpha*c.rho_fct();
      } else if (args_info.model_arg==2){
         h+=c.vol_fct()+par.alpha*c.rho_fct2();
      } 
//      else if (args_info.model_arg==2){
//         h+=area_fct2(pcell)+((double *)params)[0]*rho_fct2(pcell);
//      } else if (args_info.model_arg==3){
//         h+=area_fct2(pcell)+((double *)params)[0]*rho_fct(pcell);
//      }
//      if (args_info.rand_force_flag){
//         h+=drift_fct(pcell);
//      }
      h+=par.beta*c.thick_1_fct();
      h+=par.gamma*c.thick_fct();

      h+=par.pol_force*c.myopol_en();
   }
   *f=h;
}



