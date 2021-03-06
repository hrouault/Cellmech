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

#ifndef OPTI_HPP
#define OPTI_HPP 1


double energy(const gsl_vector *x,void *params);
void denergy(const gsl_vector *x, void *params, gsl_vector *df);
void enerdener(const gsl_vector *x, void *params, double *f, gsl_vector *df);
double energy3d(const gsl_vector *x,void *params);
void denergy3d(const gsl_vector *x, void *params, gsl_vector *df);
void enerdener3d(const gsl_vector *x, void *params, double *f, gsl_vector *df);

#endif /* OPTI_HPP */
