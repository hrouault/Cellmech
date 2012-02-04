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

#ifndef MAIN_HPP
#define MAIN_HPP 1

/*
 * =====================================================================================
 *
 *       Filename:  main.hpp
 *
 *    Description:  Header of the main file
 *
 *        Version:  1.0
 *        Created:  12.10.2009 14:21:53
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_rng.h>

#include "cmdline.h"

using namespace std;

extern struct gengetopt_args_info args_info;

extern gsl_rng * rng;

typedef vector<pthread_t> vthds;
typedef vthds::iterator ivthds;

#endif /* MAIN_HPP */
