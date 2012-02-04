# 
# Copyright (C) 2004-2012 Hervé Rouault <herve.rouault@pasteur.fr>
# 
# This file is part of Cellmech.
# 
# Cellmech is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# Cellmech is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with Cellmech. If not, see <http://www.gnu.org/licenses/>.
# 

## Process this file with automake to produce Makefile.in

bin_PROGRAMS = cellmech

BUILT_SOURCES = cmdline.c

cellmech_SOURCES = main.cpp main.hpp lattice.cpp lattice.hpp opti.cpp opti.hpp vectortypes.cpp vectortypes.hpp cmdline.c cmdline.h

cmdline.c: parameters.ggo
	gengetopt --input=$(top_srcdir)/src/parameters.ggo
