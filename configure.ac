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
dnl Process this file with autoconf to produce a configure script.

AC_INIT(Cellmech, m4_esyscmd([tools/git-version-gen .tarball-version]), [herve.rouault@pasteur.fr])
AM_INIT_AUTOMAKE()
AC_CONFIG_SRCDIR([src/main.cpp])
AC_CONFIG_HEADER([config.h])
AC_CONFIG_MACRO_DIR([m4])
m4_include([m4/ax_blas.m4])

# Checks for programs.
AC_PROG_CXX
AM_PROG_CC_C_O
AC_PROG_MAKE_SET
AC_PROG_RANLIB

# Checks for libraries.
AC_CHECK_LIB([m],[cos],,[AC_MSG_ERROR([math functions required.])])

AC_CHECK_LIB(pthread,pthread_create,,[AC_MSG_ERROR([pthread library required.])])

##AX_BLAS( [ AC_MSG_NOTICE( [ blas libs variable : x$BLAS_LIBS BLAS_LIBS ] ) ], )
AX_BLAS(LIBS="$BLAS_LIBS $LIBS $FLIBS",AC_CHECK_LIB([gslcblas],[cblas_dgemm],,AC_MSG_ERROR([blas library required.])))

#AX_BLAS([],[AC_MSG_ERROR([blas library required.])])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm],,[AC_MSG_ERROR([libgsl required.])])

# Checks for header files.
AC_CHECK_HEADERS([gsl/gsl_rng.h],[],[AC_MSG_ERROR([gsl headers required.])])


AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT()
