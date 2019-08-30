#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([pygwcoh], [0.0.1], [shallyn.liu@foxmail.com])
AC_CONFIG_SRCDIR([cohCore.c])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE

# Checks for programs.
AC_PROG_CC

# Checks for libraries.
AC_PROG_LIBTOOL

dnl python checks (you can change the required python version bellow)
AM_PATH_PYTHON(3.6)
AC_CHECK_LIB([gsl])

# Checks for header files.
AC_CHECK_HEADERS([stdint.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_INT16_T
AC_TYPE_INT64_T
AC_TYPE_UINT16_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([memset modf])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
