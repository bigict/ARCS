#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([ARCS], [0.91], [BUG-REPORT-ADDRESS])
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE
# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_LIBTOOL 

PKG_CHECK_MODULES([LOG4CXX], [liblog4cxx>=0.10.0])

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lpthread':
AC_CHECK_LIB([pthread], [main])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_CHECK_FUNCS([memset pow sqrt])

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/contiging/Makefile
                 src/contiging/32-mer/Makefile
                 src/contiging/64-mer/Makefile
                 src/contiging/96-mer/Makefile
                 src/gap_filling/Makefile
                 src/remove_repeats/Makefile
                 src/scaffolding/Makefile
                 src/scoring/Makefile])
AC_CONFIG_SUBDIRS([src/contiging/32-mer/Kmer_Hash_md5
                   src/contiging/32-mer/Kmer_Hash_old])
AC_OUTPUT
