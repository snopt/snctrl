#-------------------------------------------------------------------------------
AC_PREREQ(2.69)
AC_INIT([SNCTRL],[5.1],[optimizers@ccom.ucsd.edu],[snctrl])

AC_CONFIG_SRCDIR([src/control.f90])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([config/m4])

VERSION="5";
AC_SUBST([VERSION])

AC_PREFIX_DEFAULT([`pwd`])

user_FCFLAGS=${FCFLAGS}
user_CFLAGS=${CFLAGS}
user_CPPFLAGS=${CPPFLAGS}
user_CXXFLAGS=${CXXFLAGS}
user_LDFLAGS=${LDFLAGS}


# Libtool
LT_INIT
AC_SUBST([LIBTOOL_DEPS])


# Check compilers.
AC_PROG_INSTALL
AC_PROG_FC([gfortran nagfor ifort])
AC_PROG_CC
AC_PROG_CXX

AC_FC_LIBRARY_LDFLAGS
AC_SUBST([FCLIBS])

AC_FC_MODULE_OUTPUT_FLAG([],
  [ case $FC in
      *nagfor* ) FC_MODOUT="-mdir ";;
    esac
  ])
AC_FC_MODULE_FLAG([],
  [ case $FC in
      *nagfor* ) FC_MODINC="-mdir ";;
    esac
  ])
AC_SUBST(FC_MODOUT)
AC_SUBST(FC_MODINC)

AC_PROG_MAKE_SET
AC_SUBST([SET_MAKE])


dnl ----------------------------------------------------------------------------

# BLAS library.
AC_ARG_WITH([blas],
	    [AS_HELP_STRING([--with-blas=<lib>], [use BLAS library <lib>])],
	    [with_blas=${withval} ], [ blas=no ] )
case $with_blas in
	yes | "" | no) blas=no ;;
	-* | */* | *.a | *.so | *.so.* | *.o)
	  blas=yes; BLAS_LIBS="$with_blas" ;;
	*)
	  blas=yes; BLAS_LIBS="-l$with_blas" ;;
esac
AC_SUBST([blas])
AC_SUBST([BLAS_LIBS])

dnl ----------------------------------------------------------------------------

#AC_ARG_WITH([matlab],
#	    [AS_HELP_STRING([--with-matlab=DIR],[Provide location of MATLAB]) ],
#            [with_matlab=yes; MATLABDIR=${withval};],
#            [with_matlab=no; matlab="no"; MATLABDIR="";] )
with_matlab=no;
matlab=no;
MEX="";
if test x"${with_matlab}" == xyes; then
  AC_MSG_CHECKING([for Matlab])

  if test -e ${MATLABDIR}/bin/mex; then
    with_matlab=yes;
    matlab=yes;
    MEX="${MATLABDIR}/bin/mex";
    AC_MSG_RESULT([yes])
    AC_MSG_WARN([Please make sure to use a compiler supported by Matlab]);
  else
    with_matlab=no;
    matlab=no;
    MEX="";
    AC_MSG_RESULT([mex not found in ${MATLABDIR}/bin])
  fi
fi
AC_SUBST(matlab)
AC_SUBST(MATLABDIR)
AC_SUBST(MEX)

dnl ----------------------------------------------------------------------------

# Checking default integer size
prec=32
AC_LANG_PUSH([Fortran])
AC_CACHE_CHECK([for default size of Fortran integer], [my_cv_intsize],[
AC_RUN_IFELSE([
  AC_LANG_PROGRAM([],
     [
        use iso_fortran_env, only : file_storage_size
        integer :: nbytes, nunits
        inquire(iolength=nunits) 10
        nbytes = nunits*file_storage_size/8
        open(unit=42,file='conftestval',status='new')
        write(42,'(i1)') nbytes
        close(42)
     ])
  ],
  [ my_cv_intsize=`cat conftestval`
    if test x"${my_cv_intsize}" == x8; then
      prec=64
    fi
  ],
  [ AC_MSG_WARN([could not determine integer size -- continuing with 32-bit integers])
    prec=32 ])
])
AC_LANG_POP([Fortran])


dnl ----------------------------------------------------------------------------

# User-defined option will override above test
AC_ARG_WITH([32],
            [AS_HELP_STRING([--with-32],
                            [32-bit integer support (default is 32-bit)]) ],
 	    [prec=32;
             AC_MSG_NOTICE([32-bit integers enabled]) ], [])

AC_ARG_WITH([64],
            [AS_HELP_STRING([--with-64],
                            [64-bit integer support (default is 32-bit)]) ],
 	    [prec=64;
             AC_MSG_NOTICE([64-bit integers enabled]) ], [])

#AC_ARG_WITH([128],
#            [AS_HELP_STRING([--with-128],
#                            [quad-precision support (64-bit ints, 128-bit reals)]) ],
# 	    [prec=128;
#             AC_MSG_NOTICE([quad-precision enabled (64-bit ints, 128-bit reals)]) ], [])

AC_ARG_WITH([debug],
            [AS_HELP_STRING([--with-debug],
                            [Compile with debug flags])],
 	    [DBGFLAG="-g"; OPTFLAG="";],
	    [OPTFLAG="-O"; DBGFLAG="";])

dnl ----------------------------------------------------------------------------

# SNOPT library.
AC_ARG_WITH([snopt],
	    [AS_HELP_STRING([--with-snopt=/path/to/snopt], [Provide location of SNOPT])],
	    [with_snopt=yes; SNOPT_LIB=${withval}; ], [ with_snopt=no; SNOPT_LIB=""; ] )
if test x"${with_snopt}" == xyes; then
  sav_LDFLAGS=${LDFLAGS}
  LDFLAGS="-L${SNOPT_LIB} ${LDFLAGS}"
  AC_CHECK_LIB([snopt7],[main], [],
	       [SNOPT_LIB=""; AC_MSG_WARN([--with-snopt was given, but test for SNOPT failed])])
  LDFLAGS=${sav_LDFLAGS}
fi
AC_SUBST([SNOPT_LIB])

dnl ----------------------------------------------------------------------------

FCFLAGS="${user_FCFLAGS}"
CFLAGS="${user_CFLAGS}"
CXXFLAGS="${user_CXXFLAGS}"
LDFLAGS="${user_LDFLAGS}"

if test x"${prec}" = x64; then
   PC_FLAGS="-DD64"

   case "${FC}" in
   *gfortran*)
     PF_FLAGS="-fdefault-integer-8"
   ;;
   *nagfor*)
     PF_FLAGS="-i8"
   ;;
   *ifort*)
     PF_FLAGS="-i8"
   ;;
   esac
elif test x"${prec}" = x128; then
   PC_FLAGS=" -DD128"

   case "${FC}" in
   *gfortran*)
     PF_FLAGS=" -fdefault-integer-8"
   ;;
   *nagfor* )
     PF_FLAGS=" -i8"
   ;;
   *ifort* )
     PF_FLAGS=" -i8"
   ;;
   esac
fi

case "${FC}" in
  *nagfor* )
    PF_FLAGS+=" -w"
  ;;
esac

if test -z "${FCFLAGS}"; then
  if test -z "${DBGFLAG}"; then
    FCFLAGS="${OPTFLAG} ${PF_FLAGS}"
  else
    FCFLAGS="${DBGFLAG} ${PF_FLAGS}"
  fi
else
  FCFLAGS="${FCFLAGS} ${PF_FLAGS}"
fi

if test -z "${CFLAGS}"; then
  if test -z "${DBGFLAG}"; then
    CFLAGS="${OPTFLAG} ${PC_FLAGS}"
  else
    CFLAGS="${DBGFLAG} ${PC_FLAGS}"
  fi
else
  CFLAGS="${CFLAGS} ${PC_FLAGS}"
fi

if test -z "${CXXFLAGS}"; then
  if test -z "${DBGFLAG}"; then
    CXXFLAGS="${OPTFLAG}"
  else
    CXXFLAGS="${DBGFLAG}"
  fi
fi

if test -z "${MEXFLAGS}"; then
  if test -z "${DBGFLAG}"; then
    MEXFLAGS="${OPTFLAG}"
  else
    MEXFLAGS="${DBGFLAG}"
  fi
fi

AC_SUBST([prec])
AC_SUBST([FCFLAGS])
AC_SUBST([CFLAGS])
AC_SUBST([CXXFLAGS])
AC_SUBST([MEXFLAGS])

dnl ----------------------------------------------------------------------------

# Set modules.
if test x"${f2c_enabled}" == xno; then
  possible_modules="src examples ${api_modules}"
else
  possible_modules="src examples ${c_modules} ${api_modules}"
fi
modules=""
for module in ${possible_modules}; do
	if test -d ${srcdir}/${module}; then
		modules="${modules} ${module}"
	fi
done

if test -d ${srcdir}/ampl; then
  if test x"${with_ampl}" == xyes; then
    modules="${modules} ampl"
  fi
fi

if test -d ${srcdir}/mex; then
  if test x"${with_matlab}" == xyes; then
    modules="${modules} mex"
  fi
fi

if test -d ${srcdir}/cmex; then
  if test x"${with_cmex}" == xyes; then
    modules="${modules} cmex"
  fi
fi

AC_SUBST([modules])

dnl ----------------------------------------------------------------------------
case $FC in
  *nagfor* ) enable_shared=no;
esac

AC_SUBST([enable_shared])
AC_SUBST([enable_static])

AC_OUTPUT([GNUmakefile])

echo ""
echo "SNCTRL will compile the following modules:"
echo "    ${modules}"
echo ""
echo "  Using external BLAS library?      ${blas}    ${BLAS_LIBS} "
echo ""
#echo "  Matlab mex-executable available?  ${matlab}  ${MATLABDIR} "
#if test x"${matlab}" == xyes; then
#  echo "    (Make sure to use a compiler compatible with Matlab to create mex-files)"
#  echo ""
#fi
echo ""
echo "Libraries will be installed into ${prefix}/lib"
echo "  Shared libraries will be created: ${enable_shared}"
echo "  Static libraries will be created: ${enable_static}"
echo ""
echo ""
echo "To build all libraries              'make'          "
echo "To build the SNCTRL Fortran library 'make snctrl'    "
echo "To check the build                  'make check'    "
echo "To install the libraries            'make install'  "
#echo "To build the Matlab mex-executable  'make matlab'   "
echo "To build the examples               'make examples' "
echo ""
