AC_PREREQ([2.59])

# Initialize with name, version, and support email address
AC_INIT([FMS], [2.0-development], [])

# Find out about the host we're building on.
AC_CANONICAL_HOST

# Find out about the target we're building for
AC_CANONICAL_TARGET

# Keep libtool macros in a m4 directory
AC_CONFIG_MACRO_DIR([m4])

# Find the Fortran compiler
AX_PROG_FC_MPI([], [found_mpi="yes"], [found_mpi="no"])
if test "x$found_mpi" = "xno"; then
  AC_MSG_FAILURE([Unable to find an MPI compatible compiler.])
fi
echo "FCFLAGS = $FCFLAGS"
AC_LANG_PREPROC

# Find the install program
AC_PROG_INSTALL

# Check to see if any macros must be set to enable large (>2GB) files.
AC_SYS_LARGEFILE

# Check fo libFMS.  The Fortran modules for libFMS are required
AC_LANG_PUSH(Fortran)
AC_ARG_WITH([libFMS], AS_HELP_STRING([--with-libFMS=<dir>], [use libFMS and headers installed in <dir>.  The libFMS headers/modules are required]))
if test "x$with_libFMS" = "xno"; then
  AC_MSG_WARN([You specified --with-libFMS=no.  The headers/modules for libFMS are required.  Still checking, using compiler's default search paths.])
fi

# Check for libFMS library headers/modules
AC_MSG_CHECKING([for libFMS])
AC_LANG_CONFTEST([AC_LANG_PROGRAM([],[[
       use mpp_mod
       call mpp_init]])])
# If user specified location with --with-libFMS, check that path for the headers
if test -n "$with_libFMS" && test "x$with_libFMS" != "xno"; then
  doit='$FC -o conftest -I${with_libFMS}/include conftest.f -L${with_libFMS}/lib -lFMS'
  if AC_TRY_EVAL(doit); then
    AC_SUBST([LIBFMS_LDFLAGS], [-Lsomething -lwhat])
    AC_SUBST([LIBFMS_CPPFLAGS], [-I${with_libFMS}/include])
    libFMS_found=yes
  else
    libFMS_found=no
  fi
else
  doit='$FC -o conftest ${CPPFLAGS} ${FCFLAGS} conftest.f $LDFLAGS -lFMS'
  if AC_TRY_EVAL(doit); then
    libFMS_found=yes
  else
    libFMS_found=no
  fi
fi
AC_MSG_RESULT($libFMS_found)
if test "x$libFMS_found" = "xno"; then
  AC_MSG_ERROR([Unable to find required libFMS library.])
fi

# Set up libtool
LT_PREREQ([2.4])
LT_INIT()
LT_LANG([Fortran])

# Set up automake
AM_INIT_AUTOMAKE([foreign dist-zip])

# These files will be created when the configure script is run
AC_CONFIG_FILES([Makefile
                 t/Makefile
                ])
AC_OUTPUT()
