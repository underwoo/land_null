# build the libLandNull_GFDL library
lib_LTLIBRARIES = libland.la

# These linker flags specify libtool version info.
# See http://www.gnu.org/software/libtool/manual/libtool.html#Libtool-versioning
# for information regarding incrementing `-version-info`.
libland_la_LDFLAGS = -version-info 2:0:0

# Add the object files to the library
#libland_la_LIBADD = ${top_builddir}/land_model.la

libland_la_SOURCES = land_model.F90

# Additional FCFLAGS (set by configure to find headers/modules for libFMS
AM_CPPFLAGS = -I${LIBFMS_INCLUDE_DIR}

# Mod file depends on its o file, is build and then installed
land_model_mod.mod: land_model.$(OBJEXT)
BUILT_SOURCES = land_model_mod.mod
include_HEADERS = land_model_mod.mod

CLEANFILES = *.mod

echo_stuff:
	echo \"$(LIBFMS_INCLUDE)\"
