# automake file for land_null tests

# Find the *.mod files
AM_CPPFLAGS = -I$(top_builddir)
AM_CPPFLAGS += $(LIBFMS_CPPFLAGS)

AM_LDFLAGS = $(LIBFMS_LDFLAGS)

# Test programs to build
check_PROGRAMS = land_model_init_test

# Source code for programs
## land_moddel_init_test
land_model_init_test_SOURCES = land_model_init_test.f90
land_model_init_test_CPPFLAGS = -I$(top_builddir) $(LIBFMS_CPPFLAGS)

# Run programs
TESTS = land_model_init_test
