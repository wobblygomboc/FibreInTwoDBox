# Name of executables: The prefix "noinst" means the executables don't
# have to be installed anywhere.
noinst_PROGRAMS= \
fibre_in_two_d_box

#---------------------------------------------------------------------------

# Local sources that Joe's Poisson code depends on:
fibre_in_two_d_box_SOURCES = fibre_in_two_d_box.cc unstructured_moffatt_mesh_3_inner_boundaries.h moffat_solution.h

# Required libraries: 
fibre_in_two_d_box_LDADD = -L@libdir@  -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
#---------------------------------------------------------------------------

# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
# Automake will replace the variable @includedir@ with the actual
# include directory.
INCLUDES = -I@includedir@  

#EXTRA_DIST = boundaries_2d.lay boundaries_2d.png convergence_2d.lay run2d.bash



