#***********************************************************************
# Makefile   
#************ ***********************************************************
# This makefile may be used to create an executable 
# for the Hydro_NetCDF library
# 
FC = gfortran
#FFLAGS = -O
FFLAGS = -g -ffree-line-length-none
AD2CDF = adcirc2netcdf_v5


CDFINCL = -I$(NETCDF_INCDIR)
CDFLIBS = -L$(NETCDF_LIBDIR) -lnetcdff 
#CDFLIBS = -L$(NETCDF_LIBDIR) -lnetcdff -lnetcdf

SRCS = adcirc2netcdf_v5.f90

OBJS = ${SRCS:.f,.f90=.o}

$(AD2CDF):	$(OBJS)
		 $(FC) -o $(AD2CDF) $(FFLAGS) $(OBJS) $(CDFINCL) $(CDFLIBS) 
