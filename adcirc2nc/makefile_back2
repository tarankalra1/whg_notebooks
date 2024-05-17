#***********************************************************************
# Makefile   
#************ ***********************************************************
# This makefile may be used to create an executable 
# for the Hydro_NetCDF library
# 
FC = ifort
#FFLAGS = -O
FFLAGS = -g
AD2CDF = fort54_to_nc


CDFINCL = -I$(NETCDF_INCDIR)
CDFLIBS = -L$(NETCDF_LIBDIR) -lnetcdff 
#CDFLIBS = -L$(NETCDF_LIBDIR) -lnetcdff -lnetcdf

SRCS = tidescdf_fem.f fort54_38_to_nc.f90

OBJS = ${SRCS:.f,.f90=.o}

$(AD2CDF):	$(OBJS)
		 $(FC) -o $(AD2CDF) $(FFLAGS) $(OBJS) $(CDFINCL) $(CDFLIBS) 
