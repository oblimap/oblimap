# oblimap
OBLIMAP: a fast coupler for coupling ice sheet models with climate models

# abstract
OBLIMAP [Reerink et al., [2010](https://www.geosci-model-dev.net/3/13/2010/gmd-3-13-2010.html), 
[2016](https://www.geosci-model-dev.net/9/4111/2016/gmd-9-4111-2016.html)] is a mapping technique for exchanging climate 
fields between a Global Circulation Model (GCM) and an Ice Sheet Model (ISM).

The output fields of a GCM can be mapped to an ISM with the OBLIMAP package. Mapping in the reverse direction, from an ISM 
to a GCM, is possible as well with the OBLIMAP package. As such it acts as a coupler between grids which are based on 
geographical coordinates and grids which represent a flat surface. The OBLIMAP mapping routines can be used either stand-alone
and embedded, the latter case enables online coupling of an ISM to the GCM components of an Earth System Model (ESM).

# oblimap repository

Getting the oblimap repositry:
```
git clone https://github.com/oblimap/oblimap.git
```
See its README for compiling and immediate running of an default example. OBLIMAP is written in modern Fortran and uses netcdf.
