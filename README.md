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

# references

This software is open source (See licensing details elsewhere). In addition to the formal licensing terms, we would greatly appreciate an acknowledgement. Preferably in the form of a citation and a link to the web-page.

Citation: Reerink et al. ([2010](https://www.geosci-model-dev.net/3/13/2010/),[2016](https://www.geosci-model-dev.net/9/4111/2016/)):

Reerink, T. J., M. A. Kliphuis, R. S. W. van de Wal (2010): Mapping technique of climate fields between GCM’s and ice models, Geoscientific Model Development, 3, 13–41, doi:10.5194/gmd-3-13-2010, 
https://www.geosci-model-dev.net/3/13/2010/

Reerink, T. J., W. J. van de Berg, and R. S. W. van de Wal (2016), Oblimap 2.0: a fast climate model–ice sheet model coupler including online embeddable mapping routines, Geoscientific Model Development, 9 (11), 4111–4132, doi:10.5194/gmd-9-4111-2016, https://www.geosci-model-dev.net/9/4111/2016/


The OBLIMAP User Guide can be cited as follows:

Reerink, T. J.: [OBLIMAP User Guide, version 1.0](https://github.com/oblimap/oblimap/blob/master/documentation/oblimap-user-guide.pdf), accompanying OBLIMAP 2.0, Tech. rep., Institute for Marine and Atmospheric research Utrecht, Utrecht University, 3508 TA Utrecht, The Netherlands, 2016.
