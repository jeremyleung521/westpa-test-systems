from netCDF4 import Dataset
import numpy
rootgrp = Dataset('parent.ncrst', 'r', format="NETCDF3")
rootgrp2 = Dataset('seg.nc', 'r', format="NETCDF3")

coords = numpy.concatenate(([rootgrp.variables['coordinates'][:]], rootgrp2.variables['coordinates'][:]))

numpy.save('coords.npy', coords)
