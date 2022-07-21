from __future__ import division, print_function; __metaclass__ = type
import os, sys, math, itertools
import numpy
from westpa.core.systems import WESTSystem
from westpa.core.binning import RectilinearBinMapper

import logging
log = logging.getLogger(__name__)
log.debug('loading module %r' % __name__)


def npy_coord_loader(fieldname, coord_file, segment, single_point=False):
    print('using this loader')
    coords = numpy.load(coord_file)

    #npts = len(coord_raw)
    #assert coord_raw.shape[1] % 3 == 0
    #ngrps = coord_raw.shape[1] // 3

    #coords = numpy.empty((ngrps, npts, 3), numpy.float32)
    #for igroup in range(ngrps):
    #    for idim in range(3):
    #        coords[igroup,:,idim] = coord_raw[:,igroup*3+idim]
    # convert to Angstroms
    #coords *= 10

    print(coords.shape)
    segment.data[fieldname] = coords
