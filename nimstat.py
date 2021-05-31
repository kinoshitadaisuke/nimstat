#!/usr/pkg/bin/python3.9

#
# Time-stamp: <2021/05/31 10:26:37 (CST) daisuke>
#

#
# nimstat.py
#
#   a clone of IRAF's imstatistics task
#
#   author: Kinoshita Daisuke
#
#   version 1.0 on 31/May/2021 by Kinoshita Daisuke
#

#
# Usage:
#
#   % ./nimstat.py -r sigclip -t 3.0 -n 30 *.fits
#

# importing argparse module
import argparse

# importing pathlib module
import pathlib

# importing numpy module
import numpy
import numpy.ma

# importing astropy module
import astropy
import astropy.io.fits

# construction pf parser object
desc = 'calculating statistical information of FITS files'
parser = argparse.ArgumentParser (description=desc)

# adding arguments
list_rejection = ['NONE', 'sigclip']
parser.add_argument ('-r', '--rejection', default='NONE', \
                     choices=list_rejection, help='rejection algorithm')
parser.add_argument ('-t', '--threshold', type=float, default=4.0, \
                     help='threshold for sigma clipping (default: 4)')
parser.add_argument ('-n', '--maxiters', type=int, default=10, \
                     help='maximum number of iterations (default: 10)')
parser.add_argument ('files', nargs='+', help='FITS files')

# command-line argument analysis
args = parser.parse_args ()

# input parameters
rejection  = args.rejection
threshold  = args.threshold
maxiters   = args.maxiters
list_files = args.files

# function to read a FITS file
def read_fits (file_fits):
    # reading FITS file
    with astropy.io.fits.open (file_fits) as hdu:
        # reading header and image
        header = hdu[0].header
        image  = hdu[0].data
        # if no image in PrimaryHDU, then read next HDU
        if (header['NAXIS'] == 0):
            header = hdu[1].header
            image  = hdu[1].data
    # returning header and image
    return (header, image)

# printing information
print ("# Input parameters")
print ("#   rejection algorithm = %s" % rejection)
print ("#   threshold of sigma-clipping = %f" % threshold)
print ("#   maximum number of iterations = %d" % maxiters)

# printing header
print ("#")
print ("# %-25s %8s %8s %8s %7s %8s %8s" \
       % ('file name', 'npix', 'mean', 'median', 'stddev', 'min', 'max') )
print ("#")

# scanning files
for file_fits in list_files:
    # if the file is not a FITS file, then skip
    if not (file_fits[-5:] == '.fits'):
        continue

    # reading header and image of FITS file
    (header, image) = read_fits (file_fits)

    #
    # calculations of statistical values
    #

    # making a masked array
    image_masked = numpy.ma.array (image, mask=False)

    # sigma-clipping algorithm
    if (rejection == 'sigclip'):
        image_masked = numpy.ma.array (image, mask=False)
        # iterations
        for j in range (maxiters):
            # number of usable pixels of previous iterations
            npix_prev = len (numpy.ma.compressed (image_masked) )
            # calculation of median
            median = numpy.ma.median (image_masked)
            # calculation of standard deviation
            stddev = numpy.ma.std (image_masked)
            # lower threshold
            low = median - threshold * stddev
            # higher threshold
            high = median + threshold * stddev
            # making a mask
            mask = (image_masked < low) | (image_masked > high)
            # making a masked array
            image_masked = numpy.ma.array (image_masked, mask=mask)
            # number of usable pixels
            npix_now = len (numpy.ma.compressed (image_masked) )
            # leaving the loop, if number of usable pixels does not change
            if (npix_now == npix_prev):
                break
        
    # calculation of mean, median, stddev, min, and max
    mean   = numpy.ma.mean (image_masked)
    median = numpy.ma.median (image_masked)
    stddev = numpy.ma.std (image_masked)
    vmin   = numpy.ma.min (image_masked)
    vmax   = numpy.ma.max (image_masked)

    # number of pixels
    npix = len (image_masked.compressed () )

    # file name
    path_file = pathlib.Path (file_fits)
    filename = path_file.name

    # printing result
    print ("%-27s %8d %8.2f %8.2f %7.2f %8.2f %8.2f" \
           % (filename, npix, mean, median, stddev, vmin, vmax) )
