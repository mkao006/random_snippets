########################################################################
# Title: Maps in Python
# Date: 2015-06-15
# Source: http://matplotlib.org/basemap/users/examples.html
########################################################################

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
# set up orthographic map projection with
# perspective of satellite looking down at 50N, 100W.
# use low resolution coastlines.
map = Basemap(projection='ortho', lat_0=45, lon_0=-100, resolution='l')
# draw coastlines, country boundaries, fill continents.
map.drawcoastlines(linewidth=0.25)
map.drawcountries(linewidth=0.25)
map.fillcontinents(color='coral', lake_color='aqua')
# draw the edge of the map projection region (the projection limb)
map.drawmapboundary(fill_color='aqua')
# draw lat/lon grid lines every 30 degrees.
map.drawmeridians(np.arange(0, 360, 30))
map.drawparallels(np.arange(-90, 90, 30))
# make up some data on a regular lat/lon grid.
nlats = 73
nlons = 145
delta = 2. * np.pi / (nlons - 1)
lats = (0.5 * np.pi - delta * np.indices((nlats, nlons))[0, :, :])
lons = (delta * np.indices((nlats, nlons))[1, :, :])
wave = 0.75 * (np.sin(2. * lats)**8 * np.cos(4. * lons))
mean = 0.5 * np.cos(2. * lats) * ((np.sin(2. * lats))**2 + 2.)
# compute native map projection coordinates of lat/lon grid.
x, y = map(lons * 180. / np.pi, lats * 180. / np.pi)
# contour data over the map.
cs = map.contour(x, y, wave + mean, 15, linewidths=1.5)
plt.title('contour lines over filled continent background')
plt.show()

from mpl_toolkits.basemap import Basemap, cm
# requires netcdf4-python (netcdf4-python.googlecode.com)
from netCDF4 import Dataset as NetCDFFile
import numpy as np
import matplotlib.pyplot as plt

# plot rainfall from NWS using special precipitation
# colormap used by the NWS, and included in basemap.

nc = NetCDFFile('../../../examples/nws_precip_conus_20061222.nc')
# data from http://water.weather.gov/precip/
prcpvar = nc.variables['amountofprecip']
data = 0.01 * prcpvar[:]
latcorners = nc.variables['lat'][:]
loncorners = -nc.variables['lon'][:]
lon_0 = -nc.variables['true_lon'].getValue()
lat_0 = nc.variables['true_lat'].getValue()
# create figure and axes instances
fig = plt.figure(figsize=(8, 8))
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
# create polar stereographic Basemap instance.
m = Basemap(projection='stere', lon_0=lon_0, lat_0=90., lat_ts=lat_0,
            llcrnrlat=latcorners[0], urcrnrlat=latcorners[2],
            llcrnrlon=loncorners[0], urcrnrlon=loncorners[2],
            rsphere=6371200., resolution='l', area_thresh=10000)
# draw coastlines, state and country boundaries, edge of map.
m.drawcoastlines()
m.drawstates()
m.drawcountries()
# draw parallels.
parallels = np.arange(0., 90, 10.)
m.drawparallels(parallels, labels=[1, 0, 0, 0], fontsize=10)
# draw meridians
meridians = np.arange(180., 360., 10.)
m.drawmeridians(meridians, labels=[0, 0, 0, 1], fontsize=10)
ny = data.shape[0]
nx = data.shape[1]
lons, lats = m.makegrid(nx, ny)  # get lat/lons of ny by nx evenly space grid.
x, y = m(lons, lats)  # compute map proj coordinates.
# draw filled contours.
clevs = [0, 1, 2.5, 5, 7.5, 10, 15, 20, 30, 40, 50,
         70, 100, 150, 200, 250, 300, 400, 500, 600, 750]
cs = m.contourf(x, y, data, clevs, cmap=cm.s3pcpn)
# add colorbar.
cbar = m.colorbar(cs, location='bottom', pad="5%")
cbar.set_label('mm')
# add title
plt.title(prcpvar.long_name + ' for period ending ' + prcpvar.dateofdata)
plt.show()


"""
plot H's and L's on a sea-level pressure map
(uses scipy.ndimage.filters and netcdf4-python)
"""
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from mpl_toolkits.basemap import Basemap, addcyclic
from scipy.ndimage.filters import minimum_filter, maximum_filter
from netCDF4 import Dataset


def extrema(mat, mode='wrap', window=10):
    """find the indices of local extrema (min and max)
    in the input array."""
    mn = minimum_filter(mat, size=window, mode=mode)
    mx = maximum_filter(mat, size=window, mode=mode)
    # (mat == mx) true if pixel is equal to the local max
    # (mat == mn) true if pixel is equal to the local in
    # Return the indices of the maxima, minima
    return np.nonzero(mat == mn), np.nonzero(mat == mx)

# plot 00 UTC today.
date = datetime.now().strftime('%Y%m%d') + '00'

# open OpenDAP dataset.
# data=Dataset("http://nomads.ncep.noaa.gov:9090/dods/gfs/gfs/%s/gfs_%sz_anl" %\
#        (date[0:8],date[8:10]))
data = Dataset("http://nomads.ncep.noaa.gov:9090/dods/gfs_hd/gfs_hd%s/gfs_hd_%sz" %
               (date[0:8], date[8:10]))


# read lats,lons.
lats = data.variables['lat'][:]
lons1 = data.variables['lon'][:]
nlats = len(lats)
nlons = len(lons1)
# read prmsl, convert to hPa (mb).
prmsl = 0.01 * data.variables['prmslmsl'][0]
# the window parameter controls the number of highs and lows detected.
# (higher value, fewer highs and lows)
local_min, local_max = extrema(prmsl, mode='wrap', window=50)
# create Basemap instance.
m =\
    Basemap(llcrnrlon=0, llcrnrlat=-80, urcrnrlon=360,
            urcrnrlat=80, projection='mill')
# add wrap-around point in longitude.
prmsl, lons = addcyclic(prmsl, lons1)
# contour levels
clevs = np.arange(900, 1100., 5.)
# find x,y of map projection grid.
lons, lats = np.meshgrid(lons, lats)
x, y = m(lons, lats)
# create figure.
fig = plt.figure(figsize=(8, 4.5))
ax = fig.add_axes([0.05, 0.05, 0.9, 0.85])
cs = m.contour(x, y, prmsl, clevs, colors='k', linewidths=1.)
m.drawcoastlines(linewidth=1.25)
m.fillcontinents(color='0.8')
m.drawparallels(np.arange(-80, 81, 20), labels=[1, 1, 0, 0])
m.drawmeridians(np.arange(0, 360, 60), labels=[0, 0, 0, 1])
xlows = x[local_min]
xhighs = x[local_max]
ylows = y[local_min]
yhighs = y[local_max]
lowvals = prmsl[local_min]
highvals = prmsl[local_max]
# plot lows as blue L's, with min pressure value underneath.
xyplotted = []
# don't plot if there is already a L or H within dmin meters.
yoffset = 0.022 * (m.ymax - m.ymin)
dmin = yoffset
for x, y, p in zip(xlows, ylows, lowvals):
    if x < m.xmax and x > m.xmin and y < m.ymax and y > m.ymin:
        dist = [np.sqrt((x - x0)**2 + (y - y0)**2) for x0, y0 in xyplotted]
        if not dist or min(dist) > dmin:
            plt.text(x, y, 'L', fontsize=14, fontweight='bold',
                     ha='center', va='center', color='b')
            plt.text(x, y - yoffset, repr(int(p)), fontsize=9,
                     ha='center', va='top', color='b',
                     bbox=dict(boxstyle="square", ec='None', fc=(1, 1, 1, 0.5)))
            xyplotted.append((x, y))
# plot highs as red H's, with max pressure value underneath.
xyplotted = []
for x, y, p in zip(xhighs, yhighs, highvals):
    if x < m.xmax and x > m.xmin and y < m.ymax and y > m.ymin:
        dist = [np.sqrt((x - x0)**2 + (y - y0)**2) for x0, y0 in xyplotted]
        if not dist or min(dist) > dmin:
            plt.text(x, y, 'H', fontsize=14, fontweight='bold',
                     ha='center', va='center', color='r')
            plt.text(x, y - yoffset, repr(int(p)), fontsize=9,
                     ha='center', va='top', color='r',
                     bbox=dict(boxstyle="square", ec='None', fc=(1, 1, 1, 0.5)))
            xyplotted.append((x, y))
plt.title('Mean Sea-Level Pressure (with Highs and Lows) %s' % date)
plt.show()


from mpl_toolkits.basemap import Basemap, shiftgrid, cm
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

# read in etopo5 topography/bathymetry.
url = 'http://ferret.pmel.noaa.gov/thredds/dodsC/data/PMEL/etopo5.nc'
etopodata = Dataset(url)

topoin = etopodata.variables['ROSE'][:]
lons = etopodata.variables['ETOPO05_X'][:]
lats = etopodata.variables['ETOPO05_Y'][:]
# shift data so lons go from -180 to 180 instead of 20 to 380.
topoin, lons = shiftgrid(180., topoin, lons, start=False)

# plot topography/bathymetry as an image.

# create the figure and axes instances.
fig = plt.figure()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
# setup of basemap ('lcc' = lambert conformal conic).
# use major and minor sphere radii from WGS84 ellipsoid.
m = Basemap(llcrnrlon=-145.5, llcrnrlat=1., urcrnrlon=-2.566, urcrnrlat=46.352,
            rsphere=(6378137.00, 6356752.3142),
            resolution='l', area_thresh=1000., projection='lcc',
            lat_1=50., lon_0=-107., ax=ax)
# transform to nx x ny regularly spaced 5km native projection grid
nx = int((m.xmax - m.xmin) / 5000.) + 1
ny = int((m.ymax - m.ymin) / 5000.) + 1
topodat = m.transform_scalar(topoin, lons, lats, nx, ny)
# plot image over map with imshow.
im = m.imshow(topodat, cm.GMT_haxby)
# draw coastlines and political boundaries.
m.drawcoastlines()
m.drawcountries()
m.drawstates()
# draw parallels and meridians.
# label on left and bottom of map.
parallels = np.arange(0., 80, 20.)
m.drawparallels(parallels, labels=[1, 0, 0, 1])
meridians = np.arange(10., 360., 30.)
m.drawmeridians(meridians, labels=[1, 0, 0, 1])
# add colorbar
cb = m.colorbar(im, "right", size="5%", pad='2%')
ax.set_title('ETOPO5 Topography - Lambert Conformal Conic')
plt.show()

# make a shaded relief plot.

# create new figure, axes instance.
fig = plt.figure()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
# attach new axes image to existing Basemap instance.
m.ax = ax
# create light source object.
from matplotlib.colors import LightSource
ls = LightSource(azdeg=90, altdeg=20)
# convert data to rgb array including shading from light source.
# (must specify color map)
rgb = ls.shade(topodat, cm.GMT_haxby)
im = m.imshow(rgb)
# draw coastlines and political boundaries.
m.drawcoastlines()
m.drawcountries()
m.drawstates()
ax.set_title('Shaded ETOPO5 Topography - Lambert Conformal Conic')
plt.show()

from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset, date2index
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
date = datetime(2007, 12, 15, 0)  # date to plot.
# open dataset.
dataset = \
    Dataset('http://www.ncdc.noaa.gov/thredds/dodsC/OISST-V2-AVHRR_agg')
timevar = dataset.variables['time']
timeindex = date2index(date, timevar)  # find time index for desired date.
# read sst.  Will automatically create a masked array using
# missing_value variable attribute. 'squeeze out' singleton dimensions.
sst = dataset.variables['sst'][timeindex, :].squeeze()
# read ice.
ice = dataset.variables['ice'][timeindex, :].squeeze()
# read lats and lons (representing centers of grid boxes).
lats = dataset.variables['lat'][:]
lons = dataset.variables['lon'][:]
lons, lats = np.meshgrid(lons, lats)
# create figure, axes instances.
fig = plt.figure()
ax = fig.add_axes([0.05, 0.05, 0.9, 0.9])
# create Basemap instance.
# coastlines not used, so resolution set to None to skip
# continent processing (this speeds things up a bit)
m = Basemap(projection='kav7', lon_0=0, resolution=None)
# draw line around map projection limb.
# color background of map projection region.
# missing values over land will show up this color.
m.drawmapboundary(fill_color='0.3')
# plot sst, then ice with pcolor
im1 = m.pcolormesh(lons, lats, sst, shading='flat',
                   cmap=plt.cm.jet, latlon=True)
im2 = m.pcolormesh(lons, lats, ice, shading='flat',
                   cmap=plt.cm.gist_gray, latlon=True)
# draw parallels and meridians, but don't bother labelling them.
m.drawparallels(np.arange(-90., 99., 30.))
m.drawmeridians(np.arange(-180., 180., 60.))
# add colorbar
cb = m.colorbar(im1, "bottom", size="5%", pad="2%")
# add a title.
ax.set_title('SST and ICE analysis for %s' % date)
plt.show()


########################################################################
# Title: Play aruond with directories and path
# Date: 2015-08-20
########################################################################

import os
print __file__
print os.path.join(os.path.dirname(__file__), '..')
print os.path.dirname(os.path.realpath(__file__))
print os.path.abspath(os.path.dirname(__file__))


########################################################################
# Title: Test with array and tuples
# Date: 2015-09-20
########################################################################

import numpy as np
import scipy.stats as ss
import timeit

dimx = 100
dimy = 100
x = np.tile(np.arange(dimx) + 1, dimy)
y = np.array(sum([[val] * dimx for val in np.arange(dimy) + 1], []))
ind = zip(x, y)


def kernDist(index, center):
    distr = ss.multivariate_normal([0, 0], [[15, 0], [0, 15]])
    return distr.pdf([index[0] - center[0], index[1] - center[1]])


def eucDist(index, center):
    distr = ss.norm(0, 15**2)
    dist = np.sqrt((index[0] - center[0])**2 + (index[1] - center[1])**2)
    return distr.pdf(dist)

valKernel = np.array([kernDist(pairs, [dimx / 2, dimy / 2])
                      for pairs in ind]).reshape(dimx, dimy)
valEuc = np.array([eucDist(pairs, [dimx / 2, dimy / 2])
                   for pairs in ind]).reshape(dimx, dimy)

valKernel *= 1 / valKernel.sum()
valEuc *= 1 / valEuc.sum()

import matplotlib.pyplot as plt
plt.figure(1)
plt.subplot(211)
plt.imshow(valKernel)

plt.subplot(212)
plt.imshow(valEuc)
plt.show()

test = ['a', 'b', 'c']
'd' not in test

test = "But soft what light through yonder window breaks It is the east and Juliet is the sun Arise fair sun and kill the envious moon Who is already sick and pale with grief"


splitted = test.split(" ")
dic = ['Arise', "But"]
ind = [x in dic for x in splitted]
for ind, splitted in zip(ind, splitted):
    if ind == False:
        dic.append(splitted)

print dic.sort()

[dic.append(splitted) for ind, splitted in zip(ind, splitted) if ind == False]

fname = raw_input("Enter file name: ")
fh = open(fname)
lst = list()
for line in fh:
    texts = set(line.rstrip().split(" "))
    ind = [text in lst for text in texts]
    for ind, text in zip(ind, texts):
        if ind == False:
            lst.append(text)

lst.sort()
print lst

test = "From stephen.marquard@uct.ac.za Sat Jan  5 09:14:16 2008"

n = 0
for x in range(10):
    if x % 2 == 0:
        n += x
        print x


########################################################################
# Title: Programming test
# Date: 2015-09-30
########################################################################

def fizzBuzz(n):
    for i in [i + 1 for i in range(n)]:
        if ((i % 5 == 0) & (i % 3 == 0)):
            print "Fizz-Buzz"
        elif (i % 5 == 0):
            print "Buzz"
        elif (i % 3 == 0):
            print "Fizz"
        else:
            print i


########################################################################
# Title: Converting a number to a different base
# Date: 2015-12-23
########################################################################


def toBase(num, base):
    remainder = num
    newNum = ''
    while remainder > 0:
        newNum = str(remainder % 2) + newNum
        remainder /= base
    return int(newNum)

toBase(1000, 10)
toBase(1998, 2)


########################################################################
# Title: testing classes
# Date: 2016-02-17
########################################################################

class Parent(object):

    def __init__(self):
        self.name = "mom and dad"

    def altered(self):
        print "PARENT altered()"


class Child(Parent):

    def __init__(self):
        self.name = "son"

    def altered(self):
        print "CHILD, BEFORE PARENT altered()"
        super(Child, self).altered()
        super(Child, self).name
        print "CHILD, AFTER PARENT altered()"

dad = Parent()
son = Child()

dad.altered()
son.altered()

########################################################################
# Title: R or Python on Text Mining
# Date: 2016-11-02
# Source:
# https://datawarrior.wordpress.com/2015/08/12/codienerd-1-r-or-python-on-text-mining/
########################################################################

texts = ['I love Python.',
         'R is good for analytics.',
         'Mathematics is fun.']


# import all necessary libraries
from nltk.stem import PorterStemmer
from nltk.tokenize import SpaceTokenizer
from nltk.corpus import stopwords
from functools import partial
from gensim import corpora
from gensim.models import TfidfModel
import re

# initialize the instances for various NLP tools
tokenizer = SpaceTokenizer()
stemmer = PorterStemmer()

# define each steps
pipeline = [lambda s: re.sub('[^\w\s]', '', s),
            lambda s: re.sub('[\d]', '', s),
            lambda s: s.lower(),
            lambda s: ' '.join(filter(lambda s: not (
                s in stopwords.words()), tokenizer.tokenize(s))),
            lambda s: ' '.join(
                map(lambda t: stemmer.stem(t), tokenizer.tokenize(s)))
            ]

# function that carries out the pipeline step-by-step


def preprocess_text(text, pipeline):
    if len(pipeline) == 0:
        return text
    else:
        return preprocess_text(pipeline[0](text), pipeline[1:])

# preprocessing
preprocessed_texts = map(partial(preprocess_text, pipeline=pipeline), texts)

# converting to feature vectors
documents = map(lambda s: tokenizer.tokenize(s), texts)
corpus = [dictionary.doc2bow(document) for document in documents]
tfidfmodel = TfidfModel(corpus)
