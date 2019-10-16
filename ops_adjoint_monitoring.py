import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import os
import datetime
import cartopy.crs as ccrs

##################################      User input      ############################################

geos_ver = '525'
fp_or_fpp = 'fpp'
date = '20190826'
fields = ['u','v','tv','sphu']
level_to_plot = 'maxrms'       #Provide a level number or simply 'maxrms' to plot level of max spike

####################################################################################################


# Prepare file path and names
# ---------------------------
gmaj = geos_ver[0]
gmin = geos_ver[1:3]

y = date[0:4]
m = date[4:6]
d = date[6:8]

dtana = datetime.datetime(int(y), int(m), int(d), 0, 0, 0)
dtver = datetime.datetime(int(y), int(m), int(d)) + datetime.timedelta(hours=24)
dtr15 = datetime.datetime(int(y), int(m), int(d)) + datetime.timedelta(hours=-9)
dtr21 = datetime.datetime(int(y), int(m), int(d)) + datetime.timedelta(hours=-3)

dtanas = dtana.strftime("%Y%m%d_%H")
dtvers = dtver.strftime("%Y%m%d_%H")
dtr15s = dtr15.strftime("%Y%m%d_%H")
dtr21s = dtr21.strftime("%Y%m%d_%H")

y = dtr15s[0:4]
m = dtr15s[4:6]
d = dtr15s[6:8]

path = '/nfs3m/archive/sfa_cache01/projects/dao_ops/GEOS-'+gmaj+'.'+gmin+'/GEOSadas-'+gmaj+'_'+gmin\
       +'/f'+geos_ver+'_'+fp_or_fpp+'/prog/Y'+y+'/M'+m+'/D'+d+'/'

f15 = path+'H15/f'+geos_ver+'_'+fp_or_fpp+'.fsens_twe.eta.'+dtr15s+'z+'+dtvers+'z-'+dtanas+'z.nc4'
f21 = path+'H21/f'+geos_ver+'_'+fp_or_fpp+'.fsens_twe.eta.'+dtr21s+'z+'+dtvers+'z-'+dtanas+'z.nc4'

print("\n Looking for the following files:")
print(f15)
print(f21)

if (not os.path.exists(f15)):
    print("File for 15z does not exist, aborting")
    exit()

if (not os.path.exists(f21)):
    print("File for 21z does not exist, aborting")
    exit()

print("\n Demigrating archived files")
#os.system("dmget "+f15)
#os.system("dmget "+f21)


# Read the fields from the file
# -----------------------------

print("\n Reading fields from files")

# File handles
fh15 = Dataset(f15)
fh21 = Dataset(f21)

# Indices of arrays
nx = len(fh15.dimensions['lon'])
ny = len(fh15.dimensions['lat'])
nz = len(fh15.dimensions['lev'])
nf = len(fields)

# Arrays to hold field
f15 = np.zeros([nf, nz, ny, nx])
f21 = np.zeros([nf, nz, ny, nx])

# Arrays to hold field rms
f15_rms = np.zeros([nf, nz])
f21_rms = np.zeros([nf, nz])

# Dimension arrays
lons = np.zeros([nx])
lats = np.zeros([ny])
levs = np.zeros([nz])
lons = fh15.variables['lon'][:]
lats = fh15.variables['lat'][:]
levs = fh15.variables['lev'][:]

# Read fields from the files
for f in range(nf):
  print("  Reading "+fields[f])
  f15[f,:,:,:] = fh15.variables[fields[f]][:,:,:]
  f21[f,:,:,:] = fh21.variables[fields[f]][:,:,:]


# Compute RMS
# -----------

print("\n Computing RMS")

for f in range(nf):

  # Compute the RMS of each level
  for k in range(nz):

    f15_rms[f,k] = np.sqrt(np.mean(f15[f,k,:,:])**2)
    f21_rms[f,k] = np.sqrt(np.mean(f21[f,k,:,:])**2)


# Creat plots
# -----------

print("\n Creating plots")

for f in range(nf):

  # Make RMS plots
  # --------------
  pltfname = 'f'+geos_ver+'_'+fp_or_fpp+'.rms.'+fields[f]+'.'+dtr15s+'z+'+dtvers+'z-'+dtanas+'z.png'
  plt.figure(figsize=(6,12))
  plt.plot(f15_rms[f,:],levs,'b', label='15z sensitivity')
  plt.plot(f21_rms[f,:],levs,'r', label='21z sensitivity')
  plt.gca().invert_yaxis()
  plt.title("Root mean square of "+fields[f])
  plt.ylim(nz, 1)
  plt.legend(loc='upper right')
  plt.savefig(pltfname)


  # Make contour plots
  # ------------------
  if (level_to_plot == "maxrms"):
    pltlev15 = np.argmax(f15_rms[f,:])
    pltlev21 = np.argmax(f21_rms[f,:])
    print(pltlev15, pltlev21)
    if (f21_rms[f,pltlev21] >= f15_rms[f,pltlev15]):
      pltlev = pltlev21
    else:
      pltlev = pltlev15
  else:
    pltlev = int(level_to_plot)

  pltfname = 'f'+geos_ver+'_'+fp_or_fpp+'.contour.'+fields[f]+'.lev-'+str(pltlev)+'.'+dtr15s+'z+'+\
             dtvers+'z-'+dtanas+'z.png'

  f15plt = f15[f,pltlev,:,:]
  f21plt = f21[f,pltlev,:,:]

  clev_max = np.max([np.abs(f15plt), np.abs(f21plt)])
  clev_inc = 2*clev_max/12.
  clev = np.arange(-clev_max,clev_max,clev_inc)

  fig = plt.figure(figsize=(12, 12))
  ax = fig.add_subplot(2, 1, 1, projection=ccrs.PlateCarree())
  ax.contourf(lons, lats, f15plt, clev,
              transform=ccrs.PlateCarree(), cmap='seismic')
  ax.coastlines()
  ax.set_global()

  ax = fig.add_subplot(2, 1, 2, projection=ccrs.PlateCarree())
  ax.contourf(lons, lats, f21plt, clev,
              transform=ccrs.PlateCarree(), cmap='seismic')
  ax.coastlines()
  ax.set_global()

  fig.savefig(pltfname)