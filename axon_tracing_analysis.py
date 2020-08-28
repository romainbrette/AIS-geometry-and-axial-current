
"""
A script to measure AIS geometry from 3D axon tracing.
"""

from brian2 import *
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.image as mpimg
from skimage import io
from mpl_toolkits.mplot3d import axes3d
from scipy import interpolate
import pandas as pd


### Path to files 
path= 'data/RGC data/images/traced axons/' 

### Loading cell
day = 20200228
retina = 'A'
cell = 1
images_carac = pd.read_excel('RGC_image_caracteristics.xlsx', dtype='object')
line = images_carac.loc[(images_carac.Day == day)& (images_carac.Retina == retina ) & (images_carac.Cell == cell )]
file_name = np.array(line['File name'])[0]
print (file_name)

### Load data points coordinates of axonal profile drawn in Vaa3d and soma coordinates
# the last point corresponds to the soma
from io import StringIO
data = open(path + '%s%s%saxon3d.swc' %(str(day), retina, cell)).read().replace(',','.')
data = np.loadtxt(StringIO(data))

z_spacing = 0.5 #um

### In some cells, too much data points were drawn within the soma, so we remove them.

if (day, retina, cell) == (20190418, 'B', 1): # soma is the last point
    sorted_idx = argsort(data[:,3]) # y coords
    y = np.flip(data[:,3][sorted_idx])[1:]
    x = np.flip(data[:,2][sorted_idx])[1:]
    z = np.flip(data[:,4][sorted_idx])[1:]
    radius = np.flip(data[:,5][sorted_idx])[1:]
    
elif (day, retina, cell) == (20191115, 'B', 2): # the profile is made of 2 segments
    x = hstack((np.flip(data[:25,2]), np.flip(data[25:,2]))) # pixel coord
    y = hstack((np.flip(data[:25,3]), np.flip(data[25:,3])))
    z = hstack((np.flip(data[:25,4]), np.flip(data[25:,4])))
    radius = hstack((np.flip(data[:25,5]), np.flip(data[25:,5])))
    
elif (day, retina, cell) == (20191115, 'A', 1): # the axon profile turns out of the staining at the beginning of the axon
    sorted_idx = argsort(data[:,3]) # y coords
    y = np.delete(np.flip(data[:,3]), [ 1,2,3,4])
    x = np.delete(np.flip(data[:,2]), [ 1,2,3,4])
    z = np.delete(np.flip(data[:,4]), [ 1,2,3,4])
    radius = np.delete(np.flip(data[:,5]), [ 1,2,3,4])
    
elif (day, retina, cell) == (20191114, 'D', 1): # removing the first point that is within the soma
    # axon3d_bis
    x = hstack((np.flip(data[:35,2]), np.flip(data[35:,2]))) # pixel coord
    y = hstack((np.flip(data[:35,3]), np.flip(data[35:,3])))
    z = hstack((np.flip(data[:35,4]), np.flip(data[35:,4])))
    radius = hstack((np.flip(data[:35,5]), np.flip(data[35:,5])))
    
elif (day, retina, cell) == (20191031, 'R', 1): # the profile is made of 2 segments
    # axon3d_bis
    x = hstack((np.flip(data[:43,2]), np.flip(data[43:,2]))) # pixel coord
    y = hstack((np.flip(data[:43,3]), np.flip(data[43:,3])))
    z = hstack((np.flip(data[:43,4]), np.flip(data[43:,4])))
    radius = hstack((np.flip(data[:43,5]), np.flip(data[43:,5])))
    
elif (day, retina, cell) == (20191126, 'A', 1): # the profile is made of 2 segments
    # axon3d_bis
    x = hstack((np.flip(data[:36,2]), np.flip(data[36:,2]))) # pixel coord
    y = hstack((np.flip(data[:36,3]), np.flip(data[36:,3])))
    z = hstack((np.flip(data[:36,4]), np.flip(data[36:,4])))
    radius = hstack((np.flip(data[:36,5]), np.flip(data[36:,5])))
    
elif (day, retina, cell) == (20191126, 'B', 1): # the profile is made of 4 segments
    # axon3d_bis
    x = hstack((np.flip(data[1:29,2]), np.flip(data[40:54,2]), np.flip(data[29:40,2]), np.flip(data[54:68,2]))) # pixel coord
    y = hstack((np.flip(data[1:29,3]), np.flip(data[40:54,3]), np.flip(data[29:40,3]), np.flip(data[54:68,3])))
    z = hstack((np.flip(data[1:29,4]), np.flip(data[40:54,4]), np.flip(data[29:40,4]), np.flip(data[54:68,4])))
    radius = hstack((np.flip(data[1:29,5]), np.flip(data[40:54,5]), np.flip(data[29:40,5]), np.flip(data[54:68,5])))
    
elif (day, retina, cell) == (20190419, 'A', 1): # removing a few points after the start pt because the diameter is over-estimated (dendrite just next)
    sorted_idx = argsort(data[:,3]) # y coords
    y = np.delete(np.flip(data[:,3]), [1,2,3,5,6])#[1:]
    x = np.delete(np.flip(data[:,2]), [1,2,3,5,6])#[1:]
    z = np.delete(np.flip(data[:,4]), [1,2,3,5,6])#[1:]
    radius = np.delete(np.flip(data[:,5]), [1,2,3,5,6])#[1:]
    radius[0] = 5.
    radius[1] = radius[7]/2
    
elif (day, retina, cell) == (20190419, 'B', 1): # soma is the last point
    sorted_idx = argsort(data[:,3]) # y coords
    y = np.flip(data[:,3])[1:]
    x = np.flip(data[:,2])[1:]
    z = np.flip(data[:,4])[1:]
    radius = np.flip(data[:,5])[1:]
    
elif (day, retina, cell) == (20200213, 'B', 1): # the axon profile turns out of the staining at the beginning of the axon
    sorted_idx = argsort(data[:,3]) # y coords
    y = np.delete(np.flip(data[:,3]), [ 1,2,3])
    x = np.delete(np.flip(data[:,2]), [ 1,2,3])
    z = np.delete(np.flip(data[:,4]), [ 1,2,3])
    radius = np.delete(np.flip(data[:,5]), [ 1,2,3])
    
elif (day, retina, cell) == (20200214, 'B', 1): # the axon profile turns out of the staining at the beginning of the axon
    sorted_idx = argsort(data[:,3]) # y coords
    y = np.delete(np.flip(data[:,3]), [ 1,2,3])
    x = np.delete(np.flip(data[:,2]), [ 1,2,3])
    z = np.delete(np.flip(data[:,4]), [ 1,2,3])
    radius = np.delete(np.flip(data[:,5]), [ 1,2,3])
    
elif (day, retina, cell) == (20200220, 'A', 1): # the axon profile turns out of the staining at the beginning of the axon
    sorted_idx = argsort(data[:,3]) # y coords
    y = np.delete(np.flip(data[:,3]), [0,1])
    x = np.delete(np.flip(data[:,2]), [0,1])
    z = np.delete(np.flip(data[:,4]), [0,1])
    radius = np.delete(np.flip(data[:,5]), [0,1])
    
    z_spacing = 0.46 #um
    
elif (day, retina, cell) == (20200222, 'B', 1): # the axon profile turns out of the staining at the beginning of the axon
    sorted_idx = argsort(data[:,3]) # y coords
    y = np.delete(np.flip(data[:,3]), [1, 2, 3])
    x = np.delete(np.flip(data[:,2]), [1, 2, 3])
    z = np.delete(np.flip(data[:,4]), [1, 2, 3])
    radius = np.delete(np.flip(data[:,5]), [1, 2, 3])

elif (day, retina, cell) == (20200228, 'A', 1): # the axon profile turns out of the staining at the beginning of the axon
    sorted_idx = argsort(data[:,3]) # y coords
    y = np.delete(np.flip(data[:,3]), [0, 1])
    x = np.delete(np.flip(data[:,2]), [0, 1,])
    z = np.delete(np.flip(data[:,4]), [0, 1])
    radius = np.delete(np.flip(data[:,5]), [0, 1])
   
else:
    x = np.flip(data[:,2]) # pixel coord
    y = np.flip(data[:,3]) # pixel coord
    z = np.flip(data[:,4]) # depth coord
    radius = np.flip(data[:, 5]) 
    

### Load full 3D stack (indices: z,y,x)
ankG = io.imread('data/RGC data//images/%s/Retina %s/cell %s/%s%s%sankG.tif' %(str(file_name), retina, cell, str(day), retina, cell))

### Stack caracteristics 
tot_size = np.array(line['Image size (um)'])[0] # image size in um
tot_pix_num = ankG.shape[1] # number of pixels in a row/column
pix_size = tot_size/tot_pix_num # pixel size in um
z_n = ankG.shape[0] # number of stacks
z_depth = z_n * z_spacing

### Load and plot biocytin max projection
plt.figure('Biocytin max Z projection', figsize=(20,10))
img = mpimg.imread(path + 'MAX_%s%s%sComposite.png' %(str(day), retina, cell))
img_bio = mpimg.imread(path + 'MAX_%s%s%sbiocytin.png' %(str(day), retina, cell))

ax1 = plt.subplot(121)
ax1.imshow(img, cmap='gray')
ax1.plot(x, y,'.')
ax1.plot(x[0], y[0], 'ko', label='soma')
ax1.set_ylabel('y (pixel)')
ax1.set_xlabel('x (pixel)')
ax1.legend()

ax2 = plt.subplot(122)
ax2.plot(y, radius, 'o')
ax2.set_ylabel('Diameter (um)')

### We remove double points
double_idx = []
for i in range(len(x)):
    for j in range(i+1, len(x)):
        if x[i] == [j] and y[i] == y[j] and z[i] == z[j]:
            print ('deleting', j)
            double_idx.append(j)
            
x_corr = np.delete(x, double_idx)# a 1D array that contains the index in direction x of each pixel along the axon (index 0 corresponds to the soma)
y_corr = np.delete(y, double_idx)# a 1D array that contains the index in direction y of each pixel along the axon (index 0 corresponds to the soma)
z_corr = np.delete(z, double_idx)# a 1D array that contains the index in direction z of each pixel along the axon (index 0 corresponds to the soma)
radius_corr = np.delete(radius, double_idx)

### Interpolate data points to obtain a finer spacing (pixel size)
### 3D curve interpolation

xy_len_tot = 0 # total length of the axon profile in xy plane 
xyz_len_tot = 0 # total length of the axon profile in 3D 
for i in range(1, len(x_corr) ):
    dx = x_corr[i] - x_corr[i-1]
    dy = y_corr[i] - y_corr[i-1]
    dz = z_corr[i] - z_corr[i-1]
    di = np.sqrt(dx**2+dy**2) 
    di_z = np.sqrt((dx)**2+(dy)**2+(dz)**2) 
    xy_len_tot += di * pix_size
    xyz_len_tot += di_z 

### Interpolation and plotting
### !!! splprep returns a Value Error when to points are identical, so we should remove identical points
tck, u = interpolate.splprep([x_corr,y_corr,z_corr, radius_corr], s=3) # s controls the degree of smootheness. 
u_fine = np.linspace(0,1, int(xyz_len_tot)) # to interpolate on points spaced by one pixel in the xy plane
xf, yf, zf, radius_f = interpolate.splev(u_fine, tck)

if (day, retina, cell) == (20191031, 'R', 1): # the profile starts too far in the soma
    xf = xf[3:]
    yf = yf[3:]
    zf = zf[3:]
    radius_f = radius_f[3:]

ax1.plot(xf, yf, 'w.')
ax2.plot(yf, radius_f, '.')

### Get fluorescence intensity along the axonal profile
### we need to convert the point coordinates to indices of the image
n = len(xf)
x_pix = np.array([round(xf[i]) for i in range(n)]) 
y_pix = np.array([round(yf[i]) for i in range(n)]) 
z_pix = np.array([round(zf[i]) for i in range(n)]) 

### Length of axonal profile in um
axon_len_um = 0 # total length of the axon profile in um
axon_cum_len_um = np.zeros(n) # filling an array with the cumulative length at each point of the axonal profile
distances_um = np.zeros(n)
for i in range(1, n):
    dx = xf[i] - xf[i-1]
    dy = yf[i] - yf[i-1]
    dz = zf[i] - zf[i-1]
    di_z = np.sqrt((pix_size*dx)**2+(pix_size*dy)**2+(z_spacing*dz)**2) 
    distances_um[i] = di_z
    axon_len_um += di_z
    axon_cum_len_um[i] = axon_len_um
    
print ('Profile length:', axon_len_um, 'um')

### Intensity along axonal profile 
fi = [ankG[int(z_pix[i]), int(y_pix[i]), int(x_pix[i])] for i in range(n)] # !!! y then x

### Smoothing the intensity profile with a sliding mean
fi_slide = np.zeros(n)
d = 15 # half-window, i.e. number of pixels on each side

for i in range(n):
    if i < d: # start of the axon, not full window
        fi_slide[i] = np.mean(fi[0:i+d])
    elif i > n-d: # end of the axon, not full window
        fi_slide[i] = np.mean(fi[i-d:n])
    else: 
        fi_slide[i] = np.mean(fi[i-d:i+d])
        
fi_norm = [(fi[i]-min(fi))/(max(fi)-min(fi)) for i in range(n)]
fi_slide_norm = [(fi_slide[i]-min(fi_slide))/(max(fi_slide)-min(fi_slide)) for i in range(n)]

### By eye measurement of AIS start and end positions

coords = []
closed = False

def handle_close(evt):
    global closed
    closed = True

def waitforbuttonpress():
    while plt.waitforbuttonpress(0.2) is None:
        if closed:
            return False
    return True

def onclick(event):
    global ix, iy
    ix, iy = event.xdata, event.ydata
    print ('x = %0.2f, y = %0.2f'%(ix, iy))

    global coords
    coords.append((ix, iy))

    if len(coords) == 2:
        fig_ip.canvas.mpl_disconnect(cid)
        
    return coords

fig_ip = plt.figure('Raw intensity profile', figsize=(15,10))

cid = fig_ip.canvas.mpl_connect('button_press_event', onclick)
fig_ip.canvas.mpl_connect('close_event', handle_close)

while True:
    plt.subplot(211)
    plt.plot(axon_cum_len_um, fi, 'k', label='raw')
    plt.plot(axon_cum_len_um, fi_slide, 'k--', label='raw smoothed')
    plt.ylabel('Intensity')
    plt.legend(frameon=False)
    
    ### Normalized IP
    fi_slide_norm = [(fi_slide[i]-min(fi_slide))/(max(fi_slide)-min(fi_slide)) for i in range(n)]
    
    plt.subplot(212)
    plt.plot(axon_cum_len_um, fi_slide_norm, 'k')
    plt.xlabel('Position along axonal profile (um)')
    plt.ylabel('Normalized smoothed intensity')

    if not waitforbuttonpress():
        break
    print('.')


### AIS start
# idx along axonal profile
idx_start = np.argmin(abs(axon_cum_len_um - coords[0][0]))
x_ais_start = xf[idx_start]
y_ais_start = yf[idx_start]
idx_end = np.argmin(abs(axon_cum_len_um - coords[1][0]))
x_ais_end = xf[idx_end]
y_ais_end = yf[idx_end]

ax1.plot(x_ais_start, y_ais_start, 'ko')
ax1.plot(x_ais_end, y_ais_end, 'ko')

print( 'AIS start:', coords[0][0], 'um')
print ('AIS end:', coords[1][0], 'um')

### Axial resistance estimation
Ri = 100. * ohm * cm

sum_diameters = 0
for i in range(idx_start):
    diam = 2 * radius_f[i] * pix_size * um
    #print (i, diam)
    sum_diameters += distances_um[i] * um/diam**2

# Axial resistance
Ra = ((4*Ri)/pi)  * sum_diameters 
print ('Diam AIS start:', 2 * radius_f[idx_start] * pix_size * um)
print ('Diam AIS end:', 2 * radius_f[idx_end] * pix_size * um)

print ('Axial resistance:', Ra)

plt.show()










