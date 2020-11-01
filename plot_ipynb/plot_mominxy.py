#!/usr/bin/env python
# coding: utf-8

# # Plot mominxy*

# In[ ]:


import numpy as np
import matplotlib.pyplot as plt
import glob
import f90nml
from read_f90 import read_parameters, read_time # In-house module for GKV-diag outputs

### GKV parameters from gkvp_header.f90 ###
nxw=read_parameters("../src/gkvp_header.f90", "nxw", int)
nyw=read_parameters("../src/gkvp_header.f90", "nyw", int)

### GKV parameters from gkvp_namelist ###
nml=f90nml.read("../gkvp_namelist.001")
#print(nml)
dtout_ptn=nml["times"]["dtout_ptn"]

print("nxw =",nxw)
print("nyw =",nyw)
print("dtout_ptn =",dtout_ptn)


# In[ ]:


### Load data of mominxy* ###
filelist=sorted(glob.glob("./data/mominxy_z0000s0_t*dat"))

# Time
t=[]
for f in filelist:
    wt=read_time(f)
    t.append(wt)
t=np.array(t)
print("t =",t)

# Coordinates: xx, yy
data=np.loadtxt(filelist[0]) # Ascii data
data=data.reshape(2*nyw,2*nxw,8)
xx=data[0,:,0]
yy=data[:,0,1]
print("xx =",xx)
print("yy =",yy)    

# Values from mominxy*
dens=[]
for f in filelist:
    data=np.loadtxt(f) # Ascii data
    data=data.reshape(2*nyw,2*nxw,8)
#     x=data[0,:,0]
#     y=data[:,0,1]
#     dens=data[:,:,2]
#     upara=data[:,:,3]
#     ppara=data[:,:,4]
#     pperp=data[:,:,5]
#     qlpara=data[:,:,6]
#     qlperp=data[:,:,7]
    dens.append(data[:,:,2])
    
dens=np.array(dens)
print(dens.shape)  


# In[ ]:


it=10
fig=plt.figure()
ax=fig.add_subplot(111)
quad=ax.pcolormesh(xx,yy,dens[it,:-1,:-1],cmap="jet")
vmax=np.max(np.abs(dens[it,:-1,:-1]))
quad.set_clim(-vmax,vmax)
cbar=fig.colorbar(quad,shrink=1.0,aspect=5)
ax.set_title(r"$Density n$ (t={})".format(t[it]))
ax.set_xlabel(r"Radial direction $x$")
ax.set_ylabel(r"Poloidal direction $y$")
plt.show()


# In[ ]:


### Example of animation ###
from matplotlib.animation import FuncAnimation
from IPython.display import HTML

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlabel(r"Radial direction $x$")
ax.set_ylabel(r"Poloidal direction $y$")
title=ax.set_title(r"$Density n$ (t={})".format(t[0]))
quad=ax.pcolormesh(xx,yy,dens[0,:-1,:-1],cmap="jet")
vmax=np.max(np.abs(dens[0,:-1,:-1]))
quad.set_clim(-vmax,vmax)
cbar=fig.colorbar(quad,shrink=1.0,aspect=5)

def update_quad(i):
    title.set_text(r"$Density n$ (t={})".format(t[i]))
    quad.set_array(dens[i,:-1,:-1].flatten())
    vmax=np.max(np.abs(dens[i,:-1,:-1]))
    quad.set_clim(-vmax,vmax)

ani = FuncAnimation(fig, update_quad,
                    frames=range(0,len(t),10), interval=100)
#ani.save('advection.mp4', writer="ffmpeg", dpi=100)
#plt.close()

HTML(ani.to_jshtml())


# In[ ]:





# In[ ]:




