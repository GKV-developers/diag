#!/usr/bin/env python
# coding: utf-8

# # Plot phiinrz*

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
print("dtout_ptn =",dtout_ptn)


# In[ ]:


### Load data of mominxy* ###
filelist=sorted(glob.glob("./data/phiinrz_t*dat"))

# Time
t=[]
for f in filelist:
    wt=read_time(f)
    t.append(wt)
t=np.array(t)
print("t =",t)

# Coordinates: xx, yy
data=np.loadtxt(filelist[0]) # Ascii data
print(data.shape)
nzw=int((data.shape[0]/(2*nxw+1)-1)/2)
print(nzw)
data=data.reshape(2*nzw+1,2*nxw+1,3)
mr=data[:,:,0]
z_car=data[:,:,1]
print("mr =",mr)
print("z_car =",z_car)    

# Values from mominxy*
phi=[]
for f in filelist:
    data=np.loadtxt(f) # Ascii data
    data=data.reshape(2*nzw+1,2*nxw+1,3)
#     mr=data[:,:,0]
#     z_car=data[:,:,1]
#     phi=data[:,:,2]
#     upara=data[:,:,3]
#     ppara=data[:,:,4]
#     pperp=data[:,:,5]
#     qlpara=data[:,:,6]
#     qlperp=data[:,:,7]
    phi.append(data[:,:,2])
    
phi=np.array(phi)
print(phi.shape)  


# In[ ]:


it=40
fig=plt.figure(figsize=(6,6))
ax=fig.add_subplot(111)
quad=ax.pcolormesh(mr,z_car,phi[it,:,:],cmap="jet",shading="auto")
vmax=np.max(np.abs(phi[it,:,:]))
quad.set_clim(-vmax,vmax)
cbar=fig.colorbar(quad,shrink=1.0,aspect=5)
ax.set_title(r"Electrostatic potential $\phi$ (t={})".format(t[it]))
ax.set_aspect('equal')
ax.set_xlabel(r"Major radius $R/R_0$")
ax.set_ylabel(r"Height $Z/R_0$")
plt.show()


# In[ ]:


### Example of animation ###
from matplotlib.animation import FuncAnimation
from IPython.display import HTML

fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
ax.set_aspect('equal')
ax.set_xlabel(r"Major radius $R/R_0$")
ax.set_ylabel(r"Height $Z/R_0$")
title=ax.set_title(r"Electrostatic potential $\phi$ (t={})".format(t[0]))
quad=ax.pcolormesh(mr,z_car,phi[0,:,:],cmap="jet",shading="auto")
vmax=np.max(np.abs(phi[0,:,:]))
quad.set_clim(-vmax,vmax)
cbar=fig.colorbar(quad,shrink=1.0,aspect=5)

def update_quad(i):
    title.set_text(r"Electrostatic potential $\phi$ (t={})".format(t[i]))
    quad.set_array(phi[i,:,:].flatten())
    vmax=np.max(np.abs(phi[i,:,:]))
    quad.set_clim(-vmax,vmax)

ani = FuncAnimation(fig, update_quad,
                    frames=range(0,len(t),10), interval=100)
#ani.save('advection.mp4', writer="ffmpeg", dpi=100)
#plt.close()

HTML(ani.to_jshtml())


# In[ ]:





# In[ ]:




