#!/usr/bin/env python
# coding: utf-8

# # Plot mominkxky*

# In[ ]:


import numpy as np
import matplotlib.pyplot as plt
import glob
import f90nml
from read_f90 import read_parameters, read_time # In-house module for GKV-diag outputs

### GKV parameters from gkvp_header.f90 ###
nx=read_parameters("../src/gkvp_header.f90", "nx", int)
global_ny=read_parameters("../src/gkvp_header.f90", "global_ny", int)

### GKV parameters from gkvp_namelist ###
nml=f90nml.read("../gkvp_namelist.001")
#print(nml)
dtout_ptn=nml["times"]["dtout_ptn"]

print("nx =",nx)
print("global_ny =",global_ny)
print("dtout_ptn =",dtout_ptn)


# In[ ]:


### Load data of mominkxky* ###
filelist=sorted(glob.glob("./data/mominkxky_s0_t*dat"))

# Time
t=[]
for f in filelist:
    wt=read_time(f)
    t.append(wt)
t=np.array(t)
print("t =",t)

# Coordinates: kx, ky
data=np.loadtxt(filelist[0]) # Ascii data
data=data.reshape(global_ny+1,2*nx+1,8)
kx=data[0,:,0]
ky=data[:,0,1]
print("kx =",kx)
print("ky =",ky)    

# Values from mominkxky*
dens=[]
for f in filelist:
    data=np.loadtxt(f) # Ascii data
    data=data.reshape(global_ny+1,2*nx+1,8)
#     kx=data[0,:,0]
#     ky=data[:,0,1]
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
quad=ax.pcolormesh(kx,ky,dens[it,:-1,:-1])
cbar=fig.colorbar(quad,shrink=1.0,aspect=5)
ax.set_title(r"$<|n_k|^2>$ (t={})".format(t[it]))
ax.set_xlabel(r"Radial wavenumber $k_x$")
ax.set_ylabel(r"Poloidal wavenumber $k_y$")
plt.show()


# In[ ]:


### Example of animation ###
from matplotlib.animation import FuncAnimation
from IPython.display import HTML

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlabel(r"Radial wavenumber $k_x$")
ax.set_ylabel(r"Poloidal wavenumber $k_y$")
title=ax.set_title(r"$<|n_k|^2>$ (t={})".format(t[0]))
quad=ax.pcolormesh(kx,ky,dens[0,:-1,:-1])
vmax=np.max(dens[0,:-1,:-1])
quad.set_clim(0,vmax)
cbar=fig.colorbar(quad,shrink=1.0,aspect=5)

def update_quad(i):
    title.set_text(r"$<|n_k|^2>$ (t={})".format(t[i]))
    quad.set_array(dens[i,:-1,:-1].flatten())
    vmax=np.max(dens[i,:-1,:-1])
    quad.set_clim(0,vmax)

ani = FuncAnimation(fig, update_quad,
                    frames=range(0,len(t),10), interval=100)
#ani.save('advection.mp4', writer="ffmpeg", dpi=100)
#plt.close()

HTML(ani.to_jshtml())


# In[ ]:





# In[ ]:





# In[ ]:




