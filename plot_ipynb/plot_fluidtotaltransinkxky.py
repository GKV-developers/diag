#!/usr/bin/env python
# coding: utf-8

# # Plot fluidtotaltransinkxky*

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
dtout_eng=nml["times"]["dtout_eng"]

print("nx =",nx)
print("global_ny =",global_ny)
print("dtout_eng =",dtout_eng)


# In[ ]:


### Load data of fluidtotaltransinkxky* ###
filelist=sorted(glob.glob("./data/fluidtotaltransinkxky_s0_t*dat"))

# Time
t=[]
for f in filelist:
    wt=read_time(f)
    t.append(wt)
t=np.array(t)
print("t =",t)

# Coordinates: kx, ky
data=np.loadtxt(filelist[0]) # Ascii data
data=data.reshape(global_ny+1,2*nx+1,4)
kx=data[0,:,0]
ky=data[:,0,1]
print("kx =",kx)
print("ky =",ky)    

# Values from fluidtotaltransinkxky*
tk_es=[]
tk_em=[]
for f in filelist:
    data=np.loadtxt(f) # Ascii data
    data=data.reshape(global_ny+1,2*nx+1,4)
#     kx=data[0,:,0]
#     ky=data[:,0,1]
    tk_es.append(data[:,:,2])
    tk_em.append(data[:,:,3])

tk_es=np.array(tk_es)
tk_em=np.array(tk_em)
print(tk_es.shape)
print(tk_em.shape)  


# In[ ]:


it=10
fig=plt.figure()
ax=fig.add_subplot(111)
vmax=np.abs(tk_es[it,:,:]).max()
quad=ax.pcolormesh(kx,ky,tk_es[it,:,:],shading="auto",cmap="RdBu_r",vmin=-vmax,vmax=vmax)
cbar=fig.colorbar(quad,shrink=1.0,aspect=5)
ax.set_title(r"Electrostatic nonlinear transfer $I_E$ (t={})".format(t[it]))
ax.set_xlabel(r"Radial wavenumber $k_x$")
ax.set_ylabel(r"Poloidal wavenumber $k_y$")
plt.show()


# In[ ]:


mx=nx
my=2
fig=plt.figure()
ax=fig.add_subplot(111)
ax.plot(t,tk_es[:,my,mx],label=r"$k_x$={},$k_y$={}".format(kx[mx],ky[my]))
ax.set_xlabel("Time t")
ax.set_ylabel(r"Electrostatic nonlinear transfer $I_E$")
plt.show()


# In[ ]:




