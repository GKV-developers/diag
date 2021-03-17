#!/usr/bin/env python
# coding: utf-8

# # Plot fluiddetailtransinkxky*

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
filelist=sorted(glob.glob("./data/fluiddetailtransinkxky_x0000y0002s0_t*dat"))

# Time
t=[]
for f in filelist:
    wt=read_time(f)
    t.append(wt)
t=np.array(t)
print("t =",t)

# Coordinates: kx, ky
data=np.loadtxt(filelist[0]) # Ascii data
data=data.reshape(2*global_ny+1,2*nx+1,8)
kx=data[0,:,0]
ky=data[:,0,1]
print("kx =",kx)
print("ky =",ky)    

# Values from fluidtotaltransinkxky*
jkpq_es=[]
jpqk_es=[]
jqkp_es=[]
jkpq_em=[]
jpqk_em=[]
jqkp_em=[]
for f in filelist:
    data=np.loadtxt(f) # Ascii data
    data=data.reshape(2*global_ny+1,2*nx+1,8)
#     kx=data[0,:,0]
#     ky=data[:,0,1]
    jkpq_es.append(data[:,:,2])
    jpqk_es.append(data[:,:,3])
    jqkp_es.append(data[:,:,4])
    jkpq_em.append(data[:,:,5])
    jpqk_em.append(data[:,:,6])
    jqkp_em.append(data[:,:,7])
    
jkpq_es=np.array(jkpq_es)
jpqk_es=np.array(jpqk_es)
jqkp_es=np.array(jqkp_es)
jkpq_em=np.array(jkpq_em)
jpqk_em=np.array(jpqk_em)
jqkp_em=np.array(jqkp_em)
print(jkpq_es.shape)


# In[ ]:


it=3
fig=plt.figure()
ax=fig.add_subplot(111)
vmax=np.abs(jkpq_es[it,:,:]).max()
quad=ax.pcolormesh(kx,ky,jkpq_es[it,:,:],shading="auto",cmap="RdBu_r",vmin=-vmax,vmax=vmax)
cbar=fig.colorbar(quad,shrink=1.0,aspect=5)
ax.set_title(r"Electrostatic triad transfer $J_{E,k}^{p,q}(p_x,p_y)$ "+"(t={})".format(t[it]))
ax.set_xlabel(r"Radial wavenumber $k_x$")
ax.set_ylabel(r"Poloidal wavenumber $k_y$")
plt.show()
fig=plt.figure()
ax=fig.add_subplot(111)
# vmax=np.abs(jpqk_es[it,:,:]).max()
quad=ax.pcolormesh(kx,ky,jpqk_es[it,:,:],shading="auto",cmap="RdBu_r",vmin=-vmax,vmax=vmax)
cbar=fig.colorbar(quad,shrink=1.0,aspect=5)
ax.set_title(r"Electrostatic triad transfer $J_{E,p}^{q,k}(p_x,p_y)$ "+"(t={})".format(t[it]))
ax.set_xlabel(r"Radial wavenumber $k_x$")
ax.set_ylabel(r"Poloidal wavenumber $k_y$")
plt.show()
# fig=plt.figure()
# ax=fig.add_subplot(111)
# # vmax=np.abs(jqkp_es[it,:,:]).max()
# quad=ax.pcolormesh(kx,ky,jqkp_es[it,:,:],shading="auto",cmap="RdBu_r",vmin=-vmax,vmax=vmax)
# cbar=fig.colorbar(quad,shrink=1.0,aspect=5)
# ax.set_title(r"Electrostatic triad transfer $J_{E,q}^{k,p}(p_x,p_y)$ "+"(t={})".format(t[it]))
# ax.set_xlabel(r"Radial wavenumber $k_x$")
# ax.set_ylabel(r"Poloidal wavenumber $k_y$")
# plt.show()
# fig=plt.figure()
# ax=fig.add_subplot(111)
# # vmax=np.abs(jkpq_es[it,:,:]+jpqk_es[it,:,:]+jqkp_es[it,:,:]).max()
# quad=ax.pcolormesh(kx,ky,jkpq_es[it,:,:]+jpqk_es[it,:,:]+jqkp_es[it,:,:],shading="auto",cmap="RdBu_r",vmin=-vmax,vmax=vmax)
# cbar=fig.colorbar(quad,shrink=1.0,aspect=5)
# ax.set_title(r"Check detailed balance"+"(t={})".format(t[it]))
# ax.set_xlabel(r"Radial wavenumber $k_x$")
# ax.set_ylabel(r"Poloidal wavenumber $k_y$")
# plt.show()


# In[ ]:


fig=plt.figure()
ax=fig.add_subplot(111)
ax.plot(t,np.sum(jkpq_es[:,:,:],axis=(1,2)))
ax.set_xlabel("Time t")
ax.set_ylabel(r"Electrostatic nonlinear transfer $I_E$")
plt.show()


# In[ ]:





# In[ ]:




