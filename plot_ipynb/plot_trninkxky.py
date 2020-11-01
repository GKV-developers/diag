#!/usr/bin/env python
# coding: utf-8

# # Plot trninkxky*

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


### Load data of trninkxky* ###
filelist=sorted(glob.glob("./data/trninkxky_s0_t*dat"))

# Time
t=[]
for f in filelist:
    wt=read_time(f)
    t.append(wt)
t=np.array(t)
print("t =",t)

# Coordinates: kx, ky
data=np.loadtxt(filelist[0]) # Ascii data
data=data.reshape(global_ny+1,2*nx+1,14)
kx=data[0,:,0]
ky=data[:,0,1]
print("kx =",kx)
print("ky =",ky)    

# Values from trninkxky*
neint=[]
nmint=[]
for f in filelist:
    data=np.loadtxt(f) # Ascii data
    data=data.reshape(global_ny+1,2*nx+1,14)
#     kx=data[0,:,0]
#     ky=data[:,0,1]
#     entrpy=data[:,:,2]
#     fenegy=data[:,:,3]
#     menegy=data[:,:,4]
#     peint=data[:,:,5]
#     pmint=data[:,:,6]
#     neint=data[:,:,7]
#     nmint=data[:,:,8]
#     dcd=data[:,:,9]
#     pflux_es=data[:,:,10]
#     pflux_em=data[:,:,11]
#     eflux_es=data[:,:,12]
#     eflux_em=data[:,:,13]
    neint.append(data[:,:,7])
    nmint.append(data[:,:,8])

neint=np.array(neint)
nmint=np.array(nmint)
print(neint.shape)
print(nmint.shape)  


# In[ ]:


it=10
fig=plt.figure()
ax=fig.add_subplot(111)
quad=ax.pcolormesh(kx,ky,neint[it,:-1,:-1])
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
ax.plot(t,neint[:,my,mx],label=r"$k_x$={},$k_y$={}".format(kx[mx],ky[my]))
ax.set_xlabel("Time t")
ax.set_ylabel(r"Electrostatic nonlinear transfer $I_E$")
plt.show()


# # Plot triinkxky*

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
#print(gkvnml)
dtout_ptn=nml["times"]["dtout_ptn"]

print("nx =",nx)
print("global_ny =",global_ny)
print("dtout_ptn =",dtout_ptn)


# In[ ]:


### Load data of triinkxky* ###
filelist=sorted(glob.glob("./data/triinkxky_s0mx0000my0002_t*dat"))

# Time
t_tri=[]
for f in filelist:
    wt=read_time(f)
    t_tri.append(wt)
t_tri=np.array(t_tri)
print("t_tri =",t_tri)

# Coordinates: kx, ky
data=np.loadtxt(filelist[0]) # Ascii data
data=data.reshape(2*global_ny+1,2*nx+1,8)
kx_tri=data[0,:,0]
ky_tri=data[:,0,1]
print("kx_tri =",kx_tri)
print("ky_tri =",ky_tri)    

# Values from triinkxky*
jkpq_es=[]
jkpq_em=[]
for f in filelist:
    data=np.loadtxt(f) # Ascii data
    data=data.reshape(2*global_ny+1,2*nx+1,8)
#     kx=data[0,:,0]
#     ky=data[:,0,1]
#     jkpq_es=data[:,:,2]
#     jpqk_es=data[:,:,3]
#     jqkp_es=data[:,:,4]
#     jkpq_em=data[:,:,5]
#     jpqk_em=data[:,:,6]
#     jqkp_em=data[:,:,7]
    jkpq_es.append(data[:,:,2])
    jkpq_em.append(data[:,:,5])

jkpq_es=np.array(jkpq_es)
jkpq_em=np.array(jkpq_em)
print(jkpq_es.shape)
print(jkpq_em.shape)  


# In[ ]:


it=10
fig=plt.figure()
ax=fig.add_subplot(111)
quad=ax.pcolormesh(kx_tri,ky_tri,jkpq_es[it,:-1,:-1])
cbar=fig.colorbar(quad,shrink=1.0,aspect=5)
ax.set_title(r"Electrostatic triad transfer $J_k^{pq}$"+" (t={})".format(t_tri[it]))
ax.set_xlabel(r"Radial wavenumber $p_x$")
ax.set_ylabel(r"Poloidal wavenumber $p_y$")
plt.show()


# # Check total transfer $I_k=\sum_p\sum_q J_k^{pq}$

# In[ ]:


mx=nx
my=2
fig=plt.figure(figsize=[24,6])
ax=fig.add_subplot(111)
ax.scatter(t,neint[:,my,mx],label=r"$k_x$={},$k_y$={}".format(kx[mx],ky[my]),marker="x")
ax.scatter(t_tri,np.sum(jkpq_es,axis=(1,2)),label=r"$k_x$={},$k_y$={}".format(kx[mx],ky[my]),marker="+")
ax.set_xlabel("Time t")
ax.set_ylabel(r"Electrostatic nonlinear transfer $I_E$")
ax.set_xlim(18,22)
ax.set_ylim(-1000,500)
plt.show()


# In[ ]:





# In[ ]:




