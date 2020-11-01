#!/usr/bin/env python
# coding: utf-8

# # Plot mominz*

# In[ ]:


import numpy as np
import matplotlib.pyplot as plt
import glob
import f90nml
from read_f90 import read_parameters, read_time # In-house module for GKV-diag outputs

### GKV parameters from gkvp_header.f90 ###
global_nz=read_parameters("../src/gkvp_header.f90", "global_nz", int)

### GKV parameters from gkvp_namelist ###
nml=f90nml.read("../gkvp_namelist.001")
#print(nml)
n_tht=nml["nperi"]["n_tht"]
dtout_ptn=nml["times"]["dtout_ptn"]

print("global_nz =",global_nz)
print("n_tht =",n_tht)
print("dtout_ptn =",dtout_ptn)


# In[ ]:


### Load data of mominz* ###
#filelist=sorted(glob.glob("./data/mominz_mx0000my0001s0_t*dat"))
filelist=sorted(glob.glob("./data/mominz_connect_mx0000my0001s0_t*dat"))

# Time
t=[]
for f in filelist:
    wt=read_time(f)
    t.append(wt)
t=np.array(t)
print("t =",t)

# Coordinates: zz
data=np.loadtxt(filelist[0]) # Ascii data
zz=data[:,0]
print("zz =",zz)

# Values from mominz*
dens=[]
for f in filelist:
    data=np.loadtxt(f) # Ascii data
#     zz=data[:,0]
    re_dens=data[:,1]
    im_dens=data[:,2]
#     re_upara=data[:,3]
#     im_upara=data[:,4]
#     re_ppara=data[:,5]
#     im_ppara=data[:,6]
#     re_pperp=data[:,7]
#     im_pperp=data[:,8]
#     re_qlpara=data[:,9]
#     im_qlpara=data[:,10]
#     re_qlperp=data[:,11]
#     im_qlperp=data[:,12]
    dens.append(re_dens+1j*im_dens)
    
dens=np.array(dens)
print(dens.shape)  


# In[ ]:


it=10
fig=plt.figure()
ax=fig.add_subplot(111)
ax.plot(zz,dens[it,:].real,label=r"Re[$n_k$]")
ax.plot(zz,dens[it,:].imag,label=r"Im[$n_k$]")
ax.set_title(r"$n_k$ (t={})".format(t[it]))
ax.set_xlabel(r"Field-aligned coordinate $z$")
ax.legend()
ax.grid(ls="--")
plt.show()

# Normalized at flux-tube center zz=0
iz0=int(len(zz)/2)
dens0=dens[:,iz0]
dens_normalized=dens/dens0.reshape(len(dens0),1)

it=10
fig=plt.figure()
ax=fig.add_subplot(111)
ax.plot(zz,dens_normalized[it,:].real,label=r"Re[$n_k$]")
ax.plot(zz,dens_normalized[it,:].imag,label=r"Im[$n_k$]")
ax.set_title(r"$n_k$ (t={})".format(t[it]))
ax.set_xlabel(r"Field-aligned coordinate $z$")
ax.legend()
ax.grid(ls="--")
plt.show()


# In[ ]:


### Example of animation ###
from matplotlib.animation import FuncAnimation
from IPython.display import HTML

fig=plt.figure()
ax=fig.add_subplot(111)

def update_ax(i):
    ax.clear()
    ax.plot(zz,dens[i,:].real,label=r"Re[$n_k$]")
    ax.plot(zz,dens[i,:].imag,label=r"Im[$n_k$]")
    ax.set_title(r"$n_k$ (t={})".format(t[i]))
    ax.set_xlabel(r"Field-aligned coordinate $z$")
    ax.legend()
    ax.grid(ls="--")
    return

ani = FuncAnimation(fig, update_ax,
                    frames=range(0,len(t),10), interval=100)
#ani.save('advection.mp4', writer="ffmpeg", dpi=100)
#plt.close()

HTML(ani.to_jshtml())


# In[ ]:





# In[ ]:




