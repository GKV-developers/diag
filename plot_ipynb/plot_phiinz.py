#!/usr/bin/env python
# coding: utf-8

# # Plot phiinz*
# Alinz* also has the same format.

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


### Load data of phiinz* ###
#filelist=sorted(glob.glob("./data/phiinz_mx0000my0001_t*dat"))
filelist=sorted(glob.glob("./data/phiinz_connect_mx0000my0001_t*dat"))

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

# Values from phiinz*
phi=[]
for f in filelist:
    data=np.loadtxt(f) # Ascii data
#     zz=data[:,0]
    re_phi=data[:,1]
    im_phi=data[:,2]
    phi.append(re_phi+1j*im_phi)
    
phi=np.array(phi)
print(phi.shape)  


# In[ ]:


it=10
fig=plt.figure()
ax=fig.add_subplot(111)
ax.plot(zz,phi[it,:].real,label=r"Re[$\phi_k$]")
ax.plot(zz,phi[it,:].imag,label=r"Im[$\phi_k$]")
ax.set_title(r"$\phi_k$ (t={})".format(t[it]))
ax.set_xlabel(r"Field-aligned coordinate $z$")
ax.legend()
ax.grid(ls="--")
plt.show()

# Normalized at flux-tube center zz=0
iz0=int(len(zz)/2)
phi0=phi[:,iz0]
phi_normalized=phi/phi0.reshape(len(phi0),1)

it=10
fig=plt.figure()
ax=fig.add_subplot(111)
ax.plot(zz,phi_normalized[it,:].real,label=r"Re[$\phi_k$]")
ax.plot(zz,phi_normalized[it,:].imag,label=r"Im[$\phi_k$]")
ax.set_title(r"$\phi_k$ (t={})".format(t[it]))
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
    ax.plot(zz,phi[i,:].real,label=r"Re[$\phi_k$]")
    ax.plot(zz,phi[i,:].imag,label=r"Im[$\phi_k$]")
    ax.set_title(r"$\phi_k$ (t={})".format(t[i]))
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




