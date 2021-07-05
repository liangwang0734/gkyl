#!/usr/bin/env python
# coding: utf-8

# In[36]:


import numpy as np
import numpy as math
import matplotlib as mpl
# mpl.use('agg')
import matplotlib.pyplot as plt
import scipy.integrate as integrate


# In[37]:


import adios


# In[38]:


mu0 = 1
gasGamma = 5. / 3.

print(f'**************** Using mu0={mu0}, gasGamma={gasGamma}')


# In[39]:


filename = 'frc-equilibrium-r-theta-plane_elc_0.bp'
file = adios.file(filename)
data = file['CartGridField']

print(f'***************** Checking using fluid data in {filename}')


# In[40]:

q = file['charge'].value
m = file['mass'].value
xlo, ylo = file['lowerBounds'].value
xup, yup = file['upperBounds'].value
nx, ny = file['numCells'].value
xn = np.linspace(xlo, xup, nx + 1)
yn = np.linspace(ylo, yup, ny + 1)
xc = 0.5 * (xn[1:] + xn[:-1])
yc = 0.5 * (yn[1:] + yn[:-1])
x, y = np.meshgrid(xc, yc, indexing='ij')
r = np.sqrt(x * x + y * y)


# In[41]:


n = data[..., 0] / m
vx = data[..., 1] / n / m
vy = data[..., 2] / n / m
vz = data[..., 3] / n / m
e = data[..., 4]
p = (gasGamma - 1) * (e - 0.5 * n * m * (vx * vx + vy * vy + vz * vz))


# In[42]:


filename = 'frc-equilibrium-r-theta-plane_field_0.bp'
file = adios.file(filename)
data = file['CartGridField']

print(f'***************** Checking using field data in {filename}')

# In[43]:


Bz = data[..., 5]


# In[44]:


vr = math.sqrt(vx*vx+vy*vy)
vr.shape


# In[45]:


n_r = n[nx//2:, ny//2]
ut_r = vr[nx//2:, ny//2]
p_r = p[nx//2:, ny//2]
Bz_r = Bz[nx//2:, ny//2]
r_r = r[nx//2:, ny//2]
dr = r_r[1:] - r_r[:-1]


# In[33]:


plt.figure(figsize=(14, 4))

plt.subplot(131)
plt.plot(r_r, Bz_r, label=r'$B_z$')
plt.xlabel('r', fontsize='large')
plt.legend(fontsize='large')

plt.subplot(132)
plt.plot(r_r, -np.gradient(Bz_r, r_r) / mu0 / q, label=r'$-\partial_r Bz/\mu_0/q$')
plt.plot(r_r, ut_r, ls='--', lw=3, label=r'$u_\theta$')
plt.xlabel('r', fontsize='large')
plt.legend(fontsize='large')

plt.subplot(133)
plt.plot(r_r, np.gradient(p_r, r_r), label=r'$\partial_r p$')
plt.plot(r_r,
         n_r * q * ut_r * Bz_r - m * n_r * ut_r * ut_r / r_r,
         ls='--',
         lw=3,
         label=r'$nqu_{\theta}B_z-mnu_{\theta}^2/r$')
plt.xlabel('r', fontsize='large')
plt.legend(fontsize='large')

# plt.savefig('frc01-init.png', bbox_inches='tight')
plt.show()
