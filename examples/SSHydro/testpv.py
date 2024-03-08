#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test pv plot
============================================
test
"""

#%%
import pyvista

sphere = pyvista.Sphere()
out = sphere.plot() 

#%%
import pyvista as pv
pv.set_jupyter_backend('client')
pv.Cone().plot()
pl.show()

#%%
import pyvista as pv
pv.set_jupyter_backend('trame')

pl = pv.Plotter()
pl.add_mesh(pv.Cone())
pl.show()
