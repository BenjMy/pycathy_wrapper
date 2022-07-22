import glob
import pyvista as pv
import os
import time

os.chdir("/home/ben/Documents/CATHY/pyCATHY/notebooks/rhizo_ET_irr_PRD/vtk/")
filename = "10*.vtk"
filename0 = "100.vtk"
mesh = pv.read(filename0)

plotter = pv.Plotter(notebook=False, off_screen=True)
plotter.add_mesh(mesh, show_edges=True)

plotter.open_gif("rhizo.gif")


legend_entry = "Time=" + str(mesh["TIME"])
plotter.show_grid()
cpos = plotter.show(interactive_update=True, auto_close=False)
# plotter.add_legend(legend_entry)
plotter.add_text(legend_entry, name="time-label")


for files in glob.glob(filename):
    mesh = pv.read(files)
    array_new = mesh.get_array("pressure")
    legend_entry = "Time=" + str(mesh["TIME"])
    plotter.update_scalars(array_new, render=False)
    plotter.add_text(legend_entry, name="time-label")

    plotter.render()
    plotter.write_frame()

plotter.close()


import imageio

gif_original = "rhizo.gif"
gif_speed_down = "rhizo_new.gif"
gif = imageio.mimread(gif_original)
imageio.mimsave(gif_speed_down, gif, fps=1)

#%%

import glob
import pyvista as pv
import os
import time

os.chdir("/home/ben/Documents/CATHY/pyCATHY/notebooks/rhizo_ET_irr_PRD/vtk/")
filename = "10*.vtk"
filename0 = "100.vtk"
mesh = pv.read(filename0)

plotter = pv.Plotter(notebook=False, off_screen=False)
plotter.add_mesh(mesh, show_edges=True)

plotter.open_gif("rhizo.gif")


legend_entry = "Time=" + str(mesh["TIME"])
plotter.show_grid()
cpos = plotter.show(interactive_update=True, auto_close=False)
# plotter.add_legend(legend_entry)
plotter.add_text(legend_entry, name="time-label")


for files in glob.glob(filename):
    mesh = pv.read(files)
    array_new = mesh.get_array("pressure")
    legend_entry = "Time=" + str(mesh["TIME"])
    plotter.update_scalars(array_new, render=True)
    plotter.add_text(legend_entry, name="time-label")

    plotter.render()
    plotter.write_frame()

plotter.close()


import imageio

gif_original = "rhizo.gif"
gif_speed_down = "rhizo_new.gif"
gif = imageio.mimread(gif_original)
imageio.mimsave(gif_speed_down, gif, fps=1)


#%%

import pyvista as pv
import numpy as np

filename = "sphere-shrinking.gif"

mesh = pv.Sphere()
mesh.cell_arrays["data"] = np.random.random(mesh.n_cells)

plotter = pv.Plotter()
# Open a movie file
# plotter.open_movie(filename)
plotter.open_gif(filename)

# Add initial mesh
plotter.add_mesh(mesh, scalars="data", clim=[0, 1])
# Add outline for shrinking reference
plotter.add_mesh(mesh.outline_corners())

plotter.show(auto_close=False)  # only necessary for an off-screen movie

# Run through each frame
plotter.write_frame()  # write initial data

# Update scalars on each frame
for i in range(100):
    random_points = np.random.random(mesh.points.shape)
    mesh.points = random_points * 0.01 + mesh.points * 0.99
    mesh.points -= mesh.points.mean(0)
    mesh.cell_arrays["data"] = np.random.random(mesh.n_cells)
    # plotter.add_text(f"Iteration: {i}", name='time-label')
    plotter.add_text(str(i), name="time-label")
    plotter.write_frame()  # Write this frame

# Be sure to close the plotter when finished
plotter.close()
