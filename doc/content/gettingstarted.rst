.. _gettingstarted:

Getting Started
===============
    

.. image:: img/Francesca_instructions.png
  :width: 400
  :alt: Getting started
  
  
  
Create a project 
----------------

This will download automatically the source file from the gitbucket.

.. code:: bash

   from pyCATHY import cathy_tools
   from pyCATHY.plotters import cathy_plots as cplt
   simu = cathy_tools.CATHY(dirName=path2prj)


Update pre-processing files 
---------------------------

.. code:: bash

   simu.update_prepo_inputs()

Run preprocessor
----------------

.. code:: bash

   simu.run_preprocessor(verbose=True)
   

Run processor
-------------

.. code:: bash

   simu.run_processor(verbose=True)

Show outputs
------------

.. code:: bash

   cplt.show_vtk(unit="pressure", timeStep=1, notebook=False,
              path="./my_cathy_prj/vtk/")
   
   
