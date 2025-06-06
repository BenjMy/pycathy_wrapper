
.. DO NOT EDIT.
.. THIS FILE WAS AUTOMATICALLY GENERATED BY SPHINX-GALLERY.
.. TO MAKE CHANGES, EDIT THE SOURCE PYTHON FILE:
.. "content/SSHydro/plot_4b_pyCATHY_outputs.py"
.. LINE NUMBERS ARE GIVEN BELOW.

.. only:: html

    .. note::
        :class: sphx-glr-download-link-note

        :ref:`Go to the end <sphx_glr_download_content_SSHydro_plot_4b_pyCATHY_outputs.py>`
        to download the full example code.

.. rst-class:: sphx-glr-example-title

.. _sphx_glr_content_SSHydro_plot_4b_pyCATHY_outputs.py:


Output plots part 2
===================

Weill, S., et al. « Coupling Water Flow and Solute Transport into a Physically-Based Surface–Subsurface Hydrological Model ». 
Advances in Water Resources, vol. 34, no 1, janvier 2011, p. 128‑36. DOI.org (Crossref), 
https://doi.org/10.1016/j.advwatres.2010.10.001.

This example shows how to use pyCATHY object to plot the most common ouputs of the hydrological model.

*Estimated time to run the notebook = 5min*

.. GENERATED FROM PYTHON SOURCE LINES 17-19

Here we need to import `cathy_tools` class that control the CATHY core files preprocessing and processing
We also import `cathy_plots` to render the results

.. GENERATED FROM PYTHON SOURCE LINES 19-25

.. code-block:: Python


    from pyCATHY import cathy_tools
    from pyCATHY.plotters import cathy_plots as cplt










.. GENERATED FROM PYTHON SOURCE LINES 26-27

if you add True to verbose, the processor log will be printed in the window shell

.. GENERATED FROM PYTHON SOURCE LINES 27-45

.. code-block:: Python




    path2prj = "../SSHydro/"  # add your local path here
    simu = cathy_tools.CATHY(dirName=path2prj, 
    			prj_name="weil_exemple_outputs_plot"
    			)

    simu.run_preprocessor()
    simu.run_processor(IPRT1=2, 
                        DTMIN=1e-2,
                        DTMAX=1e2,
                        DELTAT=5,
                       TRAFLAG=0,
                       verbose=False
                       )






.. rst-class:: sphx-glr-script-out

 .. code-block:: none

    🏁 Initiate CATHY object
    🍳 gfortran compilation
    👟 Run preprocessor
    🔄 Update parm file 
    🔄 Update hap.in file
    🔄 Update dem_parameters file 
    🔄 Update dem_parameters file 
    🛠  Recompile src files [3s]
    🍳 gfortran compilation [7s]
    b''
    👟 Run processor




.. GENERATED FROM PYTHON SOURCE LINES 46-48

.. code-block:: Python

    simu.show(prop="hgsfdet")




.. image-sg:: /content/SSHydro/images/sphx_glr_plot_4b_pyCATHY_outputs_001.png
   :alt: plot 4b pyCATHY outputs
   :srcset: /content/SSHydro/images/sphx_glr_plot_4b_pyCATHY_outputs_001.png
   :class: sphx-glr-single-img





.. GENERATED FROM PYTHON SOURCE LINES 49-51

.. code-block:: Python

    simu.show(prop="dtcoupling", yprop="Atmpot-d")




.. image-sg:: /content/SSHydro/images/sphx_glr_plot_4b_pyCATHY_outputs_002.png
   :alt: plot 4b pyCATHY outputs
   :srcset: /content/SSHydro/images/sphx_glr_plot_4b_pyCATHY_outputs_002.png
   :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 .. code-block:: none

    /home/z0272571a@CAMPUS.CSIC.ES/Nextcloud/BenCSIC/Codes/BenjMy/pycathy_wrapper/pyCATHY/importers/cathy_outputs.py:322: UserWarning: Input line 3 contained no data and will not be counted towards `max_rows=236`.  This differs from the behaviour in NumPy <=1.22 which counted lines rather than rows.  If desired, the previous behaviour can be achieved by using `itertools.islice`.
    Please see the 1.23 release notes for an example on how to do this.  If you wish to ignore this warning, use `warnings.filterwarnings`.  This warning is expected to be removed in the future and is given only once per `loadtxt` call.
      dtcoupling = np.loadtxt(dtcoupling_file, skiprows=2, max_rows=2 + nstep)




.. GENERATED FROM PYTHON SOURCE LINES 52-54

.. code-block:: Python

    simu.show(prop="hgraph")




.. image-sg:: /content/SSHydro/images/sphx_glr_plot_4b_pyCATHY_outputs_003.png
   :alt: plot 4b pyCATHY outputs
   :srcset: /content/SSHydro/images/sphx_glr_plot_4b_pyCATHY_outputs_003.png
   :class: sphx-glr-single-img





.. GENERATED FROM PYTHON SOURCE LINES 55-57

.. code-block:: Python

    simu.show(prop="cumflowvol")




.. image-sg:: /content/SSHydro/images/sphx_glr_plot_4b_pyCATHY_outputs_004.png
   :alt: Cumulative flow volume
   :srcset: /content/SSHydro/images/sphx_glr_plot_4b_pyCATHY_outputs_004.png
   :class: sphx-glr-single-img





.. GENERATED FROM PYTHON SOURCE LINES 58-59

To select another time step change the value in the function argument

.. GENERATED FROM PYTHON SOURCE LINES 59-66

.. code-block:: Python

    cplt.show_vtk(
        unit="pressure",
        timeStep=1,
        notebook=False,
        path=simu.workdir + "/weil_exemple_outputs_plot/vtk/",
    )





.. image-sg:: /content/SSHydro/images/sphx_glr_plot_4b_pyCATHY_outputs_005.png
   :alt: plot 4b pyCATHY outputs
   :srcset: /content/SSHydro/images/sphx_glr_plot_4b_pyCATHY_outputs_005.png
   :class: sphx-glr-single-img




.. rst-class:: sphx-glr-script-out

 .. code-block:: none

    plot pressure




.. GENERATED FROM PYTHON SOURCE LINES 67-73

cplt.show_vtk(
    unit="saturation",
    timeStep=1,
    notebook=False,
    path=simu.workdir + "/my_cathy_prj/vtk/",
)


.. rst-class:: sphx-glr-timing

   **Total running time of the script:** (0 minutes 23.141 seconds)


.. _sphx_glr_download_content_SSHydro_plot_4b_pyCATHY_outputs.py:

.. only:: html

  .. container:: sphx-glr-footer sphx-glr-footer-example

    .. container:: sphx-glr-download sphx-glr-download-jupyter

      :download:`Download Jupyter notebook: plot_4b_pyCATHY_outputs.ipynb <plot_4b_pyCATHY_outputs.ipynb>`

    .. container:: sphx-glr-download sphx-glr-download-python

      :download:`Download Python source code: plot_4b_pyCATHY_outputs.py <plot_4b_pyCATHY_outputs.py>`

    .. container:: sphx-glr-download sphx-glr-download-zip

      :download:`Download zipped: plot_4b_pyCATHY_outputs.zip <plot_4b_pyCATHY_outputs.zip>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
