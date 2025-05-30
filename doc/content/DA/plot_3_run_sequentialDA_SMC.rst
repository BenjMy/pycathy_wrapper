
.. DO NOT EDIT.
.. THIS FILE WAS AUTOMATICALLY GENERATED BY SPHINX-GALLERY.
.. TO MAKE CHANGES, EDIT THE SOURCE PYTHON FILE:
.. "content/DA/plot_3_run_sequentialDA_SMC.py"
.. LINE NUMBERS ARE GIVEN BELOW.

.. only:: html

    .. note::
        :class: sphx-glr-download-link-note

        :ref:`Go to the end <sphx_glr_download_content_DA_plot_3_run_sequentialDA_SMC.py>`
        to download the full example code.

.. rst-class:: sphx-glr-example-title

.. _sphx_glr_content_DA_plot_3_run_sequentialDA_SMC.py:


Read SMC sensors observations to assimilate
===========================================

The notebook illustrate how to read SMC sensors dataset to be prepare for DA

*Estimated time to run the notebook = 2min*

.. GENERATED FROM PYTHON SOURCE LINES 10-17

.. code-block:: Python

    import numpy as np
    from pyCATHY.DA.cathy_DA import DA
    from pyCATHY.DA.observations import make_data_cov
    from pyCATHY.DA.cathy_DA import DA, dictObs_2pd
    from pyCATHY.DA import perturbate
    import pickle








.. GENERATED FROM PYTHON SOURCE LINES 18-19

-----------------------

.. GENERATED FROM PYTHON SOURCE LINES 19-26

.. code-block:: Python

    simuWithDA = DA(
                    dirName='./DA_with_swc',
                    prj_name='DA_SMC',
                    notebook=True,
                    )






.. rst-class:: sphx-glr-script-out

 .. code-block:: none

    🏁 Initiate CATHY object




.. GENERATED FROM PYTHON SOURCE LINES 27-34

.. code-block:: Python

    abs_data_err = 1e-1 # constant error does not vary with time
    dict_obs = {} # initiate the dictionnary

    with open('./DA_with_swc/obs_prepared_SMC.pkl', 'rb') as fp:
        dict_obs = pickle.load(fp)
    data_measure_df = dictObs_2pd(dict_obs)
        # data_measure_df.index







.. GENERATED FROM PYTHON SOURCE LINES 35-37

By default, there is no correlation between sensors
Therefore, the covariance matrices are diagonal with the error values on the diagonals

.. GENERATED FROM PYTHON SOURCE LINES 37-47

.. code-block:: Python


    _,_, stacked_data_cov = make_data_cov(
                                            simuWithDA,
                                            dict_obs,
                                            list_assimilated_obs = 'all',
                                            )
    print(np.shape(stacked_data_cov))
    simuWithDA.stacked_data_cov = stacked_data_cov
    # print(np.shape(simuWithDA.stacked_data_cov))





.. rst-class:: sphx-glr-script-out

 .. code-block:: none

    (49, 3, 3)




.. GENERATED FROM PYTHON SOURCE LINES 48-84

.. code-block:: Python

    DEM, _ = simuWithDA.read_inputs('dem')
    simuWithDA.DEM = DEM
    simuWithDA.update_dem_parameters()
    simuWithDA.update_veg_map()
    simuWithDA.update_soil()

    NENS = 5

    # ZROOT
    # -------------------
    pert_nom_ZROOT = 1
    pert_sigma_ZROOT = 0.35e-9
    minZROOT = 0
    maxZROOT = 2

    scenario = {'per_type': [None],
                 'per_name':['ZROOT'],
                 'per_nom':[pert_nom_ZROOT],
                 'per_mean':[pert_nom_ZROOT],
                 'per_sigma': [pert_sigma_ZROOT],
                 'per_bounds': [
                                {'min':minZROOT,'max':maxZROOT}
                                ],
                 'sampling_type': ['normal'],
                 'transf_type':[None],
                 'listUpdateParm': ['St. var.', 'ZROOT'],
                 'listObAss': ['SMC'],
                 }

    scenario['per_name']

    list_pert = perturbate.perturbate(simuWithDA,
                                      scenario,
                                      NENS
                                      )





.. rst-class:: sphx-glr-script-out

 .. code-block:: none

    🔄 Update dem_parameters file 
    🔄 Update hap.in file
    🔄 Update dem_parameters file 
    🔄 Update dem_parameters file 
    ─────────────────────────────────────────────────────────────────────────────── ⚠ warning messages above ⚠ ───────────────────────────────────────────────────────────────────────────────

                                The parm dictionnary is empty
                                Falling back to defaults to update CATHYH
                                This can have consequences !!
                            
    ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    🔄 Update parm file 
    🔄 Update soil
    homogeneous soil




.. GENERATED FROM PYTHON SOURCE LINES 85-86

stop

.. GENERATED FROM PYTHON SOURCE LINES 86-111

.. code-block:: Python

    import os

    var_per_dict_stacked = {}
    for dp in list_pert:
        savefig = os.path.join(
                                simuWithDA.workdir,
                                simuWithDA.project_name,
                                simuWithDA.project_name + dp['savefig']
                                )
        np.random.seed(1)
        # need to call perturbate_var as many times as variable to perturbate
        # return a dict merging all variable perturbate to parse into prepare_DA
        var_per_dict_stacked = perturbate.perturbate_parm(
                                                        var_per_dict_stacked,
                                                        parm=dp,
                                                        type_parm = dp['type_parm'], # can also be VAN GENUCHTEN PARAMETERS
                                                        mean =  dp['mean'],
                                                        sd =  dp['sd'],
                                                        sampling_type =  dp['sampling_type'],
                                                        ensemble_size =  dp['ensemble_size'], # size of the ensemble
                                                        per_type= dp['per_type'],
                                                        savefig=savefig
                                                        )





.. image-sg:: /content/DA/images/sphx_glr_plot_3_run_sequentialDA_SMC_001.png
   :alt: Histogram of ZROOT0
   :srcset: /content/DA/images/sphx_glr_plot_3_run_sequentialDA_SMC_001.png
   :class: sphx-glr-single-img





.. GENERATED FROM PYTHON SOURCE LINES 112-115

f
simuWithDA.parm
simuWithDA.read_inputs('atmbc')

.. GENERATED FROM PYTHON SOURCE LINES 115-125

.. code-block:: Python

    atmbc_times = data_measure_df.index.get_level_values(1).unique().to_list()
    simuWithDA.update_atmbc(HSPATM=1,IETO=0,
                            time=atmbc_times,
                            netValue=[0]*len(atmbc_times)
                            )

    # simuWithDA.update_parm()
    # simuWithDA.read_inputs('atmbc')






.. rst-class:: sphx-glr-script-out

 .. code-block:: none

    🔄 Update atmbc
    🔄 Update parm file 




.. GENERATED FROM PYTHON SOURCE LINES 126-127

simuWithDA.atmbc

.. GENERATED FROM PYTHON SOURCE LINES 127-138

.. code-block:: Python


    # simuWithDA.run_DA_smooth(
    #                           VTKF=2,
    #                           TRAFLAG=0,
    #                           dict_obs= dict_obs,
    #                           list_assimilated_obs='all', # default
    #                           list_parm2update= ['St. var.', 'ZROOT0'],
    #                           DA_type='enkf_Evensen2009',
    #                           dict_parm_pert=var_per_dict_stacked,
    #                         )








.. GENERATED FROM PYTHON SOURCE LINES 139-149

.. code-block:: Python


    # simuWithDA.run_DA_sequential(
    #                               VTKF=2,
    #                               TRAFLAG=0,
    #                               dict_obs= dict_obs,
    #                               list_assimilated_obs='all', # default
    #                               list_parm2update= ['St. var.', 'ZROOT0'],
    #                               DA_type='enkf_Evensen2009',
    #                               dict_parm_pert=var_per_dict_stacked,
    #                             )








.. rst-class:: sphx-glr-timing

   **Total running time of the script:** (0 minutes 0.184 seconds)


.. _sphx_glr_download_content_DA_plot_3_run_sequentialDA_SMC.py:

.. only:: html

  .. container:: sphx-glr-footer sphx-glr-footer-example

    .. container:: sphx-glr-download sphx-glr-download-jupyter

      :download:`Download Jupyter notebook: plot_3_run_sequentialDA_SMC.ipynb <plot_3_run_sequentialDA_SMC.ipynb>`

    .. container:: sphx-glr-download sphx-glr-download-python

      :download:`Download Python source code: plot_3_run_sequentialDA_SMC.py <plot_3_run_sequentialDA_SMC.py>`

    .. container:: sphx-glr-download sphx-glr-download-zip

      :download:`Download zipped: plot_3_run_sequentialDA_SMC.zip <plot_3_run_sequentialDA_SMC.zip>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
