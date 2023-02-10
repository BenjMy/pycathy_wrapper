.. pyCATHY documentation master file, created by
   sphinx-quickstart on Thu Jul 21 11:45:29 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


.. title:: Home


.. warning:: This package is in early stage of developement. The first version (0.1.0) is ready for daily use but is still changing.    
    
    
**PyCATHY** is a **Python wrapper for CATHY** :cite:p:`Camporese2010` (V1) easing the process for mesh creation, forward hydrological modeling, and output visualization of CATHY simulations.


.. hint:: The initial version (v0.1.0) includes a **Data Assimilation** extension to facilitate the integration of geophysical data into the hydrological model. You can check the list of improvements in the github releases. 


.. hint:: Beyond a scientific tool, pyCATHY aims at making the processing **interpolable**, **reusable**. 


----

.. grid:: 1 2 1 2
    :margin: 5 5 0 0
    :padding: 0 0 0 0
    :gutter: 4

    .. grid-item-card:: :octicon:`info` Getting started
        :text-align: center
        :class-title: sd-fs-5
        :class-card: sd-p-3

        New to pyCATHY? Start here!

        .. button-ref:: overview
            :ref-type: ref
            :click-parent:
            :color: primary
            :outline:
            :expand:

    .. grid-item-card:: :octicon:`comment-discussion` Need help?
        :text-align: center
        :class-title: sd-fs-5
        :class-card: sd-p-3

        Open an issue on Github.

        .. button-link:: https://github.com/BenjMy/pycathy_wrapper
            :click-parent:
            :color: primary
            :outline:
            :expand:

    .. grid-item-card:: :octicon:`file-badge` Reference documentation
        :text-align: center
        :class-title: sd-fs-5
        :class-card: sd-p-3

        A list of modules and functions.

        .. button-ref:: api
            :ref-type: ref
            :color: primary
            :outline:
            :expand:

    .. grid-item-card:: :octicon:`bookmark` Using pyCATHY for research?
        :text-align: center
        :class-title: sd-fs-5
        :class-card: sd-p-3

        Citations help support our work!

        .. button-ref:: citing
            :ref-type: ref
            :color: primary
            :outline:
            :expand:

----


.. toctree::
   :maxdepth: 1
   :caption: Getting started
   
   content/overview
   content/installing
   content/citing


.. toctree::
   :maxdepth: 3
   :caption: Gallery of examples

   content/SSHydro/index.rst   
   content/DA/index.rst


.. toctree::
   :maxdepth: 5
   :caption: Reference documentation
   
   content/api/api.rst
   content/change_log
   content/cathy_list_citation
   content/references


.. toctree::
   :maxdepth: 2
   :caption: Getting help and contributed
   
   content/contribute
   

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
