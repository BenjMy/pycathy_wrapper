.. _installing:

Installing
==========   
    
**PyCATHY** is a **Python wrapper for CATHY** :cite:p:`Camporese2010` (V1) easing the process for mesh creation, forward hydrological modeling, and output visualization of CATHY simulations.


.. hint:: On all platforms, we recommend to install pyCATHY via the conda package manager contained in the Anaconda distribution. For details on how to install Anaconda, we refer to: (https://docs.anaconda.com/anaconda/install/). 


Quick use
---------

To avoid conflicts with other packages, we recommend to install pyCATHY in a separate environment. Here we call this environment pyCATHY, but you can give it any name. Note that this environment has to be created only once.

Open a terminal (Linux & Mac) or the Anaconda Prompt (Windows) and type::

	conda create --name pyCATHY python=3.10
	conda activate pyCATHY
	pip install pyCATHY


Usage with Spyder or JupyterLab::

	conda install -c conda-forge spyder
	
Or alternatively, the web-based IDE JupyterLab (https://jupyterlab.readthedocs.io)::

	conda install -c conda-forge jupyterlab


Staying up-to-date
------------------

Install using setup.py::

    git clone https://github.com/BenjMy/pycathy_wrapper
    cd pycathy_wrapper
    python setup.py develop|install
    import pyCATHY


Dependencies
------------
- numpy, scipy, and matplotlib beyond default python packages

.. warning:: pyCATHY for Data Assimilation relies on others libraries 

   - [pyDA](https://github.com/hickmank/pyda) for the enkf and SIR algoritms


.. warning:: pyCATHY for ERT Data Assimilation relies on others libraries 

   - [pyGIMLI](https://github.com/gimli-org/gimli) 
   - OR [Resipy](https://gitlab.com/hkex/resipy)
		



How to run tests
----------------

Download one of the exemple and run python *.py script


