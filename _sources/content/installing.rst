.. _installing:

Installing
==========   
    
**PyCATHY** is a **Python wrapper for CATHY** :cite:p:`Camporese2010` (V1) easing the process for mesh creation, forward hydrological modeling, and output visualization of CATHY simulations.


.. hint:: On all platforms, we recommend to install pyCATHY via the conda package manager contained in the Anaconda distribution. For details on how to install Anaconda, we refer to: (https://docs.anaconda.com/anaconda/install/). 


Set up the environment
----------------------

To avoid conflicts with other packages, we recommend to install pyCATHY in a separate environment. Here we call this environment pyCATHY, but you can give it any name. Note that this environment has to be created only once.

Open a terminal (Linux & Mac) or the Anaconda Prompt (Windows) and type::

	conda create --name pyCATHY python=3.10
	conda activate pyCATHY

Usage with Spyder or JupyterLab::

	conda install -c conda-forge spyder
	
Or alternatively, the web-based IDE JupyterLab (https://jupyterlab.readthedocs.io)::

	conda install -c conda-forge jupyterlab

	
Install pyCATHY
---------------

Pin to the latest released pip version::

    pip install pyCATHY
	
Staying up-to-date using setup.py::

    git clone https://github.com/BenjMy/pycathy_wrapper
    cd pycathy_wrapper
    python setup.py develop|install


.. hint::  pyCATHY already includes CATHY cores, and thus does not require a separate installation of CATHY before using pyCATHY
 
Dependencies
------------
- NumPy, SciPy, and Matplotlib beyond default python packages

- Fortran and the algebra libraries BLAS and LAPACK.

.. tab-set::




    .. tab-item:: Ubuntu 22

        .. code:: bash

           sudo apt-get update
           sudo apt-get install gfortran
           sudo apt install libopenblas-dev
           
    .. tab-item:: Ubuntu 20

        .. code:: bash

           sudo apt-get update
           sudo apt-get install gfortran
           sudo apt-get install blas-dev lapack-dev

    .. tab-item:: Fedora

        .. code:: bash

            sudo dnf upgrade --refresh
            sudo dnf install gcc-gfortran
	    sudo dnf install blas-devel
	    sudo dnf install lapack-devel


    .. tab-item:: Windows

        We recommand to install a **virtual machine**:

        .. code:: bash

            install WSL from Microsoft Store app (search WSL)
	    install linux from Microsoft Store app (search ubuntu lts)
	    turn Windows features on or off:
	      - enable Windows subsystem for Linux
	      - enable Virtual machine platform
	     restart PC 


.. warning:: pyCATHY for Data Assimilation relies on others libraries 

   - [pyDA](https://github.com/hickmank/pyda) for the enkf and SIR algoritms


.. warning:: pyCATHY for ERT Data Assimilation relies on others libraries 

   - [pyGIMLI](https://github.com/gimli-org/gimli) 
   - OR [Resipy](https://gitlab.com/hkex/resipy)
		


How to run tests
----------------

Start jupyterlab, select and exemple and run it. 

OR 

Download one of the exemple and **run in the terminal**:

    conda activate pyCATHY
    python myexemple.py


