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

- A Fortran compiler and the algebra libraries BLAS and LAPACK. **They are installed in the
  conda/mamba environment: nothing has to be installed system-wide** (no ``apt install gfortran``,
  no WSL on Windows).

Linux, macOS and Windows (conda or mamba)::

    conda create --name pyCATHY python=3.10
    conda activate pyCATHY
    conda install -c conda-forge fortran-compiler libblas liblapack
    pip install pyCATHY

or, from the repository, ``conda env create -f environment.yml``.

pyCATHY looks for the compiler in this order: ``$PYCATHY_FC``, ``$FC``, the gfortran of the active
conda environment, then ``gfortran`` / ``flang`` on the ``PATH``. The compiled CATHY executables are
cached (content hash of the Fortran sources, ``CATHY.H``, compiler and flags): changing input files
such as the vegetation parameters does not trigger a recompilation. Set ``$PYCATHY_CACHE`` to choose the
cache directory and ``$PYCATHY_FFLAGS`` to replace the default compiler flags.

.. tab-set::

    .. tab-item:: System compiler (advanced)

        If you prefer a compiler installed by your system, make sure ``gfortran``, BLAS and LAPACK are
        on the path of the linker, for example:

        .. code:: bash

           sudo apt-get install gfortran libopenblas-dev     # Ubuntu 22
           sudo dnf install gcc-gfortran blas-devel lapack-devel   # Fedora

    .. tab-item:: Windows

        Install Miniforge/Anaconda for Windows and follow the conda instructions above (native
        Windows support is being validated, see ``REFACTOR_GFORTRAN.md``).


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


