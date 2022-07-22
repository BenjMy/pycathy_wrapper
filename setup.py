import setuptools

with open("README.md", "r") as f:
    long_description = f.read()

setuptools.setup(
    name="pyCATHY",
    version="0.0.1",
    author="B. Mary",
    author_email="bmary@lbl.gov",
    description="Python package for managing CATHY processing, and visualization at the project level",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/BenjMy/pycathy_wrapper",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU Lesser General Public License " + "v3 (LGPLv3)",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.9',
    install_requires=["numpy", "pandas", "pyvista", "numba"]
) 

