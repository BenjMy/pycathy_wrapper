[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1Zl-VpMbrESu9MbeNpgUvXS6nfluno_dG#scrollTo=be0ef8fe-51bb-4d86-b00a-fa027c286ecf)

## Python package for hydrogeophysical data modeling

pyCATHY is a **Python wrapper for CATHY** (V1, and plant model version) allowing mesh creation, forward and inverse modeling, and simple output visualization of CATHY simulations

### Dependencies

pyCATHY for Data Assimilation relies on others libraries: 
- [pyGIMLI](https://github.com/gimli-org/gimli) or [Resipy](https://gitlab.com/hkex/resipy) for the ERT simulation
- [pyDA](https://github.com/hickmank/pyda) for the enkf and SIR algoritms


### Installation and testing

* Install using setup.py
* Or with pip: "pip install pyCATHY" [**NOT YET DEPLOYED**]
* Dependencies: numpy, scipy, and matplotlib beyond default python packages
* How to run tests: Test with ??\??*.py scripts

---

### Licence ###
This is free software: you can redistribute it and/or modify it under the terms of the BSD 3-clause License. A copy of this license is provided in LICENSE.md.

### References ###

- Camporese, M., Paniconi, C., Putti, M., Orlandini, S., 2010. Surface-subsurface flow modeling with path-based runoff routing, boundary condition-based coupling, and assimilation of multisource observation data: SURFACE-SUBSURFACE FLOW MODELING. Water Resour. Res. 46. https://doi.org/10.1029/2008WR007536



