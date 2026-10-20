Localisation
============

With a small ensemble, the covariance between a state far from an observation and that
observation is mostly sampling noise. *Localisation* limits an observation to the states it
can plausibly inform. pyCATHY has two unrelated mechanisms; choose by the kind of data.

.. contents::
   :local:
   :depth: 2

Which mechanism for which data
------------------------------

.. list-table::
   :header-rows: 1
   :widths: 26 26 48

   * - Data
     - Domain localisation (``localisation=``)
     - Covariance localisation (matrix)
   * - maps of soil moisture, ET, EM (observed on every surface node)
     - yes: ``"nodes"``, ``"zones"``, ``"veg_map"``
     - yes
   * - **ERT** (apparent resistivity)
     - **no**: raises ``Cannot localise with apparent resistivity data.``
     - yes, the only option
   * - point sensors (tensiometers, a few probes)
     - **no**: fails in ``_DA_analysis``
     - yes (one column per sensor)

Domain localisation
-------------------

``AnalysisConfig(localisation=...)`` (or ``DAConfig``) takes one of ``"nodes"``, ``"zones"``,
``"veg_map"``; a list with more than one entry is an error. The domain is cut into groups of
surface nodes, one analysis is run per group, and each analysis uses only the observations of
its group.

.. list-table::
   :header-rows: 1
   :widths: 14 28 30 28

   * - Mode
     - Group
     - Updated in the group
     - Notes
   * - ``"nodes"``
     - one surface node
     - the state of the column below it; **all** parameters of ``listUpdateParm`` (the code logs
       ``loc_shortcut``)
     - one analysis per node: expensive. How the parameter results of the many analyses are
       combined has **not been verified**
   * - ``"zones"``
     - one zone of the zone map
     - the state and the parameter of the zone
     - zone numbers 1 to *K*, nodata excluded
   * - ``"veg_map"``
     - one vegetation class
     - the state and the parameter of the class (for instance ``ZROOT0`` for class 1)
     - class numbers 1 to *K*, nodata excluded

Things to know:

* A surface observation is extended to **every depth** of its column: the group of a zone
  contains all the nodes below the surface nodes of the zone. This is natural for soil
  moisture maps and wrong for data whose sensitivity changes with depth.
* A zone or class without any observation is skipped silently; nodes that no group analysed keep
  their forecast state (a warning is logged).
* **Pairing of parameters.** With ``"zones"`` or ``"veg_map"`` the parameter ``<name><k>`` of
  ``listUpdateParm`` belongs to zone or class number ``k + 1`` (``ZROOT0`` to class 1), whatever
  its position in the list, and the state entries (``"St. var."``) are kept for every zone. A
  zone with no parameter of its own updates the state only. Names without a numeric suffix
  fall back to the old positional rule (the parameter at position ``zone + 1``, after a leading
  state entry).

.. code-block:: python

   from pyCATHY.DA.cathy_DA import DAConfig, AnalysisConfig

   config = DAConfig(analysis=AnalysisConfig(localisation="veg_map"))
   simu.run_DA_sequential(
       parallel=True, dict_obs=dict_obs, list_assimilated_obs="all",
       list_parm2update=["St. var.", "ZROOT0", "ZROOT1", "ZROOT2"],
       DA_type="enkf_Evensen2009_Sakov", dict_parm_pert=var_per_dict_stacked,
       config=config,
   )

Covariance localisation (a matrix)
----------------------------------

The covariance of each state with each observation is multiplied element by element by a
weight in [0, 1] that decreases with their distance. It is available with
``DA_type="enkf_analysis_localized_with_inflation"``, and it acts on the **states only**: the
parameters are not localised.

*The weight.* ``localisation.gaspari_cohn(r, L)`` is the compactly supported function of
Gaspari and Cohn (1999). It is 1 at ``r = 0``, decreases monotonically, and is **exactly zero
beyond ``2 L``**: ``L`` is half the support, not the cut-off radius. ``L`` must be positive.

*The matrix.* ``localisation.build_localization_matrix(obs_pos, grid_pos, L)`` returns an array
of shape ``(n_grid, n_obs)``: ``loc[i, j]`` weights state ``i`` with observation ``j``. Positions
are the rows of ``obs_pos`` and ``grid_pos``, in any column order provided both use the same
one (``grid3d['mesh3d_nodes']`` is ordered ``[y, x, z]``). Add the ``z`` column to localise in
depth.

*Giving it to a run.* Either one matrix per DA cycle in ``simu.localisationMatrix`` (a list
indexed by cycle), or ``AnalysisConfig(localisation_matrix=...)`` for all cycles.

.. code-block:: python

   import numpy as np
   from pyCATHY.DA import localisation

   nodes = simu.grid3d["mesh3d_nodes"]                 # (n_nodes, 3), [y, x, z]
   # one position per observation: for ERT, for instance, the centre of the four electrodes
   obs_pos = quadrupole_centres                         # (n_obs, 3) in the same [y, x, z] order
   loc = localisation.build_localization_matrix(obs_pos, nodes, L=2.0)   # (n_nodes, n_obs)
   simu.localisationMatrix = [loc] * n_cycles

**With ERT.** The quadrupoles flagged unusable (see ``mask_bad_ert``) are removed from the
analysis, and so are their columns of the matrix: give the matrix for all the quadrupoles of
the survey and it is cut to the kept ones; a matrix that already has the width of the kept
data is used as is; any other width logs a warning and is rejected by the filter. Before this
was handled, one excluded quadrupole made the localised filter fail with a shape error.

Choosing
--------

* Soil moisture or ET maps, several soil units or vegetation classes with their own parameter:
  ``"zones"`` or ``"veg_map"``.
* Maps with spatially varying parameters and a large ensemble budget: ``"nodes"``.
* ERT, or any data that is not a surface map: a covariance localisation matrix, with ``z`` in
  the positions when the sensitivity depends on depth.
* Few data and a homogeneous domain: no localisation is needed.
