Perturbing and updating parameters by layer, zone and vegetation class
======================================================================

A data assimilation run needs an *ensemble*: many copies of the model whose parameters and
initial state differ. The ``scenario.toml`` file describes how those differences are drawn
(the *perturbation*), and which of the perturbed parameters the observations are allowed to
change (the *update*).

.. contents::
   :local:
   :depth: 2

How a scenario is written
-------------------------

A scenario is a table of lists. Every list is aligned with ``per_name``: the *i*-th entry of
each list describes the *i*-th parameter.

.. list-table::
   :header-rows: 1
   :widths: 22 20 58

   * - Key
     - Type
     - Meaning
   * - ``per_name``
     - list of str
     - Parameters to perturb (see the parameter catalogue below). A name may appear only once.
   * - ``per_nom``
     - number, or list
     - Nominal value (the value of the unperturbed model).
   * - ``per_mean``
     - number, or list
     - Mean of the distribution the ensemble is drawn from (usually the nominal value).
   * - ``per_sigma``
     - number >= 0, or list
     - Standard deviation (same unit as the parameter). 0 gives no spread.
   * - ``per_type``
     - ``"None"``, ``"additive"``, ``"multiplicative"``
     - ``"None"``: the draw *is* the parameter. ``"additive"`` / ``"multiplicative"``: the draw is
       added to / multiplies the nominal value. A single value per Feddes parameter, a list per
       entry for soil parameters controlled by layer, zone or class.
   * - ``sampling_type``
     - ``"normal"``, ``"lognormal"``, ``"uniform"``
     - Distribution of the draw.
   * - ``transf_type``
     - ``"None"``, ``"log"``, ``"log-ratio"``, ``"hyperbolic"``
     - Transform applied before the analysis and inverted after it, to keep a parameter in its
       physical range (``"log"`` for a positive one, ``"log-ratio"`` for a bounded one).
   * - ``per_bounds``
     - optional
     - Physical bounds, see *Bounds* below.
   * - ``pert_control_name``
     - optional
     - What an entry applies to: ``"layers"``, ``"zone"``, ``"root_map"``, or ``"None"`` (the
       whole domain). TOML has no null: write the string ``"None"``.
   * - ``listUpdateParm``
     - list of str
     - What the observations may change: ``"St. var."`` (the state) and the **keys** of the
       perturbed parameters (``"Ks0"``, ``"ZROOT0"``...).

What a control applies to
-------------------------

The soil table is indexed by ``(zone, layer)``. A parameter perturbed with a control becomes
one entry per layer, zone or class, and each entry writes some rows of that table.

.. list-table::
   :header-rows: 1
   :widths: 18 30 52

   * - Control
     - Keys
     - Rows of the soil table written by entry ``k``
   * - none
     - ``Ks`` (one key)
     - every row
   * - ``layers``
     - ``Ks0``, ``Ks1``...
     - layer ``k``, in every zone
   * - ``zone``
     - ``Ks0``, ``Ks1``...
     - zone number ``k + 1``, in every layer
   * - ``root_map``
     - ``Ks0``, ``Ks1``...
     - the zones that contain vegetation class ``k + 1``, in every layer

.. warning::
   The index ``k`` in a key is **0-based**, but zone and class numbers are **1-based**:
   ``Ks0`` is zone 1. Zone and class numbers must run from 1 to *K* without gaps, and the
   nodata value ``-9999`` of a map is not a zone.

.. warning::
   Entries are written in order into the same table. When a vegetation class spans several
   zones, the entries of different classes overlap and **the last one written wins**: with
   classes 1 to 3 over two zones where class 2 covers both, ``Ks0`` is overwritten by
   ``Ks1`` in zone 1, and perturbing or updating ``Ks0`` then has no effect.
   ``scenario.describe_ensemble`` lists these overlaps.

Parameter catalogue and supported layouts
-----------------------------------------

Not every parameter accepts every control. The table below was measured on the real
``perturbate()`` (``test_scenario_layouts.py``); a layout marked *fails* raises an obscure error
inside the sampler, and ``validate_scenario`` reports it by name before the run.

.. list-table::
   :header-rows: 1
   :widths: 24 19 19 38

   * - Parameter
     - Writes
     - Controls
     - Layout of ``per_nom``, ``per_mean``, ``per_sigma``
   * - ``Ks``
     - ``PERMX``, ``PERMY``, ``PERMZ`` together (isotropic)
     - none, layers, zone, root_map
     - scalar without control; a list with one entry per layer, zone or class with one
   * - ``porosity``
     - ``POROS``
     - none, layers, zone (**not** root_map)
     - as ``Ks``
   * - ``n_VG``, ``thetar_VG``, ``SATCELL_VG``
     - ``VGNCELL``, ``VGRMCCELL``, ``VGPSATCELL``
     - none, layers, zone, root_map
     - as ``Ks``
   * - ``alpha_VG``
     - nothing (no soil column is mapped to it)
     - -
     - **not usable**: it builds and samples, but writing it into the member soil files fails
       with ``Unknown parameter``; perturb ``SATCELL_VG`` instead. The validator rejects it
   * - ``ic``
     - initial pressure head
     - none, layers, zone (**not** root_map)
     - as ``Ks``; the key is ``ic`` (scalar) or ``ic0``, ``ic1``...
   * - ``WTPOSITION``
     - water table position
     - none only
     - scalar
   * - ``ZROOT``, ``PCREF``
     - Feddes table
     - ignored (rows are vegetation classes)
     - scalar if the Feddes table has **one** row, otherwise a list with **one entry per row**
   * - ``PCWLT``, ``OMGC``, ``PZ``
     - Feddes table
     - ignored
     - **scalar only**: one nominal for every row, sampled independently per row (a list fails)
   * - ``VGP``
     - nothing
     - -
     - **not implemented**: its perturbation is commented out, the validator rejects it

A list without a control (``Ks`` with ``per_nom = [1e-4, 1e-5]`` and no ``pert_control_name``)
fails for every soil parameter and for ``ic``.

**Keys.** The key of a perturbed entry is what ``listUpdateParm`` must name: ``Ks`` for a
scalar ``Ks``; ``Ks0``, ``Ks1``... for lists; ``ZROOT0``, ``ZROOT1``... for the rows of the
Feddes table (``ZROOT0`` even when there is one row). For Feddes parameters the index is the
**row position** in the Feddes table (vegetation class order, 0-based), which is *not* the
rule of the soil parameters (``veg_map == k + 1``).

Bounds
------

``per_bounds`` is optional and takes one entry per parameter: ``"None"``, a table
``{ min = a, max = b }`` (either key may be omitted, the missing side is open) or a pair
``[a, b]``. For a parameter with several entries the entry may be one specification shared by
all of them, or a list with one specification per entry.

Bounds act twice: the prior is drawn from the **truncated** distribution, and every analysed
value is **clipped** to the bounds after the update (members are kept, a warning lists the
clipped values; rejecting them would shrink the ensemble). A malformed bound, or ``min >= max``,
is an error: it used to be ignored silently, leaving the parameter unbounded.

Examples
--------

The examples assume a model with 2 layers, 2 zones and 3 vegetation classes; the line
``# model: feddes_rows=N`` gives the number of rows of the Feddes table.

**A homogeneous soil** (one value of ``Ks`` for the whole domain, one rooting depth):

.. code-block:: toml

   # model: feddes_rows=1
   [scenario.homogeneous]
   per_name = ["Ks", "ZROOT"]
   per_type = ["None", "None"]
   per_nom = [1e-4, 0.4]
   per_mean = [1e-4, 0.4]
   per_sigma = [1e-5, 0.02]
   per_bounds = [{ min = 1e-6, max = 1e-3 }, { min = 0.1, max = 1.5 }]
   sampling_type = ["normal", "normal"]
   transf_type = ["None", "None"]
   listUpdateParm = ["St. var.", "ZROOT0"]

**Two layers**, each with its own conductivity and its own bounds; only the upper layer is
updated:

.. code-block:: toml

   # model: feddes_rows=1
   [scenario.layers]
   per_name = ["Ks"]
   per_type = [["None", "None"]]
   per_nom = [[1e-4, 1e-5]]
   per_mean = [[1e-4, 1e-5]]
   per_sigma = [[1e-5, 1e-6]]
   per_bounds = [[{ min = 1e-6, max = 1e-3 }, { min = 1e-7, max = 1e-4 }]]
   sampling_type = ["normal"]
   transf_type = ["None"]
   pert_control_name = ["layers"]
   listUpdateParm = ["St. var.", "Ks0"]

**Two zones** (``Ks0`` is zone 1, ``Ks1`` is zone 2), one bound shared by both:

.. code-block:: toml

   # model: feddes_rows=1
   [scenario.zones]
   per_name = ["Ks"]
   per_type = [["None", "None"]]
   per_nom = [[1e-4, 1e-5]]
   per_mean = [[1e-4, 1e-5]]
   per_sigma = [[1e-5, 1e-6]]
   per_bounds = [{ min = 1e-7, max = 1e-3 }]
   sampling_type = ["normal"]
   transf_type = ["None"]
   pert_control_name = ["zone"]
   listUpdateParm = ["St. var.", "Ks0", "Ks1"]

**Three vegetation classes**: a soil parameter by class (``root_map``) and the rooting depth of
each class (a Feddes parameter: one list entry per row of the Feddes table):

.. code-block:: toml

   # model: feddes_rows=3
   [scenario.vegetation]
   per_name = ["Ks", "ZROOT"]
   per_type = [["None", "None", "None"], "None"]
   per_nom = [[1e-4, 1e-5, 1e-6], [0.3, 0.5, 0.7]]
   per_mean = [[1e-4, 1e-5, 1e-6], [0.3, 0.5, 0.7]]
   per_sigma = [[1e-5, 1e-6, 1e-7], [0.02, 0.02, 0.02]]
   per_bounds = ["None", { min = 0.05, max = 1.5 }]
   sampling_type = ["normal", "normal"]
   transf_type = ["None", "None"]
   pert_control_name = ["root_map", "None"]
   listUpdateParm = ["St. var.", "ZROOT0", "ZROOT1", "ZROOT2"]

Updating parameters
-------------------

``listUpdateParm`` selects what an analysis may change. ``"St. var."`` is the state (the
pressure head); every other entry must be the **key** of a parameter that was built, not a
name from ``per_name``: ``["Ks0"]`` for a scalar ``Ks`` is an error (the key is ``Ks``), and
``perturbate()`` says so, with a suggestion. A parameter that is perturbed but not listed
stays at its sampled value for the whole run.

After each analysis, in this order: the analysed parameters are transformed back
(``transf_type``), **clipped** to their bounds, checked for physical validity (a negative
rooting depth, one deeper than the mesh, a run that CATHY could not complete, rejects the
member, with the reason recorded in ``<project>_rejections.json``), and written to the input
files of each member.

Be careful with what the data can *tell*: a parameter to which the observations are
insensitive (a rooting depth for a short ERT survey of bare soil) is moved only by sampling
noise, so updating it is a source of error, not information.

Checking a scenario before running it
--------------------------------------

.. code-block:: python

   from pyCATHY.DA.scenario import validate_scenario, describe_ensemble

   errors, warnings = validate_scenario(
       scenario, nlayers=2, zone_raster=zones, veg_raster=veg, n_feddes_rows=3,
       raise_on_error=False)

   table = simu.describe_perturbation(scenario["listUpdateParm"])
   print(table[["key", "control", "target", "updated"]])
   print(table.attrs["warnings"])          # overlapping writes

``perturbate()`` runs the validator itself (``validate_scenario=False`` skips it) and checks
``listUpdateParm`` against the keys it built. ``run_DA_sequential`` logs the same table, as
``Perturbed parameters:``, at the start of a run.
