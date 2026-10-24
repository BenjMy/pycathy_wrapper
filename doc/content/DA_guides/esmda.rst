ES-MDA (smoother)
=================

*Ensemble Smoother with Multiple Data Assimilation* (Emerick & Reynolds, 2013) estimates
**parameters** from all the observations of a period at once. Unlike the sequential filters
(EnKF, PF), which update cycle after cycle and carry the analysed state forward, ES-MDA
assimilates the same data ``n_iter`` times and **re-runs the whole simulation from t = 0**
between iterations. It suits strongly non-linear parameter problems, where one Kalman update
is not enough, and it never needs to restart CATHY from an analysed (possibly non-physical)
state.

.. contents::
   :local:
   :depth: 2

.. warning::

   The CATHY adapter (``DA.run_DA_esmda``) is **experimental**. It is unit tested against a
   stub of the simulation object but has not been validated on a real CATHY run yet: start
   with few members and a short period, and read `What to check on a first run`_.
   The numerics (``pyCATHY.DA.esmda``) are validated against the analytical Kalman posterior.

The method
----------

For member :math:`j` at iteration :math:`i = 1 \dots N_a`:

.. math::

   d_j = d_{obs} + \sqrt{\alpha_i}\, R^{1/2} z_j, \qquad z_j \sim \mathcal{N}(0, I)

   m_j^{i+1} = m_j^{i} + C_{MD}\,\left(C_{DD} + \alpha_i R\right)^{-1}\left(d_j - g(m_j^{i})\right)

where :math:`g` is the forward model run over the **whole** period, :math:`C_{MD}` the
covariance between parameters and predictions and :math:`C_{DD}` the covariance of the
predictions. The inflation coefficients must satisfy

.. math::

   \sum_{i=1}^{N_a} \frac{1}{\alpha_i} = 1

so that, for a linear Gaussian problem, the result equals one Kalman update with the original
error :math:`R`. The observations are perturbed again at every iteration and for every member.

Try it without CATHY
--------------------

``examples/DA/plot_6_ESMDA_synthetic.py`` runs the whole workflow on a toy model: alpha
schedule, forward-model callable, result object and plots. The generic loop is
``pyCATHY.DA.esmda.run_esmda``:

.. code-block:: python

   from pyCATHY.DA.esmda import ESMDAConfig, run_esmda
   from pyCATHY.plotters.esmda_plots import plot_esmda_diagnostics

   def forward(param, valid):          # param: (n_param, n_ens)
       prediction = my_model(param)    # (n_obs, n_ens), all observation times stacked
       ok = np.ones(param.shape[1], dtype=bool)   # False for a member that failed
       return prediction, ok

   result = run_esmda(prior, forward, data, data_var,
                      ESMDAConfig(n_iter=4, alpha_schedule="geometric", seed=1))
   plot_esmda_diagnostics(result, names=["Ks", "porosity"])

Run it on a CATHY ensemble
--------------------------

Prepare the observations and the perturbed parameters exactly as for the sequential DA
(see ``plot_3_run_sequentialDA_SMC.py``: ``dict_obs``, ``stacked_data_cov``, the atmospheric
boundary conditions over the **full** period, ``perturbate``), then:

.. code-block:: python

   from pyCATHY.DA.esmda import ESMDAConfig

   result = simuWithDA.run_DA_esmda(
       dict_obs=dict_obs,
       list_parm2update=["Ks"],                 # exact keys or roots of numbered keys; "St. var." is ignored
       dict_parm_pert=var_per_dict_stacked,
       list_assimilated_obs=["tensio"],         # "tensio", "swc" and/or "ERT" (or "all")
       config=ESMDAConfig(n_iter=4, alpha_schedule="geometric", seed=1),
   )
   # posterior parameters: var_per_dict_stacked["Ks"]["update_nb1"]
   # prior is kept under "ini_perturbation"; history: simuWithDA.esmda_result

What it does each iteration: writes the current parameters in every member folder, runs the
simulation from t = 0, reads the pressure head / saturation at **every** observation time,
stacks the predictions in the order of the data and updates the parameters. After the last
update the model is run once more, so ``result.prediction_history[-1]`` holds the posterior
predictions.

Which parameters are estimated
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ES-MDA updates **parameters only**, so ``list_parm2update`` must name at least one. The ``listUpdateParm`` of a
scenario written for the sequential filters is ``"St. var."`` plus the parameters updated alongside the states, and
it does **not** list the initial condition ``ic`` (the state update corrects it): used as is, a scenario such as
``nc_ERT_ic_VGP_pert`` leaves nothing to estimate for the smoother. Add ``"ic"`` (the homogeneous initial
pressure head, written as a parameter at every iteration) to the list:

.. code-block:: python

   list_parm2update = list(scenario["listUpdateParm"])          # ["St. var."]
   list_parm2update.append("ic")                                # ["St. var.", "ic"]: ic is estimated

Parameters that are perturbed but not listed keep their prior sample in every member (model uncertainty, not
estimated). A *layered* initial condition (``pert_control_name = "layers"``) is not supported: the library writes it
from another key than the one the smoother updates, so it raises instead of silently ignoring the updates.

Check the setup first: ``dry_run``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``run_DA_esmda(..., dry_run=True)`` creates the ensemble folders, writes the printing times, stacks the
observations (mask and log10 included), resolves the parameter names and runs the checks, then returns a
summary ``dict`` **without simulating anything**: the parameters, the number of data (raw, kept, masked), for
each survey the stamp and the state time it is compared with, the printing times, the alphas, and the cost
(``n_model_runs``, ``n_member_simulations``, ``n_ert_forwards``). A wrong parameter name, an incomplete
observation entry, a stamp that is not an atmbc time or a covariance of the wrong size shows up in seconds
instead of after the first iteration. Run it first on any new setup.

Supported today
~~~~~~~~~~~~~~~


.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Item
     - Status
   * - ``tensio`` and ``swc`` observations
     - supported
   * - ``ERT`` observations
     - supported, see `ERT observations`_ (needs ``parallel=True`` and the Archie parameters)
   * - EM, ET observations
     - ``NotImplementedError`` (their operators read only the last state of a window)
   * - parameters
     - any parameter ``perturbate`` builds, given by exact key (``"ic"``, ``"SATCELL_VG0"``) or by the root of
       numbered keys (``"Ks"`` -> ``"Ks0"``, ``"Ks1"``, ...), never by substring (``"ic"`` does not select
       ``"a_Archie"``); an entry that matches nothing raises. Log transforms and the scenario bounds
       (``per_bounds``) are respected after every update (``param_bounds_mode = "reflect"`` is honoured)
   * - state update (``"St. var."``)
     - not done: the state is never updated, an initial-condition *parameter* (``ic``) is updated like any other
   * - localisation
     - covariance localisation through ``ESMDAConfig(localisation_matrix=...)`` (parameters x observations)

Time convention
---------------

The smoother compares each observation with the state of one run over the whole period, so it
must know *which time* an observation stamped ``t`` refers to. ``run_DA_sequential`` has used two
conventions:

.. list-table::
   :header-rows: 1
   :widths: 24 38 38

   * - ``obs_time_convention``
     - Stamp ``t`` is compared with the state at
     - Matches ``run_DA_sequential``
   * - ``"state_time"`` (default)
     - ``t`` itself
     - since commit ``d37036f9`` (2026-10-19): the analysis is made when the window **ends** at ``t``
   * - ``"window_start"``
     - the atmbc time that **follows** ``t``
     - up to commit ``083f0ac``: the analysis was made after the window that **starts** at ``t``;
       also what notebook helpers such as ``obs_window_start`` assume (stamp = atmbc time just before
       the acquisition)

Use the one your observations were stamped for. Stamped consistently, both pair a survey with the
state at its acquisition time; stamped with the wrong one, every survey is compared with a state one
model step off, and the misfit cannot be trusted. ``pyCATHY.DA.time_alignment`` still describes the
older sequential pairing (``t`` compared with the next time): it does not apply to the smoother.
With ``"state_time"`` the stamps must be CATHY printing times; the adapter prints every state time.

ERT observations
----------------

Each iteration maps the saturation of every member at every survey time through the Archie law and
the pygimli forward (the sequential ERT operator, called with the row of the member ``sw`` output
that belongs to the survey time instead of the last row). The data then follow the same steps as in
the sequential analysis:

- unusable quadrupoles are dropped from data **and** predictions
  (``DAConfig.mask_bad_ert``, ``ERT_checks.quadrupole_quality_mask``);
- data and predictions are put in log10 and R is propagated to log10
  (``DAConfig.use_log_transformed_obs`` / ``log_transform_cov``), so the misfit is in log10 units;
- predictions that are not finite (a member whose forward failed for a quadrupole) leave that
  quadrupole out of the update.

Requirements: ``parallel=True``; ``set_Archie_parm`` done; the printing of the vtk output
(``VTKF=2`` in ``parm``, as in the sequential notebooks); the forward must return the same number
of quadrupoles as the data file (otherwise ``... the forward gives N data, the observations hold M``).
The cost of one iteration is one simulation of every member over the whole period plus one ERT forward
per member and survey.

Options
-------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - ``ESMDAConfig`` option
     - Meaning
   * - ``n_iter``
     - Number of assimilations (default 4). More iterations = smaller, safer updates but more runs.
   * - ``alpha_schedule``
     - ``"constant"`` (``alpha = n_iter``), ``"geometric"`` (large first, then smaller) or ``"user"``
       with ``user_alphas``. The constraint ``sum(1/alpha) = 1`` is checked.
   * - ``geometric_ratio``
     - Ratio between successive alphas of the geometric schedule (default 0.7).
   * - ``seed``
     - Seed of the observation perturbations, for reproducible runs.
   * - ``inflation_params``
     - Multiplicative inflation of the parameter anomalies after each update (default 1: none).
   * - ``localisation_matrix``
     - ``(n_param, n_obs)`` taper on the parameter-observation covariance.
   * - ``obs_time_convention``
     - Which state a stamp ``t`` is compared with: ``"state_time"`` (default) or ``"window_start"``,
       see `Time convention`_.
   * - ``update_initial_conditions``, ``cleanup_between_iterations``
     - Not supported yet: ``True`` raises ``NotImplementedError``.

Reading the result
------------------

``ESMDAResult`` (returned and kept as ``simuWithDA.esmda_result``):

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Field
     - Content
   * - ``param``, ``valid``
     - Final ensemble ``(n_param, n_ens)`` and the members still active
   * - ``param_history``
     - ``n_iter + 1`` snapshots: prior, then after every update
   * - ``misfit_history``
     - Reduced chi-square of every forward run, about **1** when the data are fitted
       within their error; it should fall from the first to the last run
   * - ``prediction_history``, ``valid_history``, ``alphas``
     - Predictions of every run (NaN for unusable members), active members at every update,
       coefficients used

Members that fail (rejected by the mass-balance check or unreadable output) are dropped for
the rest of the run with a warning. Observations with a non-finite prediction are left out of
that update. The run stops if fewer than two members remain.

What the run tells you
----------------------

Before the first model run (and also with ``dry_run=True``) the strategy is printed and logged (INFO) as one block:

- the estimated parameters with their prior mean and standard deviation (physical units);
- the iterations: number of updates, the alphas with the schedule and the ``sum(1/alpha) = 1`` check, inflation,
  localisation, seed;
- the members and the water-balance threshold above which a member is dropped, the bounds mode;
- the observations: sensors, number of assimilation windows, number of data (masked out of raw), log10;
- the time convention and the **simulated period** (first and last atmbc time in seconds, hours and days, the number
  of atmbc times and atmbc windows);
- the cost: model runs, member simulations, ERT forwards;
- a table of the **assimilation windows**: for each one the stamp, the state time it is compared with, where the window
  starts (the first atmbc time, then the previous state time), its length and the number of atmbc times in it.

The model is **not** restarted at these windows, as it is in the sequential DA: the whole period is run from the start
at every iteration, and the windows are the instants at which the state of that one run is compared with the data.
The same information is returned by ``dry_run=True`` as a ``dict`` (``period``, ``windows``, ``n_atmbc_times``,
``alphas``, ``n_data``, ...).

During the run, one INFO line per iteration gives the number of valid members, the misfit (reduced chi-square) of the
forward run and the alpha used for the update, then the mean and standard deviation of the first updated parameters;
the last line gives the misfit of the final run with the posterior.

Console log
-----------

As ``run_DA_sequential``, ``run_DA_esmda`` copies everything printed during the run (``print``, ``logging``
messages including the INFO progress lines and the warnings about dropped members or excluded quadrupoles, and the
rich console of the CATHY object) to ``<workdir>/<project>/console_log.txt``, live, appending one
``===== run started ... =====`` block per run. The output of the compiled Fortran processor itself and rich output
inside Jupyter are not captured.

That file only exists once the run starts: an error while the script prepares the model (setup, survey filtering,
the pre-run checks) is not in it. To log the whole script, start the same tee at its top (the ERT notebook does it
with ``LOG_CONSOLE = True``, writing next to the project folder, which is deleted and recreated):

.. code-block:: python

   import atexit
   from pyCATHY.DA.runtime import _tee_run_output

   tee = _tee_run_output("../simulations/console_log_my_script.txt")
   tee.__enter__()
   atexit.register(tee.__exit__, None, None, None)   # closes the log, also after an uncaught traceback

What to check on a first run
----------------------------

The adapter assumes, and a first real run should confirm, that:

- the observation times of ``dict_obs`` are on the same time axis as the atmospheric
  boundary conditions and are CATHY **printing times** (the adapter writes the printing times
  into every member's ``parm`` file);
- the member folders hold the full period after ``prepare_DA`` and CATHY overwrites its
  ``psi`` / ``sw`` outputs on every run;
- the mass-balance check of a member uses the end of the window written in ``parm``;
- for ``swc`` observations the 3d mesh is needed: if it is not loaded yet the adapter does the same
  dry run of CATHY (``IPRT1=3``) as the sequential DA before the first run (not done for ``tensio``);
- ``misfit_history`` decreases. If it does not, start from a smaller ``n_iter``, a larger
  prior spread, or check the data error.

Troubleshooting
---------------

.. list-table::
   :header-rows: 1
   :widths: 34 36 30

   * - Message
     - Cause
     - Fix
   * - ``sum(1/alpha) must be 1, got ...``
     - ``user_alphas`` do not satisfy the constraint.
     - Use ``"constant"`` / ``"geometric"``, or correct the list (e.g. ``[2, 4, 4]``).
   * - ``No output at t=...: printed times are [...]``
     - An observation time is not a CATHY printing time.
     - Make the observation times printing times (same time axis as the atmbc).
   * - ``ES-MDA on CATHY supports only ('tensio', 'swc') observations, got 'ERT1'``
     - An unsupported sensor is selected.
     - Restrict ``list_assimilated_obs`` to ``tensio`` / ``swc``.
   * - ``Observation 'swc1' at t=... needs 'data' and 'mesh_nodes'``
     - The observation entry is incomplete (it was not prepared for the DA).
     - Build ``dict_obs`` as in the observation preparation examples (each sensor needs ``data`` and ``mesh_nodes``).
   * - ``member N: cannot read the output (...)`` (warning)
     - The output of a member is missing or truncated (CATHY stopped early); the member is dropped.
     - Look at the CATHY log of that member. Other errors (``KeyError``, ``AttributeError``) are not
       treated as a member failure and stop the run.
   * - ``No parameter to update for list_parm2update=['St. var.']: ES-MDA updates parameters``
     - The list only names the state (a scenario made for the sequential filters).
     - Add the parameters to estimate, ``"ic"`` for the initial condition (see `Which parameters are estimated`_).
   * - ``'ic': a layered initial condition (pert_control_name='layers') is not supported by ES-MDA yet``
     - The layered ic writer reads ``ini_ic_withLayers``, not the key the smoother rewrites.
     - Use a global ic (``control = global`` in the scenario table).
   * - ``list_parm2update entries [...] are not perturbed parameters``
     - An entry is neither a key of ``dict_parm_pert`` nor the root of numbered keys (no substring match).
     - Use the exact key (``"ic"``) or the root (``"SATCELL_VG"`` for ``SATCELL_VG0..``).
   * - ``survey stamped ... s: the state at ... s is not an atmbc time (nearest ... s)``
     - The times of ``sim.atmbc`` are not those the stamps were built from. Typical cause:
       ``sim.read_inputs('atmbc')`` **replaces** ``sim.atmbc`` by the vector read back from the file, which is float32
       (9240.001 s comes back as 9240.0009765625 s, ~0.06 s off at 1e6 s); ``prepare_DA`` copies ``sim.atmbc['time']``.
     - Keep the exact vector and restore it after the call (``exact = sim.atmbc; sim.read_inputs('atmbc'); sim.atmbc = exact``).
       ``check_state_times`` compares with a relative tolerance, so float32 rounding alone is not reported.
   * - ``obs_time_convention='window_start': t=... is not an atmbc time``
     - Stamps meant as window starts are not atmbc times (or the last one).
     - Stamp with atmbc times, or use ``"state_time"``.
   * - ``ERT ... the forward gives N data, the observations hold M``
     - The ERT forward and the data file use different quadrupoles.
     - Same sequence file for the data and the forward; check the quadrupole filtering of the data files.
   * - ``No parameter to update among [...]``
     - ``list_parm2update`` matches no key of ``dict_parm_pert`` (``"St. var."`` is ignored).
     - Name at least one perturbed parameter.
   * - ``forward returned predictions of shape (a, b), expected (c, d)``
     - The stacked predictions do not match the stacked data (different sensors or nodes per time).
     - Check that every selected sensor is present at every assimilation time.
   * - ``Only N valid member(s) left at iteration ...``
     - Too many members failed.
     - Check the rejection ledger (``<project>_rejections.json``), the prior ranges and the time step.
   * - ``No usable observation at iteration ...``
     - Every prediction is non-finite.
     - The runs probably failed: read the CATHY logs of a member.
   * - ``update_initial_conditions=True is not supported yet``
     - Option not implemented.
     - Leave it ``False``.

References
----------

Emerick, A. A., and Reynolds, A. C. (2013). Ensemble smoother with multiple data assimilation.
*Computers & Geosciences*, 55, 3-15.
