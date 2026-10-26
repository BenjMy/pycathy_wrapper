.. _changes:

Changelog
=========

Unreleased
----------

Data assimilation outputs (EnKF, PF and ES-MDA now share one layout):

* ``<project>_states.nc`` is always written (``DAConfig.write_states=True``); the ``<project>_states/`` chunks are only
  crash safety (``save_results``). ES-MDA stores the states of the posterior forward run and ``<project>_params.nc``.
* ``<project>_issues.log``: light log of the problems met (rejected members, ensemble collapse, PF degeneracy, ...).
  ``console_log.txt`` is now opt-in (``DAConfig(console_log="full")``).
* ``<project>_run.json``: manifest written at the end of every run, also a failed one (``DA.write_run_summary``): outputs,
  issue and rejection counts, options (``perturb_obs``, ``sakov``, ``inflation``, ``localisation``, ``threshold_rejected``,
  ``psi_bounds``, updated parameters, ``state_only``), ``wall_time_s``, ``states_nc_bytes`` and package versions.
* ``python -m pyCATHY.DA.run_status``: verdict per run (``ok``, ``partial``, ``no states``, ``no results``, ``failed``).
* ES-MDA also writes the ``<project>_states/`` chunks (with ``save_results``); states carry a ``states_kind`` attribute
  (``forecast_before_update`` for EnKF/PF, ``smoothed_posterior`` for ES-MDA) that also survives a chunk-only reload.
* The issue log now warns about options the chosen filter ignores or that degrade the ensemble: ``perturb_obs=False`` (default, the
  posterior variance is underestimated), inflation ignored by ``enkf_dual_analysis*`` / ``pf``; the inflation convention of every
  DA type is recorded in ``<project>_run.json`` (``options.inflation_convention``).
* ``<project>_run.json`` has a ``parameters`` block (realised prior and final ensemble of every updated parameter); a parameter
  ensemble that collapsed (spread contracted by >= 95 %) or piled on a bound (>= 25 % of the members) is logged in the issue log.
* ``load_backup`` prints one line per run instead of one error block per missing file.

Behaviour changes (round S19), read before re-running an old script:

* **``AnalysisConfig.perturb_obs`` now defaults to ``True``** (was ``False``). Every member assimilates the data plus a draw of
  N(0, R), which restores the ``K R K^T`` term of the posterior covariance (without it the posterior variance was 20-50 % of the
  Kalman one in a linear test and the ensemble collapsed over the cycles). Runs that use the default config (``analysis=None``) with
  ``enkf_Evensen2009*`` or the localized filter give different numbers; ``AnalysisConfig(perturb_obs=False)`` reproduces the old ones.
  The ``*_inflation*``, ``*_dual_*``, ``pf`` and ``esmda`` filters are bit-identical. The low-level ``enkf.enkf_analysis`` keeps
  ``perturb_obs=False``. ``<project>_run.json`` records ``perturb_obs`` (effective value) and ``perturb_obs_applied``.
* Lognormal priors: new option ``lognormal_center`` (``"median"`` = historical default, bit-identical; ``"mean"``: the nominal
  value becomes the ensemble MEAN, ``ln x ~ N(ln(nominal) - sd^2/2, sd^2)``). Set it per call (``sampling_dist(..., center="mean")``),
  per parameter entry (``parm["lognormal_center"]``) or for a block of code (``with perturbate.using_lognormal_center("mean"):``).
  The median-centred ensemble mean is higher than the nominal value by ``exp(sd^2/2)`` (+5 % at sd 0.3, +28 % at 0.7, +65 % at 1.0).
  Not yet readable from the scenario file. With clip bounds the law is truncated and the mean is only approximately the nominal value.
* ``read_observations``: a second sensor of another type at the same time no longer REPLACES the first one (the name was bumped
  from the last key only); ``discharge`` is read; a file name for ``swc`` / ``scale`` / ``tensio`` no longer dies on an undefined
  variable; an unknown ``data_type``, a missing ``tA`` or a missing file reader raise an explicit error; ``data_format`` is always
  stored for ERT / EM and ``meta`` is optional.
* ``make_data_cov``: ``list_assimilated_obs='swc'`` (a string, as in the documentation) selected every sensor whose name contains
  one of the LETTERS; ``['all']`` selected nothing. The selection is now the rule of the DA loop (``observations.select_sensors``,
  also used by ``DA._obs_key_select``), the sensors are stacked in the order of ``dict_obs[tA]`` (the order of the data vector; a
  different order at another time used to permute the variances silently), a sensor with n data gets n variances, resipy ERT uses
  ``resist`` (pygimli: ``rhoa``) for a relative error, and EM coils are read from the entry.

Version 0.1.2
-------------

Documentation:

* New instruction to install Fortran and the algebra libraries BLAS and LAPACK on Fedora (see (`#PR19 <https://github.com/BenjMy/pycathy_wrapper/pull/19>`__))
* New instruction to install a virtual machine for windows (see [#PR?](?))

This release contains contributions from:

* Luca Peruzzo (@Peruz)
* Francesca Zanetti
