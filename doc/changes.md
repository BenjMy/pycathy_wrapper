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
* ``load_backup`` prints one line per run instead of one error block per missing file.

Version 0.1.2
-------------

Documentation:

* New instruction to install Fortran and the algebra libraries BLAS and LAPACK on Fedora (see (`#PR19 <https://github.com/BenjMy/pycathy_wrapper/pull/19>`__))
* New instruction to install a virtual machine for windows (see [#PR?](?))

This release contains contributions from:

* Luca Peruzzo (@Peruz)
* Francesca Zanetti 
