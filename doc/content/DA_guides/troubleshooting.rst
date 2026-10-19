Troubleshooting
===============

Messages of the scenario checks, the run and the diagnostics, with their cause and the fix.
The first column is the text to search for in the log.

Scenario checks (before the run)
--------------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 36 34

   * - Message
     - Cause
     - Fix
   * - ``unknown parameter 'Zroot' (did you mean 'ZROOT'?)``
     - The name is not one ``perturbate()`` builds; it would be ignored silently.
     - Correct the name in ``per_name`` (the names are case-sensitive).
   * - ``'per_sigma' has 1 entries but per_name has 2``
     - A list is not aligned with ``per_name``.
     - Every list needs one entry per parameter, in the same order.
   * - ``listUpdateParm names that were not built: 'Ks0' ... did you mean 'Ks'?``
     - The update names a key that does not exist: a scalar ``Ks`` is built as ``Ks``,
       per-entry parameters as ``Ks0``, ``Ks1``...
     - Use the keys of the perturbed parameters (``simu.describe_perturbation()`` lists them).
   * - ``pert_control_name 'zone' needs one entry per zone``
     - A control was set but the values are scalars.
     - Give lists with one entry per layer, zone or class (see the guide).
   * - ``per_nom is a list but pert_control_name is not set``
     - A list without a control fails inside numpy for soil parameters and ``ic``.
     - Set ``pert_control_name`` to ``"layers"``, ``"zone"`` or ``"root_map"``.
   * - ``pert_control_name 'root_map' is not supported for this parameter``
     - ``porosity`` and ``ic`` do not support ``root_map``; ``WTPOSITION`` supports no control.
     - Use another control, or another parameter.
   * - ``the Feddes table has 3 rows ... need a list with 3 entries``
     - ``ZROOT`` / ``PCREF`` with several vegetation classes need one entry per row.
     - Give a list as long as the Feddes table (scalar only for one row).
   * - ``must be scalars (one nominal applied to every vegetation class...)``
     - ``PCWLT``, ``OMGC``, ``PZ`` accept one scalar only.
     - Give a scalar.
   * - ``VGP: its perturbation is commented out``, ``alpha_VG: ... Unknown parameter``
     - The parameter cannot be perturbed or written (see the catalogue).
     - Remove it; perturb ``SATCELL_VG`` instead of ``alpha_VG``.
   * - ``zone map classes are [1, 3]: they must be numbered 1..2``
     - Gap in the numbering: entry ``k`` is written to class ``k + 1``.
     - Renumber the map from 1 without gaps.
   * - ``per_bounds ...: unknown bound key(s)``, ``min (2) must be lower than max (1)``
     - Malformed bounds (used to be ignored silently).
     - Use ``{ min = a, max = b }`` or ``[a, b]``, with ``a < b``.
   * - ``perturbation layout: Ks0, Ks1 all write PERMX/PERMY/PERMZ of row (1, 0): the last one written wins``
     - Entries overlap (a vegetation class spans several zones).
     - Rethink the entries: the earlier perturbation is overwritten and has no effect.

During the run
--------------

.. list-table::
   :header-rows: 1
   :widths: 30 36 34

   * - Message
     - Cause
     - Fix
   * - ``observation times: N of M observation time(s) will never be assimilated``
     - An observation time is not exactly a window start (a rounded float, or the final time).
       The k-th analysis uses the k-th observation, so **every later observation is paired
       with the wrong analysis**.
     - Make the observation times exactly equal to the atmospheric-forcing times, and drop
       the observation at the final time.
   * - ``observation(s) are assimilated after the window that STARTS at their time``
     - An observation stamped ``t`` is compared with the model state at ``t`` plus one window
       (the last time of the member output).
     - Check that the data were generated with the same convention.
   * - ``ERT ...: N/M quadrupoles excluded from the analysis {...}``
     - Quadrupoles with an unusable geometric factor, non-positive rhoa or an invalid flag.
       They stay in the data and in the reports; only the filter skips them.
     - ``<project>_ERT_checks.json`` lists the indices and reasons.
   * - ``number of rejected ensemble is too high: ... Rejected members by stage``
     - More members failed than ``threshold_rejected`` allows.
     - ``<project>_rejections.json`` has the cycle, stage (``after_run`` or ``after_analysis``)
       and reason of every rejection.
   * - ``... analysed value(s) clipped to [a, b]``
     - The analysis moved a parameter outside its bounds.
     - Expected from time to time; many members at a bound mean the data push the parameter
       out of its range.
   * - ``positive pressure heads in the analysed ensemble``
     - The update itself produced positive pressure heads (the forecast ones are reported
       separately).
     - Check the observation error and the gain; positive heads are reported, not clipped.
   * - ``mbeconv ... has no time step``
     - The member run stopped before writing its first time step.
     - Read the CATHY log of that member.
   * - ``localisation matrix has N columns but M observations were available``
     - The matrix matches neither the observations nor the ones kept by the ERT mask.
     - Build it with one column per observation of the survey.
   * - ``Cannot localise with apparent resistivity data.``
     - Domain localisation was requested with ERT.
     - Use a covariance localisation matrix.
   * - ``Wrong stacked_data_cov shape``
     - The stored covariance does not match the number of data.
     - Recompute it with ``make_data_cov``.
