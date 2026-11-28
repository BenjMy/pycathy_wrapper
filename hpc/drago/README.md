# pyCATHY on Drago (CSIC HPC)

Recipe to install pyCATHY and run it with Slurm on Drago. **Written from the public Drago
documentation and tested only locally (Linux + gfortran), not on Drago itself.** Points marked
**[VERIFY]** could not be confirmed from the public docs: check them once, or ask
soporteaic@csic.es.

Files in this folder: `setup_env.sh` (install), `run_case.py` (template case / ensemble member),
`pycathy_single.sbatch` (one run, e.g. the monthly LAI loop), `pycathy_warmcache.sbatch` +
`pycathy_array.sbatch` + `submit_ensemble.sh` (DA ensemble / sweeps).

## Facts used (Drago docs)

- Login: `ssh usuario@drago.csic.es` (SSH key, `ed25519`). Do **not** compute on login nodes.
- Slurm partitions: `generic` (CPU), `gpu`, `special`. Interactive: `srun -p generic -N1 -n1 --pty bash -i`.
- Storage is Lustre: `/lustre/home/<centro>/<usuario>`; quota: `lfs quota -h $HOME`; node-local job dir: `$SLURM_TMPDIR`.
- Modules (Lmod): `module avail | spider | load | purge`. Big transfers with `rsync -aP`.
- Compile in the queue enabled for compilation (or avoid using all cores).

## 1. Install (once)

```bash
ssh usuario@drago.csic.es
srun -p generic -N1 -n1 --mem=8G -t 0-01:00 --pty bash -i     # not on the login node
bash setup_env.sh        # Miniforge + env "pycathy" + conda gfortran + pyCATHY
```

What it does: conda env with `fortran-compiler` (gfortran from conda-forge) and
`pip install pyCATHY` from the refactor branch. **No system gfortran, no BLAS/LAPACK, no
`module load` of a compiler is needed.** Everything lives in `$HOME/pycathy` (override with
`PYCATHY_BASE`).

- **[VERIFY]** `module spider conda` / `module spider miniforge`: if the site provides one, use it and skip the Miniforge download (the script already reuses a `conda` found on the `PATH`).
- **[VERIFY]** internet from the compute node / interactive job (GitHub, conda-forge, PyPI). If blocked, run the installer on the login node (light job) or ask for a proxy/mirror.
- **[VERIFY]** node architecture: `uname -m` (the script picks the matching Miniforge installer).
- A conda env is many small files on Lustre: check `lfs quota -h $HOME`.

## 2. Single run (e.g. monthly LAI loop)

1. Put your model set-up in `run_case.py` (marked block). The LAI → root depth/Feddes loop only
   rewrites input files (`soil`, `root_map`): **no recompilation** as long as `MAXVEG` and the other
   `CATHY.H` values stay fixed (they are, between months).
2. `sbatch pycathy_single.sbatch`

The job runs in `$SLURM_TMPDIR` (CATHY writes many files; local disk is much faster than Lustre) and
`rsync`s everything to `results_<jobid>/` on exit, including on failure or at the time limit.
CATHY is serial: `-n 1`. Adjust `--mem` / `-t` with `seff <jobid>` after the first run.

## 3. DA ensemble / parameter sweep

```bash
bash submit_ensemble.sh 50 10      # 50 members, 10 running at once
```

1. `pycathy_warmcache.sbatch` compiles CATHY once (parallel build) into the build cache.
2. `pycathy_array.sbatch` (job array, `afterok`) runs one member per task; each finds the executable in
   the cache and starts immediately. Cache writes are atomic, so concurrent tasks are safe.

The cache directory is `$PYCATHY_CACHE` (default in the scripts: `$HOME/.cache/pycathy`; use the group
directory to share it with colleagues). It is keyed by sources + `CATHY.H` + compiler + flags: members with
the same dimensions share one binary; a member with a different `CATHY.H` compiles its own (about 5 s with
4+ cores, longer on one core).

## 4. Interactive test

```bash
srun -p generic -N1 -n2 --mem=8G -t 0-01:00 --pty bash -i
conda activate pycathy && export PYCATHY_CACHE=$HOME/.cache/pycathy
python run_case.py $SLURM_TMPDIR/test --build-only
```

## 5. Troubleshooting

| Symptom | Check |
|---|---|
| `No Fortran compiler found` | `conda activate pycathy`; `conda install -c conda-forge fortran-compiler`; or `export PYCATHY_FC=/path/to/gfortran` |
| `CATHY()` fails cloning sources | the project creation downloads `CATHY_src_org` from GitHub: needs internet on the node (see [VERIFY]); create the project in an interactive/login session if not |
| Job killed, memory | `seff <jobid>`; raise `--mem` (bare number = MB) |
| Slow I/O | run in `$SLURM_TMPDIR`, not in Lustre (as in the templates) |
| Always recompiling | `echo $PYCATHY_CACHE` must point to the same directory in every job; `CATHY.H` must be identical |
| Want the real LAPACK | `export PYCATHY_BLAS=system` (needs `conda install libblas liblapack`) |

Sources: [Drago docs](https://docaic.rstools.csic.es/computacion),
[Slurm](https://docaic.rstools.csic.es/computacion/slurm),
[modules](https://docaic.rstools.csic.es/computacion/drago_modules),
[cheat sheet](https://docaic.rstools.csic.es/cheatsheet/cheat_sheet_drago.pdf),
[good practices](https://docaic.rstools.csic.es/general/buenas_practicas).
