"""Template: one pyCATHY case (or one DA ensemble member) inside a Slurm job.

    python run_case.py <workdir> [--member N] [--build-only]

Replace the marked block by your own model set-up (DEM, soil, vegetation / LAI loop...).
"""
import argparse
import os

from pyCATHY import cathy_tools


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("workdir")
    ap.add_argument("--member", type=int, default=0)
    ap.add_argument("--build-only", action="store_true", help="only compile (fills the build cache)")
    args = ap.parse_args()

    os.makedirs(args.workdir, exist_ok=True)
    c = cathy_tools.CATHY(
        dirName=args.workdir, prj_name=f"member_{args.member:03d}", verbose=False, clear_src=True
    )

    # ---- your model set-up goes here ------------------------------------------------------
    # e.g. c.update_prepo_inputs(...), c.update_soil(...), c.update_veg_map(...), perturbations
    # of member `args.member`, ... Everything that only touches input files never recompiles.
    # ---------------------------------------------------------------------------------------

    c.run_preprocessor(verbose=False)
    if args.build_only:
        c.recompileSrc()
        assert c.last_build.ok, c.last_build.report()
        print("build OK (cache key %s, from cache: %s)" % (c.last_build.key[:12], c.last_build.from_cache))
        return

    # Monthly LAI loop: do it here, in the same job. Same CATHY.H => the executable comes from
    # the build cache, nothing is recompiled between months.
    c.run_processor(IPRT1=2, TRAFLAG=0, verbose=False)
    assert c.last_build.ok, c.last_build.report()


if __name__ == "__main__":
    main()
