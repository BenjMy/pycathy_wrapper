#!/usr/bin/env python
"""Prove that a change did not alter behaviour: compare the AST of every .py file between
two git revisions (or the working tree). Whitespace and comment changes leave the AST identical.

Usage:  python verify_ast_identical.py [REV_BEFORE] [REV_AFTER]
        (default: HEAD against the working tree; run from the repository root)
Exit status 0 = every file has an identical AST, 1 = at least one differs.
"""
import ast
import subprocess
import sys
from pathlib import Path

before_rev = sys.argv[1] if len(sys.argv) > 1 else "HEAD"
after_rev = sys.argv[2] if len(sys.argv) > 2 else None  # None = working tree


def show(rev, path):
    if rev is None:
        p = Path(path)
        return p.read_text() if p.exists() else None
    r = subprocess.run(["git", "show", f"{rev}:{path}"], capture_output=True, text=True)
    return r.stdout if r.returncode == 0 else None


names = subprocess.run(
    ["git", "diff", "--name-only", before_rev] + ([after_rev] if after_rev else []),
    capture_output=True, text=True, check=True,
).stdout.split()
different, same = [], 0
for path in (n for n in names if n.endswith(".py")):
    a, b = show(before_rev, path), show(after_rev, path)
    if a is None or b is None:
        different.append((path, "added or removed"))
        continue
    if ast.dump(ast.parse(a)) == ast.dump(ast.parse(b)):
        same += 1
    else:
        different.append((path, "AST differs"))
print(f"{same} file(s) with identical AST")
for path, why in different:
    print(f"DIFFERENT: {path} ({why})")
sys.exit(1 if different else 0)
