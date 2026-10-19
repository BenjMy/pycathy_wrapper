"""Every library module must compile with warnings as errors (syntax errors, invalid escapes)."""
import pathlib
import warnings

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
SOURCES = sorted(
    p for p in ROOT.rglob("*.py")
    if not {"IA_files", "tests"} & set(p.relative_to(ROOT).parts)
)


@pytest.mark.parametrize("path", SOURCES, ids=lambda p: str(p.relative_to(ROOT)))
def test_module_compiles_cleanly(path):
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        compile(path.read_text(encoding="utf-8"), str(path), "exec")


def test_no_mutable_default_arguments_in_DA_and_ERT():
    """A default list/dict/set is shared between calls (latent bug)."""
    import ast

    offenders = []
    for path in [p for p in SOURCES if p.relative_to(ROOT).parts[0] in ("DA", "ERT")]:
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                defaults = node.args.defaults + [d for d in node.args.kw_defaults if d is not None]
                if any(isinstance(d, (ast.List, ast.Dict, ast.Set)) for d in defaults):
                    offenders.append(f"{path.relative_to(ROOT)}:{node.lineno} {node.name}")
    assert not offenders, offenders
