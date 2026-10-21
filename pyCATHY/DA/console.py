"""Readable console reports of the parameter update (``[DA-PARAM]`` lines).

Uses ``rich`` (colours, tables, panels) when installed and falls back to plain text with
the same content otherwise. Nothing here changes the run: it only prints.
"""
import sys

try:  # optional dependency
    from rich.console import Console
    from rich.markup import escape
    from rich.panel import Panel
    from rich.table import Table

    _HAS_RICH = True
except ImportError:  # pragma: no cover
    _HAS_RICH = False

TAG = "[DA-PARAM]"


def _console():
    # created per call: follows the current sys.stdout (notebooks, pytest capture, log files)
    return Console(file=sys.stdout, width=140, highlight=False)


def _fmt(x, spec=".4g"):
    try:
        return format(float(x), spec)
    except (TypeError, ValueError):
        return str(x)


def settings(mode, material_lines, n_upd, nens, warn):
    """Banner printed once at the start of a run."""
    rows = [f"out-of-bounds handling after analysis : {mode.upper()}",
            f"parameters updated per cycle          : {n_upd}   (ensemble size {nens})"]
    rows += [f"MATERIAL control : {m}" for m in material_lines]
    if warn:
        rows.append(f"WARNING: {n_upd} parameters >= {nens} members: rank-deficient update, "
                    "expect spurious correlations")
    if not _HAS_RICH:
        print("\n".join([f"{TAG} ===== parameter update settings ====="] + [f"{TAG} {r}" for r in rows]),
              flush=True)
        return
    body = "\n".join(
        ("[bold red]" + escape(r) + "[/]") if r.startswith("WARNING") else
        ("[cyan]" + escape(r) + "[/]") if r.startswith("MATERIAL") else escape(r)
        for r in rows)
    _console().print(Panel(body, title=f"{TAG} parameter update settings", border_style="blue",
                           title_align="left"))


def bounds_alert(cycle, name, n_out, n, lo, hi, mode, raw, final, flag):
    """One line per parameter that left its bounds (red when the spread collapsed)."""
    msg = (f"{TAG} cycle {cycle} {name}: {n_out}/{n} out of [{lo}, {hi}] -> {mode.upper()} | "
           f"raw mean {_fmt(raw[0])} sd {_fmt(raw[1], '.3g')} -> final mean {_fmt(final[0])} "
           f"sd {_fmt(final[1], '.3g')}")
    if flag:
        msg += f"  <<< {flag}"
    if not _HAS_RICH:
        print(msg, flush=True)
        return
    style = "bold white on red" if flag else "yellow"
    _console().print(f"[{style}]{escape(msg)}[/]")


def cycle_table(cycle, rows):
    """Table of one DA cycle. ``rows``: (name, mean_b, mean_a, sd_b, sd_a)."""
    if not rows:
        return
    if not _HAS_RICH:
        for name, mb, ma, sb, sa in rows:
            print(f"{TAG} cycle {cycle} UPDATE {name}: mean {_fmt(mb)} -> {_fmt(ma)} | "
                  f"sd {_fmt(sb, '.3g')} -> {_fmt(sa, '.3g')}", flush=True)
        return
    t = Table(title=f"{TAG} cycle {cycle}: parameter update", title_justify="left",
              header_style="bold", border_style="dim")
    for col in ("parameter", "mean before", "mean after", "shift", "sd before", "sd after", "spread"):
        t.add_column(col, justify="left" if col == "parameter" else "right")
    for name, mb, ma, sb, sa in rows:
        shift = (ma - mb) / abs(mb) * 100 if mb not in (0, None) else float("nan")
        ratio = sa / sb if sb else float("nan")
        spread_style = "red" if ratio < 0.1 else "yellow" if ratio < 0.5 else "green"
        shift_style = "red" if abs(shift) > 50 else "yellow" if abs(shift) > 20 else "green"
        t.add_row(escape(name), _fmt(mb), _fmt(ma), f"[{shift_style}]{shift:+.1f}%[/]",
                  _fmt(sb, ".3g"), _fmt(sa, ".3g"), f"[{spread_style}]x{ratio:.2f}[/]")
    _console().print(t)
