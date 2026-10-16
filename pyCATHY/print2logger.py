"""AST-based print / self.console.* -> logger conversion (single file). Usage: print2logger.py file [--apply]"""
import ast, re, sys
path = sys.argv[1]; apply = "--apply" in sys.argv
src = open(path, encoding="utf-8").read()
lines = src.splitlines(keepends=True)
offs = [0]
for l in lines: offs.append(offs[-1] + len(l.encode("utf-8")))
b = src.encode("utf-8")
def seg(n): return b[offs[n.lineno-1]+n.col_offset: offs[n.end_lineno-1]+n.end_col_offset].decode()

EMOJI = [":checkered_flag:", ":worried_face:", ":inbox_tray:", ":cooking:", ":athletic_shoe:", ":x:",
         ":white_check_mark:", ":arrows_counterclockwise:", ":sponge:", ":warning:", ":hammer_and_wrench:", ":pensive_face:", ":orange_square:"]
def clean(t):
    t = re.sub(r"\[/?(?:b|red|green|yellow|bold)\]", "", t)
    for e in EMOJI: t = t.replace(e + " ", "").replace(e, "")
    return t
ERR = ("Error", "error", "Failed", "failed", ":x:")
WARN = ("not existing", "non existing", "Warning:", "⚠", "[SKIP]", "WARNING", ":warning:", "not found", "could not", "not recognized", "missing", "unknown", "impossible",
        "⚠️", "cant ", "Falling back", "not installed", "Issue here", "must be", "must start", "Inconsistent", "Cannot", "cannot", "not equal", "No root_map", "no proxy", "not defined", "not a multiple", "more than one outlet", "no file", "Unsupported",
        "not found", "no valid", "NOTE:")
INFO = ("saved", "{:20s} {:>7s}", "{:20s} {:=7.3f}", "created successfully", "Initiate", "Run preprocessor", "Run processor", "gfortran compilation", "Fetch", "fetch", "Compilation successful", "src files")
def level(text, kind):
    if kind == "rule": return "debug"
    if "injected" in text: return "debug"
    if re.search(r"stdout|stderr|output\.decode", text) : return "debug"
    if any(k in text for k in ERR): return "error"
    if any(k in text for k in WARN): return "warning"
    if any(k in text for k in INFO): return "info"
    return "debug"

edits = []
class V(ast.NodeVisitor):
    def visit_Call(self, n):
        self.generic_visit(n)
        f = ast.unparse(n.func)
        if f not in ("print", "self.console.print", "self.console.rule"): return
        kind = "rule" if f.endswith("rule") else "print"
        text = " ".join(seg(a) for a in n.args)
        if "eval(" in text: return
        bad = [k.arg for k in n.keywords if k.arg not in ("style", "end", "sep")]
        assert not bad, (n.lineno, bad)
        if kind == "rule":
            if not n.args or (isinstance(n.args[0], ast.Constant) and not n.args[0].value):
                new = "pass"; lvl = "-"
            else:
                a0 = clean(seg(n.args[0]))
                if "warning messages above" in a0: new = "pass"; lvl = "-"
                else:
                    lvl = "warning" if ":warning:" in seg(n.args[0]) else "debug"; new = f"logger.{lvl}({a0})"
        else:
            if not n.args: new = 'logger.debug("")'; lvl = "debug"
            else:
                lvl = level(text, kind)
                if any(k.arg == "style" and ast.unparse(k.value) == "'yellow'" for k in n.keywords) or "[yellow]" in text: lvl = "warning"
                if len(n.args) == 1: arg = clean(seg(n.args[0]))
                else: arg = '" ".join(str(a) for a in (' + ", ".join(clean(seg(a)) for a in n.args) + ",))"
                new = f"logger.{lvl}({arg})"
        edits.append((n, new, lvl))
V().visit(ast.parse(src))
for n, new, lvl in edits:
    print(f"{n.lineno:5d} {lvl:8s} {new.splitlines()[0][:100]}")
print(len(edits), "calls")
if apply:
    for n, new, lvl in sorted(edits, key=lambda e: (e[0].lineno, e[0].col_offset), reverse=True):
        s = offs[n.lineno-1] + n.col_offset; e = offs[n.end_lineno-1] + n.end_col_offset
        b = b[:s] + new.encode() + b[e:]
    open(path, "w", encoding="utf-8").write(b.decode())
