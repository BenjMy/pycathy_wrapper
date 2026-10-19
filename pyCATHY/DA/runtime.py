"""Run-time helpers: output tee, working-directory handling, logging set-up.
Moved out of ``cathy_DA`` without change; every name is still importable from
``pyCATHY.DA.cathy_DA``."""
from __future__ import annotations

from typing import Callable, Iterator
import contextlib
import functools
import logging
import os
import sys

logger = logging.getLogger(__name__)


# format of the default handler (same layout as the pyGIMLi log lines)
_LOG_FORMAT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"


_LOG_DATEFMT = "%d/%m/%y - %H:%M:%S"


def _configure_logging(verbose: bool) -> None:
    """Set the verbosity of this module's logger (called when a DA run starts).

    INFO (progress, what is being assimilated) by default, DEBUG (diagnostic
    dumps of ranges, shapes, dictionaries) when ``verbose``; warnings and errors
    are always shown. A stdout handler (timestamp, logger, level) is attached only when the
    application has not configured logging itself, so the messages still appear
    for scripts and notebooks that never call ``logging.basicConfig``.
    """
    level = logging.DEBUG if verbose else logging.INFO
    package_logger = logging.getLogger("pyCATHY")  # parent of every pyCATHY module logger
    package_logger.setLevel(level)
    logger.setLevel(level)
    if not package_logger.handlers and not logger.handlers and not logging.getLogger().handlers:
        handler = logging.StreamHandler(sys.stdout)
        handler.setFormatter(logging.Formatter(_LOG_FORMAT, _LOG_DATEFMT))
        package_logger.addHandler(handler)
        package_logger.propagate = False


class _Tee:
    """File-like object writing to several streams (live, line flushed)."""

    def __init__(self, *streams):
        self.streams = streams

    def write(self, data):
        for s in self.streams:
            try:
                s.write(data)
                s.flush()
            except Exception:  # a broken stream must never stop the run
                pass
        return len(data)

    def flush(self):
        for s in self.streams:
            try:
                s.flush()
            except Exception:
                pass

    def isatty(self):
        return False

    def __getattr__(self, name):  # encoding, fileno, ... of the first stream
        return getattr(self.streams[0], name)


@contextlib.contextmanager
def _tee_run_output(path: str, console=None) -> Iterator[None]:
    """Duplicate everything printed during the block into the log file ``path``.

    Captured live (the file is line buffered):

    - ``print`` and ``rich.print`` (``sys.stdout`` / ``sys.stderr``),
    - ``logging`` messages of handlers that write to stdout/stderr,
    - the rich ``console`` of the CATHY object. When it is not verbose it writes to
      ``os.devnull``: it is temporarily redirected to the log file only, so nothing
      new is displayed on screen.

    Not captured: output written by compiled code directly to the OS file
    descriptors (the Fortran processor) and rich output inside Jupyter.
    """
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    logf = open(path, "a", encoding="utf-8", buffering=1)
    logf.write(
        "\n===== run started {} =====\n".format(
            __import__("datetime").datetime.now().strftime("%d/%m/%y %H:%M:%S")
        )
    )
    out, err = sys.stdout, sys.stderr
    tee_out, tee_err = _Tee(out, logf), _Tee(err, logf)
    sys.stdout, sys.stderr = tee_out, tee_err

    # logging handlers created before the tee still hold the original streams
    rebound = []
    seen = set()
    for lg in (logging.getLogger("pyCATHY"), logging.getLogger(), logger):
        for h in lg.handlers:
            if id(h) in seen or not isinstance(h, logging.StreamHandler):
                continue
            seen.add(id(h))
            if isinstance(h, logging.FileHandler):
                continue
            if h.stream is out:
                h.setStream(tee_out)
                rebound.append((h, out))
            elif h.stream is err:
                h.setStream(tee_err)
                rebound.append((h, err))

    # non-verbose rich console: devnull -> log file only
    old_console_file = None
    if console is not None and getattr(console.file, "name", None) == os.devnull:
        old_console_file = console.file
        console.file = logf

    try:
        yield
    finally:
        if old_console_file is not None:
            console.file = old_console_file
        for h, stream in rebound:
            h.setStream(stream)
        sys.stdout, sys.stderr = out, err
        logf.write("===== run ended =====\n")
        logf.close()


@contextlib.contextmanager
def _working_directory(path: str) -> Iterator[None]:
    """Change the working directory inside the ``with`` block, restore it afterwards (even on error).

    CATHY files are addressed relative to the directory of an ensemble member,
    and some pyCATHY routines read ``os.getcwd()``; this is the only way
    ``cathy_DA`` changes the directory of the process.
    """
    previous = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)


def _restores_cwd(method: Callable) -> Callable:
    """Decorator: the working directory is restored when the method returns."""

    @functools.wraps(method)
    def wrapper(*args, **kwargs):
        with _working_directory(os.getcwd()):
            return method(*args, **kwargs)

    return wrapper
