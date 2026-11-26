import os
import stat

from pyCATHY.cathy_tools import _rmtree


def test_rmtree_removes_read_only_files(tmp_path):
    target = tmp_path / "tmp_src" / ".git" / "objects" / "pack"
    target.mkdir(parents=True)
    pack = target / "pack-x.idx"
    pack.write_text("data")
    os.chmod(pack, stat.S_IREAD)  # read-only, as git writes packs (blocks rmtree on Windows)
    _rmtree(tmp_path / "tmp_src")
    assert not (tmp_path / "tmp_src").exists()
