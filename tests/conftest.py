import os
import shutil

import pytest

from diffquant import config as C

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


@pytest.fixture()
def example_cfg(tmp_path):
    cfg = C.load(os.path.join(ROOT, "examples", "config.yaml"))
    cfg["output"] = str(tmp_path / "out")
    cfg["jobs"] = 1
    return cfg
