import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent))


@pytest.fixture
def rng() -> np.random.Generator:
    return np.random.default_rng(0)
