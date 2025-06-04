import subprocess as sp
from pathlib import Path

def test_salmon():
    reference = Path("tests/unit/salmon/expected/SRR1039508/quant.sf")
    output = Path("results/salmon/SRR1039508/quant.sf")
    with open(reference, "r") as f1, open(output, "r") as f2:
        assert sum(1 for _ in f1) == sum(1 for _ in f2)