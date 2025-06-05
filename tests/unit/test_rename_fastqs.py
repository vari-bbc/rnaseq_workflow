import subprocess as sp
import os
import pytest
from pathlib import Path

@pytest.mark.skipif(not os.path.isdir("results"), reason = "No results dir")
def test_rename_fastqs():
    reference = Path("tests/unit/rename_fastqs/expected/SRR1039508_0_R1.fastq.gz")
    output = Path("results/rename_fastqs/SRR1039508_0_R1.fastq.gz")
    try:
        sp.check_output(["cmp", reference, output])
    except sp.CalledProcessError:
        assert False, "It doesn't match expected output."