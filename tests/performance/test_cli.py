import argparse
import os
from unittest import mock

import pyKVFinder

DATADIR = os.path.join(os.path.dirname(pyKVFinder.__file__), "data")

def _standard_args():
    return argparse.Namespace(
        input=os.path.join(DATADIR, "tests", "1FMO.pdb"),
        verbose=False,
        base_name=None,
        output_directory=os.path.join(DATADIR, "tests", "output"),
        model=None,
        nthreads=os.cpu_count() - 1,
        dictionary=os.path.join(DATADIR, "vdw.dat"),
        step=0.6,
        probe_in=1.4,
        probe_out=4.0,
        volume_cutoff=5.0,
        removal_distance=2.4,
        surface="SES",
        ignore_backbone=False,
        depth=False,
        plot_frequencies=False,
        hydropathy=False,
        box=None,
        ligand=None,
        ligand_cutoff=5.0,
    )

def test_cli_standard_mode_perf(benchmark):
    with mock.patch(
        "argparse.ArgumentParser.parse_args",
        return_value=_standard_args(),
    ):
        result = benchmark(pyKVFinder.main.cli)

    assert result == 0

def test_cli_depth_mode_perf(benchmark):
    args = _standard_args()
    args.depth = True

    with mock.patch(
        "argparse.ArgumentParser.parse_args",
        return_value=args,
    ):
        benchmark(pyKVFinder.main.cli)

def test_cli_ligand_mode_perf(benchmark):
    args = _standard_args()
    args.ligand = os.path.join(DATADIR, "tests", "ADN.pdb")

    with mock.patch(
        "argparse.ArgumentParser.parse_args",
        return_value=args,
    ):
        benchmark(pyKVFinder.main.cli)

def test_cli_hydropathy_mode_perf(benchmark):
    args = _standard_args()
    args.hydropathy = "EisenbergWeiss"

    with mock.patch(
        "argparse.ArgumentParser.parse_args",
        return_value=args,
    ):
        benchmark(pyKVFinder.main.cli)
