pyKVFinder.run_workflow
=======================

.. autofunction:: pyKVFinder.run_workflow(input: Union[str, pathlib.Path], ligand: Optional[Union[str, pathlib.Path]] = None, vdw: Optional[Union[str, pathlib.Path]] = None, box: Optional[Union[str, pathlib.Path]] = None, step: Union[float, int] = 0.6, probe_in: Union[float, int] = 1.4, probe_out: Union[float, int] = 4.0, removal_distance: Union[float, int] = 2.4, volume_cutoff: Union[float, int] = 5.0, ligand_cutoff: Union[float, int] = 5.0, include_depth: bool = False, include_hydropathy: bool = False, hydrophobicity_scale: Union[str, pathlib.Path] = 'EisenbergWeiss', surface: str = 'SES', ignore_backbone: bool = False, model: Optional[int] = None, nthreads: Optional[int] = None, verbose: bool = False) -> pyKVFinder.pyKVFinderResults
