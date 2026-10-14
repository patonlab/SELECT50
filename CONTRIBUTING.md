# Contributing SELECT50 results

The [data guide](results/README.md) describes the benchmark files, methods and conventions. Use the [reference XYZ files](results/reference/starting_structures/) as inputs and retain their case IDs and fixed `major`/`minor` roles.

1. Add optimized XYZ files under `optimized-geometries/<method-version>/<structure_id>/<structure_id>.xyz`, retaining the reference atom order.
2. Place raw electronic and Gibbs energy CSVs in `results/energetics/`. Identify each structure and state energy units explicitly.
3. Describe the model/checkpoint, software versions, charge and spin treatment, optimization settings and convergence criteria. For Gibbs energies, state the temperature and thermochemistry settings. Include hardware and timing definitions when reporting calculation times.
4. Record calculation failures and exclusions explicitly. Keep failed predictions blank in the pair-energy tables.
5. Register new methods in the validator, update the corresponding summary tables and documentation, and run:

```sh
python3 scripts/validate_data.py
python3 -m unittest discover -s tests
```

NumPy is required for the RMSD tests. Preserve the signed minor-minus-major convention and cite the calculation or experimental sources associated with each dataset.
