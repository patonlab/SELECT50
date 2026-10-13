# Contributing SELECT50 results

Download the [100 current starting XYZ files](results/reference/starting_structures/) and check [repository_status.md](repository_status.md) for current coverage, fixed electronic-energy roles and remaining provenance gaps.

1. Use the reference XYZ coordinates as inputs and retain current case IDs and `major`/`minor` roles. Document deviations from the reference protocol.
2. Add optimized structures under `optimized-geometries/<method-version>/<structure_id>/<structure_id>.xyz`. Keep original benchmark data unchanged.
3. Include a method `README.md` stating model/checkpoint, software versions, charge and spin treatment, optimization and convergence settings, hardware if relevant, and run date.
4. Include starting-structure IDs, raw outputs, an energy table with explicit units, and failures or excluded cases. Distinguish electronic energies from free energies and state temperature and thermochemistry settings where applicable.
5. Run `python3 scripts/validate_data.py` and `python3 -m unittest discover -s tests` before submitting a pull request. Register a new method's table label and folder in the validator when adding it to summary tables, and update documented counts.

The [data guide](results/README.md) describes current files. Preserve the recorded aliases and fixed E/G roles. Update the validator’s required columns for intentional schema changes and document any changed failure/exclusion policy. Missing calculations must remain explicit; do not fill them with zero or a substitute geometry.
