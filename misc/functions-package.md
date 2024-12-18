## Description of main functions in the Python package

To install: `pip install isweep`.

To get more details: 'help(your_function)'.

All the code is in `src/isweep/`.

### Main functions:
---
In `coalescent.py`
- `walk_variant_backward` : Wright-Fisher process backward in time (with selection).
- `walk_variant_forward` : Wright-Fisher process forward in time (with selection).
- `simulate_ibd_isweep` : generate IBD segments around a locus for a selective sweep.
    - Set `p0 = 0` for neutral setting (no sweep).
    - In no sweep, and constant Ne setting, `simulate_ibd_constant` is a really fast option.
- There are `_tv` versions of the `walk*` and `simulate_ibd*` functions for time-varying selection.
In `inference.py`
- `read_ibd_file` : input file is the output style of the Browning Lab `hap-ibd.jar`
- `chi2_isweep` : use this as function to optimize with `scipy.optimize`.
    - See example usage in `vignettes/` folder.
- `bootstrap_*` : different bootstrap interval functions (standard, hall, efron).
In `utilities.py`
- `read_Ne` : load in IBDNe input style file as a dictionary
    - Other functions `*_Ne*` create or modify Ne dictionaries.
- `bin_ibd_segments` : split into bins by IBD segment length.
    - Output is used in `chi2_isweep` for estimation. See example usage in `vignettes/` folder.
In `outgroups.py`
- `make_ibd_graph` : form IBD graph defined in Temple, Waples, Browning (AJHG, 2024). Input is output of `read_ibd_file`.


### Unlikely to custom use
---

- `slow.py` should never be used. It is for comparing algorithm speeds.
- `outgroups.py` is incorporated into the modeling pipeline.
- `favoredalleles.py` is incorporated into the modeling pipeline.