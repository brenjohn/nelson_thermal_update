# Nelson & Langer thermal update

This project uses the Nelson & Langer (1999) heating and cooling model to advance gas temperatures with fixed chemical abundances, density, and radiation field. It also computes how the updated temperature depends on each input abundance.

## Run the notebook

Install Python and GCC. From this directory, run the following commands and open [nelson_thermal_update.ipynb](nelson_thermal_update.ipynb) in JupyterLab.

```sh
python -m pip install jupyterlab numpy scipy matplotlib
python -m jupyterlab
```

The Python modules build their shared libraries from the generated C sources in `c/` on import if the libraries are missing.

## Use your own tracer data

The temperature and derivative examples run without tracer data. Skip the tracer section or set the data path, `rows_per_tracer`, `rows_per_segment`, and `UVREF1` for your tracer dataset. Update the column mapping if your file uses a different layout.
