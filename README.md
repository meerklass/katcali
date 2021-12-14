# katcali: single dish calibration pipeline for MeerKAT

Look at the "examples" folder for information on how to access the data.

## Installation

On Ilifu (https://www.ilifu.ac.za/)

Install to per-user directory:
```
srun --pty bash
```
```
singularity shell /data/exp_soft/containers/katcal.sif
```
```
python setup.py install --user
```

Not on Ilifu

Install to system-wise directory:
```
python setup.py install
```

Install to per-user directory:
```
python setup.py install --user
```

