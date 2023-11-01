# KATcali: single dish calibration pipeline for MeerKAT

Look at the "examples" folder for information on how to access the data.

## Installation
### Step I:
#### New user:
```
git clone git@github.com:meerklass/katcali.git
```

#### Old user:
```
cd katcali/
```
```
git pull
```
### Step II:
#### Python3 version:
Congrats you can skip this step

#### Python2 version:
if you want to switch to Python2 version (mostly for 2019 data) then one more line is needed
```
git checkout py2
```
### Step III:
#### On Ilifu (https://www.ilifu.ac.za/)

Install to per-user directory:
```
srun --pty bash
```
```
singularity shell /data/exp_soft/containers/katcal.sif
```
```
python setup.py install --user

Note: you might need to enforce that Python3 is used: python3 setup.py install --user
```

#### Not on Ilifu

Install to system-wise directory:
```
python setup.py install
```

Install to per-user directory:
```
python setup.py install --user
```

