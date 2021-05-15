# katcali: single dish calibration pipeline for MeerKAT

Look at the "examples" folder for information on how to access the data.

## Installation
Install to system-wise directory:
```
python setup.py install
```

Install to per-user directory:
```
python setup.py install --user
```
## Set the Container katcal
```
cp -r katcali/Katcal_container ~/.local/share/jupyter/kernels/
```
```
cd ~/.local/share/jupyter/kernels/Katcal_container
```
```
emacs kernel.json
```
change $username to your own username

