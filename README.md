# katcali: single dish calibration pipeline for MeerKAT


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
cp -r katcali/katcal_container ~/.local/share/jupyter/kernels/
```
cd ~/.local/share/jupyter/kernels/katcal_container
```
emacs kernel.json
change $username to your own username

