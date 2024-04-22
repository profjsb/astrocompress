
## Getting the data from sources

### Setup

(First time only) Clone the environment:

```
mamba env create --file astrocompress.yaml
```

Activate the environment:

```
conda activate astrocompress
```

### Download the data

```
cd sources
python download.py
```
This will do a trial run. To pull all the data:

```
python download.py --full
```


### Keck Observatory Archive (KOA)

Link: https://koa.ipac.caltech.edu/

Description: this is ground based imaging a spectroscopy data from a variety of telescopes on the world's largest optical/infrared telescope.

Acknowledgements: This research has made use of the Keck Observatory Archive (KOA), which is operated by the W. M. Keck Observatory and the NASA Exoplanet Science Institute (NExScI), under contract with the National Aeronautics and Space Administration.