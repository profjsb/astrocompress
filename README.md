
## Getting the data from sources

### Setup

You'll need a modern unix-like environment (with `git` and `git-lfs`) where Anaconda can be installed:

```bash
wget https://repo.anaconda.com/archive/Anaconda3-2024.02-1-Linux-x86_64.sh
sh Anaconda3-2024.02-1-Linux-x86_64.sh
```

(First time only) Get this repo and clone the conda environment:

```bash
git clone https://github.com/profjsb/astrocompress.git
cd astrocompress
mkdir -p data/integer/imaging/sci/
conda env create --file astrocompress.yaml
```

Activate the environment:

```
conda activate astrocompress
```

### Download the SDSS Data Cubes

You'll need to get the large CSV containing the list of the SDSS data fields:

```bash
cd sources/SDSS
git lfs pull --include "Stripe82Fields.csv"
```

Now get a sample of cubes. For example:

```bash
python ./make_sdss_cubes.py  -m None --size="(800,800)" --random_state=42 -n 2
```
See the documentation in `make_sdss_cubes.py` for command-line arguments.

### Download the Keck data

```bash
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