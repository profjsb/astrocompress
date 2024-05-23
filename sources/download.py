"""
Download KOA data
"""

import argparse
from pathlib import Path

from pykoa.koa import Koa
from astropy.io import ascii

koa_source_dir = Path("KOA")
data_dir = Path("../data/integer/")
#instruments = ["esi", "hires", "lris", "nirc1", "nirc2"]
instruments = ["lris"]
#modes = ["imaging", "spectroscopy"]
modes = ["imaging"]
testing = True
#im_types = ["sci", "calib"]
im_types = ["sci"]

def download_data(testing=False):

    for instrument in instruments:
        for mode in modes:
            for im_type in im_types:
                out_dir = Path(f"../data/integer/{mode}/{im_type}/")
                out_dir.mkdir(parents=True, exist_ok=True)
                fname = (
                    koa_source_dir
                    / f"raw_{instrument}_{mode}{'_calib' if im_type == 'calib' else ''}.tbl"
                )
                if fname.exists():
                    data = ascii.read(
                        fname, format="ipac", guess=False, fast_reader=False
                    )
                    nrows = len(data)
                    print(fname, nrows)
                    Koa.download(
                        fname._str,
                        "ipac",
                        out_dir._str,
                        start_row=0,
                        end_row=0 if testing else nrows,
                        calibdir=0,
                    )


if __name__ == "__main__":

    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Download KOA data")
    parser.add_argument(
        "--full", action="store_true", default=False, help="Download full data"
    )
    args = parser.parse_args()
    full = args.full
    download_data(testing=not full)
    print("Done")
