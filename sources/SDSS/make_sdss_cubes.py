from functools import partial
import warnings
from pathlib import Path

import pandas as pd
import numpy as np
import requests

from astropy.io import fits
from astropy.io.fits.verify import VerifyWarning
from astropy.wcs import FITSFixedWarning
from astropy.wcs import WCS
from astropy.nddata import Cutout2D


import reproject
from reproject.mosaicking import find_optimal_celestial_wcs, reproject_and_coadd
from reproject import reproject_interp
from SciServer import CasJobs

warnings.simplefilter('ignore', category=VerifyWarning)
warnings.simplefilter('ignore', category=FITSFixedWarning)

print(f"reproject version={reproject.__version__}")


def get_most_recent_stripe_rerun(fname="dr7.csv", use_stripe=82):
    """
    Get the most recent rerun for a given stripe from the SDSS

    fname : str
        The file to read in
    use_stripe : int
        The stripe to use, default is 82
    """
    a = pd.read_csv(fname, comment="#")
    str82 = a[a["stripe"]==use_stripe]
    latest = str82.loc[str82.groupby(["strip", "run", "MJD"])["rerun"].idxmax()]
    return latest

def clean_stripe82_fields(fname="Stripe82Fields.csv", quality=[3]):

    """
    Clean the fields in the Stripe82Fields.csv file
    and return the fields that are of the desired quality
    including the URL to download the fits files

    fname : str
        The file to read in, default is Stripe82Fields.csv
    quality : list
        The quality flags to use, default is [3]. 
        Use [2, 3] for more good/ok fields
    """
    df_all = pd.read_csv(fname, index_col=None)
    df_clean = df_all[df_all["quality"].isin(quality)]
    df_clean["url"] = "http://das.sdss.org/raw/" + df_clean["run"].astype(str) + "/" + \
                   df_clean["rerun"].astype(str) + "/corr/" +  \
                    df_clean["camcol"].astype(str) + "/fpC-" + \
                    df_clean["run"].apply(lambda x: '{:06}'.format(x)) + \
                    df_clean["camcol"].apply(lambda x: '-{{band}}{camcol:1}'.format(camcol=x)) + \
                    df_clean["field"].apply(lambda x: '-{:04}.fit.gz'.format(x))
    return df_clean


def get_sdss_image(url, data_dir=Path("../../data/integer/imaging/sci/"), verbose=False):

    """
    Get an SDSS image from a URL

    url : str
        The URL to download the image from
    data_dir : Path
        The directory to save the image in
    verbose : bool
        Whether to print out information about the download
    """
    im_name = data_dir / url.split("/")[-1]

    if not im_name.exists():
        r = requests.get(url, stream=True)
        if r.status_code != 200:
            print(f"Failed to download {url}")
            return None
        with open(im_name, 'wb') as fd:
            for chunk in r.iter_content(chunk_size=128):
                fd.write(chunk)
    else:
        if verbose:
            print(f"{im_name} already exists")

    return(fits.open(im_name))

def make_mosaics_around_field_center(df_clean, loc, bands=["r"], \
    data_dir=Path("../../data/integer/imaging/sci/")):

    """
    Make mosaics around the center of a field. This downloads
    the images, coadds them and saves the mosaics using the field
    centers and flanking fields

    df_clean : pd.DataFrame
        The cleaned dataframe with the fields
    loc : int
        The row in the dataframe to use
    bands : list
        The bands to use, default is ["r"]
    """
    source = df_clean.iloc[loc]
    ra_center = (source.raMax - source.raMin)/2 + source.raMin
    dec_center = (source.decMax - source.decMin)/2 + source.decMin
    run, rerun, camcol, field = source.run, source.rerun,  source.camcol, source.field

    ## get a source near the center of this specific frame
    sql = \
        f"""
        SELECT  top 1 p.ra, p.dec
        FROM PhotoObjAll p
        JOIN dbo.fGetNearbyObjEq({ra_center:0.6f},{dec_center:0.6f},1) n ON n.objID = p.objID
        WHERE n.objID = p.objID 
        and p.run =  {run} and p.rerun = {rerun} and p.camcol = {camcol} 
        and p.field = {field} 
        order by n.distance
        """
    df = CasJobs.executeQuery(sql=sql, context="stripe82")

    # get the images that match this position
    sql1 = f"""
        SELECT
        p.run, p.rerun, p.camcol, p.field
        FROM PhotoObjAll p
        JOIN dbo.fGetNearbyObjEq({df.iloc[0].ra},{df.iloc[0].dec}5,0.03) n ON n.objID = p.objID
        WHERE n.objID = p.objID 
        and p.rerun >= 40 and p.rerun <= 42
        order by n.distance
        """
    df = CasJobs.executeQuery(sql=sql1, context="stripe82")

    # add the download links
    df["url"] = "http://das.sdss.org/raw/" + df["run"].astype(str) + "/" + \
                   df["rerun"].astype(str) + "/corr/" +  \
                    df["camcol"].astype(str) + "/fpC-" + \
                    df["run"].apply(lambda x: '{:06}'.format(x)) + \
                    df["camcol"].apply(lambda x: '-{{band}}{camcol:1}'.format(camcol=x)) + \
                    '-{field:04}.fit.gz'

    center_url=   df.iloc[0]["url"]
    if df.iloc[0].run != run or df.iloc[0].rerun != rerun or df.iloc[0].camcol != camcol \
        or df.iloc[0].field != field:
        print("Center not in the same field, skipping")
        return None

    out_mosaics = []

    for band in bands:
        band_mosaics = []
        for i, row in df.iterrows():
            hdus = []
            used_fields = []
            for field_add in [0, 1, -1]:
                rez = get_sdss_image(row["url"].format(band=band, field=row.field+ field_add), data_dir=data_dir)
                if rez is not None:
                    hdus.append(rez)
                    used_fields.append(row.field+ field_add)

            mosaic_filename = data_dir / f"mosaic_run{row.run}_rerun{row.rerun}_camcol{row.camcol}_f{'.'.join([str(x) for x in used_fields])}_{band}.fits"
            if mosaic_filename.exists():
                print(f"{mosaic_filename} already exists. Skipping.", flush=True)
                band_mosaics.append(mosaic_filename)
                continue

            if len(hdus) == 0:
                print(f"No images found for {band} in {run} {rerun} {camcol} {field}")
                continue

            # get the optimal wcs and make an output array
            wcs_out, shape_out = find_optimal_celestial_wcs(hdus)
            outarray = np.zeros(shape=shape_out, dtype="uint16")

            # make a new header including the WCS from wcs_out and 
            # the data from the first image
            new_header = hdus[0][0].header.copy()

            # remove WCS information
            for key in ["CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2", "CD1_1", "CD1_2", "CD2_1", "CD2_2",
                        "LONGPOLE", "LATPOLE", "RADESYS", "EQUINOX", "CTYPE1", "CTYPE2", "CUNIT1", "CUNIT2",
                        "MJDREF"]:
                try:
                    new_header.remove(key)
                except KeyError:
                    pass

            # add new wcs information from wcs_out
            wcshead = wcs_out.to_header()
            for key in wcshead.keys():
                new_header[key] =wcshead[key]

            # add some metadata
            new_header["BAND"] = band
            new_header["MRUN"] = run
            new_header["MRERUN"] = rerun
            new_header["MCAMCOL"] = camcol
            new_header["CEN_RA"] = ra_center
            new_header["CEN_DEC"] = dec_center
            new_header["MURL"] = center_url
            
            for fn, fields in enumerate(used_fields):
                new_header[f"FIELD{fn}"] = fields

            # reproject and coadd the images
            ri = partial(reproject_interp, order="nearest-neighbor")
            _ = reproject_and_coadd(hdus,wcs_out, shape_out=shape_out,
                                            reproject_function=ri, combine_function="first",
                                            output_array=outarray)
            newhdu = fits.PrimaryHDU(data=outarray, header=new_header)
            newhdu.writeto(mosaic_filename, overwrite=True)
            print(f"Saved {mosaic_filename}")
            band_mosaics.append(mosaic_filename)
        out_mosaics.append(band_mosaics)
    return out_mosaics

def make_data_cube_around_loc(source, mosaic_list, size=(512,512), 
                              data_dir=Path("../../data/integer/imaging/sci/"), 
                              bands=["r"]):
    """
    Make a data cube around a location in the SDSS data

    source : pd.Series
        The source to use as the center of the location
    mosaic_list : list of lists
        The list of lists of mosaic files to use by band
    size : tuple
        The size of the cube to make, default is (512, 512)
    bands : list
        The bands to use, default is ["r"]
    """
    
    run, camcol, field = source.run, source.camcol, source.field

    center_fname =  data_dir / f'fpC-{run:06}-{bands[0]}{camcol:1}-{field:04}.fit.gz'
    if not center_fname.exists():
        print(f"Center file {center_fname} does not exist. Skipping")
        return None

    main_im = fits.open(center_fname)
    im_size = main_im[0].data.shape
    wcs = WCS(main_im[0].header)

    center = (im_size[1]/2, im_size[0]/2)  # switch axes
    cutout = Cutout2D(main_im[0].data, position=center, size=size, wcs=wcs, mode='strict')
    sky_center = wcs.pixel_to_world(*center)

    # make the data cube
    max_mos = max([len(x) for x in mosaic_list])
    data_cube = np.zeros((max_mos, len(bands), size[0], size[1]), dtype="uint16")
    data_hdu = fits.PrimaryHDU(data=data_cube)
    data_hdu.header = main_im[0].header
    data_hdu.header.update(cutout.wcs.to_header())

    for j, band in enumerate(bands):
        data_hdu.header[f"BAND{j}"] = band
        time_sorted_mosaics = sorted([(fits.open(x)[0].header["TAI"], x) for x in mosaic_list[j]])
        for i, (tai, mosaic) in enumerate(time_sorted_mosaics):
            im = fits.open(mosaic)
            data_hdu.header[f"TAI{i}"] = tai
            data_cube[i, j, :, :] = Cutout2D(im[0].data, position=sky_center, size=size, 
                                             wcs=WCS(im[0].header), mode='strict').data.T[::-1,:]

    # data_hdu.scale('int16', bzero=32768)
    cutout_filename = "example_cutout_1.fits"
    data_hdu.writeto(cutout_filename, overwrite=True)
    print(cutout_filename)

def main(df_clean, loc, bands=["r"], data_dir=Path("../../data/integer/imaging/sci/")):

    mosaic_list = make_mosaics_around_field_center(df_clean, loc, bands=bands, \
        data_dir=data_dir)

    source = df_clean.iloc[loc]
    make_data_cube_around_loc(source, mosaic_list, size=(512,512), 
                              data_dir=data_dir, 
                              bands=bands)