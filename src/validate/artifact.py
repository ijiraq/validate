"""
The observation class holds all the meta information about an observation.
"""
import os
import inspect
import logging
from dataclasses import dataclass, field
from astropy import units
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Row
from astropy.units import Quantity

from mocas.cadc_archive_tap_client import CADCArchiveTapClient

from .downloaders import VOSpaceFITSCutoutDownloader
band_filter_map = {'gri.MP9605': 'w'}

WAITING_IMAGE_PATH = os.path.join(os.path.dirname(inspect.getfile(inspect.currentframe())),
                                  'data',
                                  'waiting.fits')


@dataclass
class ApertureCorrection:
    observation_id: str
    extname: str
    ccd_num: int
    apcor_uri: str = field(init=False)


def key_from_mjd(mjdate):
    return f"{mjdate:.3f}"

@dataclass
class DBImagesArtifact:
    frame: str
    predicted_coordinate: SkyCoord
    ephemeris_uncertainty: (Quantity, Quantity, Quantity)
    obs_date: float
    key: str = None
    # cutout: HDUList = field(default_factory=HDUList)
    observation_id: str = field(init=False)
    dbimages_uri: str = field(init=False)
    cutout_radius: Quantity = field(init=False)
    hdulist: fits.HDUList = None
    apcor_uri: str = None
    comparison: object = None
    obs_record: object = None
    downloaded: bool = False

    def __post_init__(self):
        self.observation_id = self.frame.split('p')[0]
        self.dbimages_uri = f"dbimages:{self.observation_id}/{self.observation_id}p.fits"
        self.cutout_radius = 3*max(self.ephemeris_uncertainty[:2])
        self.key = key_from_mjd(self.obs_date)
        with fits.open(WAITING_IMAGE_PATH, lazy_load_hdus=False) as waiting_image:
            self.hdulist = fits.HDUList([fits.PrimaryHDU(),
                                         fits.ImageHDU(data=waiting_image[1].data,
                                                       header=waiting_image[1].header)])

    def get_apcor_uri(self, extname):
        """
        Get the apcor uri for the given artifact and extension name.
        """
        ccd_num = int(self.hdulist[extname].header['EXTVER'])
        return f"dbimages:{self.observation_id}/{extname}/{self.observation_id}p{ccd_num:02d}.apcor"

    def set_hdulist(self, hdulist):
        self.downloaded = True
        self.hdulist = hdulist


class ComparisonImageFinder:

    def __init__(self, tap_client: CADCArchiveTapClient, fits_cutout_downloader: VOSpaceFITSCutoutDownloader,
                 available_observation_ids: list):
        self.tap_client = tap_client
        self.vospace_fits_cutout_downloader = fits_cutout_downloader
        self.available_observation_ids = available_observation_ids

    def __call__(self, artifact: DBImagesArtifact) -> DBImagesArtifact:
        """
        Find an image that is not artifact that overlaps in sky area.
        """
        overlaping_observations = self.get_list_of_overlapintg_observations_in_cadc_archive(
            artifact.predicted_coordinate.ra.degree,
            artifact.predicted_coordinate.dec.degree,
            artifact.observation_id)
        for row in overlaping_observations:
            if row['observationID'] not in self.available_observation_ids:
                continue
            try:
                artifact = DBImagesArtifact(row['observationID'],
                                            artifact.predicted_coordinate,
                                            artifact.ephemeris_uncertainty,
                                            artifact.obs_date)
                artifact.hdulist = self.vospace_fits_cutout_downloader(artifact.dbimages_uri,
                                                                       artifact.predicted_coordinate,
                                                                       artifact.cutout_radius)
                return artifact
            except Exception as ex:
                logging.warning(f"Failed to download comparison image {row['observationID']}: {ex}")
        return None

    def get_list_of_overlapintg_observations_in_cadc_archive(self, ra: float, dec: float, observation_id: str):
        query = (f"SELECT observationID, (p.time_bounds_upper+p.time_bounds_lower)/2.0 as mjd  "
                 f"FROM caom2.Plane p JOIN caom2.Observation o ON p.obsID=o.obsID "
                 f"WHERE INTERSECTS ( CIRCLE ('ICRS', {ra}, {dec}, 0.02), p.position_bounds ) = 1 "
                 f"AND observationID != '{observation_id}'")
        return self.tap_client.get_table(query)


def dbimages_artifact_from_ssois_row(row: Row):
    return DBImagesArtifact(row['Image'],
                            SkyCoord(row['Object_RA'], row['Object_Dec'], unit='deg'),
                            (row['dra'] * units.arcsec,
                             row['ddec'] * units.arcsec,
                             row['PA'] * units.deg),
                            row['MJD'])
