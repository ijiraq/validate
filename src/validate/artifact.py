"""
The observation class holds all the meta information about an observation.
"""
import logging
from dataclasses import dataclass, field
from astropy import units
from astropy.coordinates import SkyCoord
from astropy.io.fits import HDUList
from astropy.table import Row
from astropy.units import Quantity
from mocas.cadc_archive_tap_client import CADCArchiveTapClient

from .downloaders import VOSpaceDownloader, VOSpaceFITSCutoutDownloader
band_filter_map = {'gri.MP9605': 'w'}


@dataclass
class ApertureCorrection:
    observation_id: str
    extname: str
    ccd_num: int
    apcor_uri: str = field(init=False)


@dataclass
class DBImagesArtifact:
    frame: str
    predicted_coordinate: SkyCoord
    ephemeris_uncertainty: (Quantity, Quantity, Quantity)
    # cutout: HDUList = field(default_factory=HDUList)
    observation_id: str = field(init=False)
    dbimages_uri: str = field(init=False)
    cutout_radius: Quantity = field(init=False)
    hdulist: HDUList = None
    apcor_uri: str = None
    comparison: object = None

    def __post_init__(self):
        self.observation_id = self.frame.split('p')[0]
        self.dbimages_uri = f"dbimages:{self.observation_id}/{self.observation_id}p.fits"
        self.cutout_radius = max(self.ephemeris_uncertainty[:2])

    def get_apcor_uri(self, extname):
        """
        Get the apcor uri for the given artifact and extension name.
        """
        ccd_num = int(self.hdulist[extname].header['EXTVER'])
        return f"dbimages:{self.observation_id}/{extname}/{self.observation_id}p{ccd_num:02d}.apcor"


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
        query = (f"SELECT observationID FROM caom2.Plane p JOIN caom2.Observation o ON p.obsID=o.obsID WHERE "
                 f"INTERSECTS "
                 f"( CIRCLE('ICRS', "
                 f"{artifact.predicted_coordinate.ra.degree}, "
                 f"{artifact.predicted_coordinate.dec.degree}, "
                 f"0.02), "
                 f"p.position_bounds ) = 1 "
                 f"AND observationID != '{artifact.observation_id}'")
        result = self.tap_client.get_table(query)
        comparison_id = None
        for row in result:
            if row['observationID'] in self.available_observation_ids:
                comparison_id = row['observationID']
                try:
                    artifact = DBImagesArtifact(comparison_id, artifact.predicted_coordinate, artifact.ephemeris_uncertainty)
                    artifact.hdulist = self.vospace_fits_cutout_downloader(artifact.dbimages_uri,
                                                                           artifact.predicted_coordinate,
                                                                           artifact.cutout_radius)
                    return artifact
                except Exception as ex:
                    logging.warning(f"Failed to download comparison image {comparison_id}: {ex}")
        return None

def dbimages_artifact_from_mocas_row(row: Row):
    return DBImagesArtifact(row['Image'],
                            SkyCoord(row['Object_RA'], row['Object_Dec'], unit='deg'),
                            (row['dra'] * units.arcsec,
                             row['ddec'] * units.arcsec,
                             row['PA'] * units.deg))
