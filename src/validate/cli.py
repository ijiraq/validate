import asyncio
import glob
import re

from mocas.cadc_archive_tap_client import CADCArchiveTapClient
import multiprocessing
from mp_ephem import EphemerisReader
import vos
from astropy import units
from astropy.io import fits
import logging
import argparse
from mocas.ssois import SSOIS
from .app import ObservationValidator
from .artifact import dbimages_artifact_from_ssois_row, DBImagesArtifact, ComparisonImageFinder, key_from_mjd
from .downloaders import AsyncVOSpaceFITSCutoutDownloader, VOSpaceApcorDownloader, VOSpaceFITSCutoutDownloader
import inspect
import os.path


def main() -> None:
    """
    The validate app
    :return:
    """

    parser = argparse.ArgumentParser(description='Display an image from the MOCAS database.')
    parser.add_argument('ast_file', help='Starting astromerty file.')
    parser.add_argument('start_date', help='Start date of the time interval to search.')
    parser.add_argument('end_date', help='End date of the time interval to search.')
    parser.add_argument('--dbimages', help='The DBImages vos root. default: %(default)s',
                        default='vos:Kozai/MOP/dbimages/')
    parser.add_argument('--skip-previous', help='Skip previously measured images.',
                        action='store_true', default=False)
    parser.add_argument('--log-level', default='ERROR', help='The logging level. default: %(default)s',
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'])


    args = parser.parse_args()
    logging.basicConfig(level=args.log_level)
    current_obs_records = EphemerisReader().read(args.ast_file)
    obs_records = {}
    for obs_record in current_obs_records:
        obs_records[key_from_mjd(obs_record.date.mjd)] = obs_record

    ssois_search = SSOIS(current_obs_records)
    ssois_table = ssois_search(args.start_date, args.end_date)
    available_observation_ids = vos.Client().listdir(args.dbimages)

    artifacts = get_list_of_artifacts_to_examine(ssois_table, available_observation_ids)

    # possible skip previously measured images
    artifacts_to_use = []
    for artifact in artifacts:
        matching_obs_record = lookup_obs_record(obs_records, artifact)
        if args.skip_previous and matching_obs_record is not None:
            continue
        artifact.obs_record = lookup_obs_record(obs_records, artifact)
        artifacts_to_use.append(artifact)
    artifacts = artifacts_to_use

    # start download processes.
    mp_download(artifacts, args.dbimages)

    # create some utility service for the main application.
    apcor_downloader = VOSpaceApcorDownloader(dbimages=args.dbimages)
    fits_cutout_downloader = VOSpaceFITSCutoutDownloader(dbimages=args.dbimages,
                                                         minimum_cutout_radius=10 * units.arcsec)
    cadc_archive_tap_client = CADCArchiveTapClient()
    comparison_finder = ComparisonImageFinder(tap_client=cadc_archive_tap_client,
                                              fits_cutout_downloader=fits_cutout_downloader,
                                              available_observation_ids=available_observation_ids)

    # Launch the application and update the obs_records with new measurements.
    obs_records.update(ObservationValidator(artifacts=artifacts,
                                            apcor_downloader=apcor_downloader,
                                            comparison_finder=comparison_finder,
                                            obs_records=obs_records).create_obs_records())


    # store the result.
    out_filename, ext = os.path.splitext(args.ast_file)
    out_filename = remove_trailing_number(out_filename)
    count = len(glob.glob(out_filename + '*' + ext))
    out_filename += f".{count}{ext}"
    with open(out_filename, 'w') as out_file:
        for obs_record in obs_records.values():
            out_file.write(obs_record.to_string() + '\n')


def remove_trailing_number(filename: str) -> str:
    """
    Remove the trailing number from a filename.
    """
    return re.sub(r'\.[0-9]+$', '', filename)


def lookup_obs_record(obs_records: dict, artifact: DBImagesArtifact):
    return obs_records.get(artifact.key, None)

def get_list_of_artifacts_to_examine(ssois_Table, available_observation_ids):
    artifacts = []
    for mocas_row in ssois_Table:
        artifact = dbimages_artifact_from_ssois_row(mocas_row)
        if artifact.observation_id not in available_observation_ids:
            continue
        artifacts.append(artifact)
    return artifacts


def mp_download(artifacts: list[DBImagesArtifact], dbimages: str):
    """
    Download cutouts from the artifacts located in dbimages on VOSpace.

    :param artifacts: list of artifact objects that we want to download image data for.
    :param dbimages: location in VOSpace where data is stored. (e.g. vos:Kozai/MOP/dbimages/)
    :return:
    """
    fits_cutout_downloader = VOSpaceFITSCutoutDownloader(dbimages=dbimages,
                                                         minimum_cutout_radius=10 * units.arcsec)
    pool = multiprocessing.Pool(10)
    for artifact in artifacts:
        artifact.async_download = pool.apply_async(fits_cutout_downloader,
                                                   args=(artifact.dbimages_uri,
                                                         artifact.predicted_coordinate,
                                                         artifact.cutout_radius),
                                                   callback=artifact.set_hdulist)


async def async_download(artifacts: list[DBImagesArtifact], dbimages: str):
    """
    Asynchronously download cutouts from the artifacts located in dbimages on VOSpace.

    :param artifacts: list of artifact objects that we want to download image data for.
    :param dbimages: location in VOSpace where data is stored. (e.g. vos:Kozai/MOP/dbimages/)
    :return:
    """
    async_fits_cutout_downloader = AsyncVOSpaceFITSCutoutDownloader(dbimages=dbimages,
                                                                    minimum_cutout_radius=10 * units.arcsec)
    hdu_lists = await asyncio.gather(*[async_fits_cutout_downloader(uri=artifact.dbimages_uri,
                                                                    centre_point=artifact.predicted_coordinate,
                                                                    cutout_radius=artifact.cutout_radius)
                                       for artifact in artifacts])
    for artifact, hdulist in zip(artifacts, hdu_lists):
        artifact.hdulist = hdulist


