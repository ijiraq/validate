import asyncio
from mocas.cadc_archive_tap_client import CADCArchiveTapClient
import vos
from astropy import units
import logging
import argparse
from astropy.table import Table

from .app import ObservationValidator
from .artifact import dbimages_artifact_from_mocas_row, DBImagesArtifact, ComparisonImageFinder
from .downloaders import AsyncVOSpaceFITSCutoutDownloader, VOSpaceApcorDownloader, VOSpaceFITSCutoutDownloader


def main() -> None:
    """
    The validate app
    :return:
    """

    parser = argparse.ArgumentParser(description='Display an image from the MOCAS database.')
    parser.add_argument('mocas_table', help='The MOCAS table to display.')
    parser.add_argument('--dbimages', help='The DBImages vos root. default: %(default)s',
                        default='vos:Kozai/MOP/dbimages/')
    parser.add_argument('--log-level', default='ERROR', help='The logging level. default: %(default)s',
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'])
    args = parser.parse_args()
    logging.basicConfig(level=args.log_level)

    mocas_table = Table.read(args.mocas_table, format='ascii.csv', delimiter='\t')
    available_observation_ids = vos.Client().listdir(args.dbimages)

    artifacts = get_list_of_artifacts_to_examine(mocas_table, available_observation_ids)

    asyncio.run(async_download(artifacts, args.dbimages))

    apcor_downloader = VOSpaceApcorDownloader(dbimages=args.dbimages)
    fits_cutout_downloader = VOSpaceFITSCutoutDownloader(dbimages=args.dbimages,
                                                         minimum_cutout_radius=10 * units.arcsec)
    cadc_archive_tap_client = CADCArchiveTapClient()
    comparison_finder = ComparisonImageFinder(tap_client=cadc_archive_tap_client,
                                              fits_cutout_downloader=fits_cutout_downloader,
                                              available_observation_ids=available_observation_ids)
    obs_records = ObservationValidator(artifacts=artifacts,
                                       apcor_downloader=apcor_downloader,
                                       comparison_finder=comparison_finder).create_obs_records()
    for artifact_id in obs_records:
        print(obs_records[artifact_id].to_string())


def get_list_of_artifacts_to_examine(mocas_table, available_observation_ids):
    artifacts = []
    for mocas_row in mocas_table:
        artifact = dbimages_artifact_from_mocas_row(mocas_row)
        if artifact.observation_id not in available_observation_ids:
            continue
        artifacts.append(artifact)
    return artifacts


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


