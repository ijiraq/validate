import asyncio
import contextlib
import logging
from io import BytesIO
import vos
from astropy.io import fits
from requests import ReadTimeout
import time
from urllib.parse import urlparse
from astropy import units


class VOSpaceDownloader:

    def __init__(self, dbimages: str,
                 scheme: str = 'dbimages',
                 retry_limit: int = 10):
        self.root_node = dbimages
        self.scheme = scheme
        self.retry_limit = retry_limit
        self._vos_client = None

    def __call__(self, *args, **kwargs):
        self.check_schema(kwargs['uri'])
        result = self.get_bytes_from_vospace(kwargs['uri'])
        return result

    @staticmethod
    def _cast_result(result):
        return result

    def get_bytes_from_vospace(self, uri: str, cutout_string=None, attempt: int = 0):
        attempt += 1
        logging.debug(f"Getting {uri} using client {self.vos_client} with root node {self.root_node}")
        try:
            with contextlib.closing(self.open(uri, cutout=cutout_string)) as vos_obj:
                result = BytesIO(vos_obj.read())
            return self._cast_result(result)
        except ConnectionError as e:
            logging.error(f"Connection error: {e}")
            if not isinstance(e.args[0], ReadTimeout) or attempt > self.retry_limit:
                raise e
            logging.warning(f"Retrying in 2 seconds (attempt {attempt}) of {self.retry_limit}")
            time.sleep(2)
            return self.get_bytes_from_vospace(uri, cutout_string, attempt=attempt)

    def check_schema(self, uri: str):
        parts = urlparse(uri)
        if parts.scheme != self.scheme:
            raise ValueError(f"URI {uri} scheme {parts.scheme} not supported: expected scheme {self.scheme}.")

    @property
    def vos_client(self) -> vos.Client:
        if self._vos_client is None:
            self._vos_client = vos.Client(root_node=self.root_node)
        return self._vos_client

    def open(self, uri, cutout=None) -> vos.VOFile:
        parts = urlparse(uri)
        view = 'data'
        if cutout is not None:
            view = 'cutout'
        return self.vos_client.open(parts.path, view=view, cutout=cutout)


class VOSpaceApcorDownloader(VOSpaceDownloader):

    @staticmethod
    def _cast_result(result):
        apcor = result.read().decode('utf8').split('\n')[0].split()
        assert len(apcor) == 4
        apcor = [float(x) for x in apcor]
        return {'rin': apcor[0], 'rout': apcor[1], 'apcor': apcor[2], 'apcor_err': apcor[3]}


class AsyncVOSpaceFITSCutoutDownloader(VOSpaceDownloader):
    """
    Information on a FITS image cutout from the VOSpace storage system.
    """

    def __init__(self, minimum_cutout_radius: units.Quantity = 10*units.arcsec, **kwargs):
        super().__init__(**kwargs)
        self.minimum_cutout_radius = minimum_cutout_radius

    async def __call__(self, uri, centre_point, cutout_radius):
        self.check_schema(uri)
        cutout_string = self.cutout_string(centre_point, cutout_radius)
        result = await asyncio.to_thread(self.get_bytes_from_vospace, uri, cutout_string)
        return result

    def cutout_string(self, centre_point, cutout_radius) -> str:
        radius = max(cutout_radius, self.minimum_cutout_radius)
        return (f"CIRCLE ICRS {centre_point.ra.to('degree').value:.4f} "
                f"{centre_point.dec.to('degree').value:.4f} {radius.to('degree').value:0.4f}")

    @staticmethod
    def _cast_result(fits_obj) -> fits.HDUList:
        _hdulist = fits.open(fits_obj)
        if len(_hdulist) == 1:
            _hdulist = fits.HDUList([fits.PrimaryHDU(),
                                     fits.ImageHDU(data=_hdulist[0].data,
                                                   header=_hdulist[0].header)])
        return _hdulist


class VOSpaceFITSCutoutDownloader(AsyncVOSpaceFITSCutoutDownloader):
    def __call__(self, uri, centre_point, cutout_radius):
        self.check_schema(uri)
        cutout_string = self.cutout_string(centre_point, cutout_radius)
        return self.get_bytes_from_vospace(uri, cutout_string)