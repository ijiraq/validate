"""
The validate app
"""
import math
import sys

from astropy.time import Time
from mp_ephem import ObsRecord
from photutils.aperture import aperture_photometry, CircularAperture, CircularAnnulus, ApertureStats
from .display import Display
import logging
from .wcs import WCS

band_filter_map = {'gri.MP9605': 'w'}


class ObservationValidator:

    def __init__(self, artifacts, apcor_downloader, comparison_finder):
        self.comparison_finder = comparison_finder
        self.apcor_downloader = apcor_downloader
        self.artifacts = artifacts
        self.display = Display()
        self.frames = {}
        self._cursor = None
        self.obs_records = {}

    def display_images(self):
        for artifact in self.artifacts:
            self.frames[self.display.load(artifact.hdulist)] = artifact
            self.display.mark_ellipse(artifact.predicted_coordinate,
                                      *artifact.ephemeris_uncertainty,
                                      colour='yellow')

    def measure_current_source(self):
        recentroid = (self.cursor['key'] == 'a')
        frameno = self.cursor['frameno']
        x = self.cursor['x']
        y = self.cursor['y']
        extname = self.cursor['extname']
        artifact = self.frames[frameno]
        hdu = artifact.hdulist[extname]
        apcor = self.apcor_downloader(artifact.get_apcor_uri(extname))
        x, y, mag, merr = do_phot(hdu.data,
                                  hdu.header['PHOTZP'], x, y,
                                  recentroid=recentroid, **apcor)
        self.display.mark_annuli(x, y, annuli=[apcor['rin'], apcor['rout'] + apcor['rin'],
                                               apcor['rout'] + 2 * apcor['rin']])
        ra, dec = WCS(hdu.header).xy2sky(x, y, usepv=True)
        self.obs_records[artifact.observation_id] = (
            build_obs_record(hdu.header, mag=mag, merr=merr, x=x, y=y, ra=ra, dec=dec))

    def display_help(self):
        for key in self.commands:
            sys.stderr.write(f"{key}  -> {self.commands[key][0]}\n")

    @property
    def commands(self):
        return {'q': ['Quit', ],
                'a': ['Accept, recentroid', self.measure_current_source],
                'm': ['Accept, do not recentroid', self.measure_current_source],
                'n': ['Next Frame', self.display.next_frame],
                'p': ['Previous Frame', self.display.previous_frame],
                'b': ['Load blank', self.load_blank],
                '?': ['Help', self.display_help]
                }

    def load_blank(self):
        frameno = self.cursor['frameno']
        artifact = self.frames[frameno]
        if artifact.comparison is None:
            artifact.comparison = self.comparison_finder(artifact)
            if artifact.comparison is None:
                logging.warning(f"Could not find comparison for {artifact.observation_id}")
                return
            self.frames[self.display.load(artifact.comparison.hdulist)] = artifact.comparison
        self.goto_frame(artifact.comparison)

    def goto_frame(self, artifact):
        for frameno in self.frames:
            if self.frames[frameno] == artifact:
                self.display.goto_frame(frameno)
                return
        raise ValueError("Could not find frame for artifact.")

    def invalid_key(self):
        logging.warning(f"Unknown command: {self.cursor['key']}")
        self.display_help()

    @property
    def cursor(self):
        if self._cursor is None:
            self._cursor = self.display.get_cursor()
        return self._cursor

    def reset_cursor(self):
        self._cursor = None

    def create_obs_records(self) -> list[ObsRecord]:
        self.display_images()
        while self.cursor['key'] != 'q':
            self.commands.get(self.cursor['key'], [None, self.invalid_key])[1]()
            self.reset_cursor()
        return self.obs_records


def build_obs_record(header, mag, merr, x, y, ra, dec):
    """
    Build an observation record from the given artifact and photometry results.
    """
    date = Time((header['MJDEND'] + header['MJD-OBS']) / 2.0, format='mjd', precision=6).mpc
    frame = f"{header['EXPNUM']}p{header['EXTVER']:02d}"
    plate_uncertainty = header.get('ASTERR', None)
    astrometric_level = header.get('ASTLEVEL', None)
    band = band_filter_map.get(header['FILTER'], 'r')
    return ObsRecord(provisional_name=header['OBJECT'],
                     ra=ra,
                     dec=dec,
                     xpos=x,
                     ypos=y,
                     mag=mag,
                     band=band,
                     mag_err=merr,
                     date=date,
                     frame=frame,
                     plate_uncertainty=plate_uncertainty,
                     astrometric_level=astrometric_level,
                     survey_code='O',
                     observatory_code=568
                     )


def do_phot(data, zero_point, x, y, recentroid, rin, rout, apcor, apcor_err):
    """
    Perform aperture photometry on the image at the given coordinates.
    """
    source_mag = source_mag_err = None
    positions = [x-1, y-1]
    data_err = data**0.5
    try:
        sky_annulus = CircularAnnulus(positions, r_in=rout+rin, r_out=rout+2*rin)
        sky_counts = ApertureStats(data, sky_annulus).mode

        if recentroid:
            source_aperture = CircularAperture(positions, r=2 * rin)
            positions = ApertureStats(data-sky_counts, source_aperture).centroid

        source_aperture = CircularAperture(positions, r=rin)
        phot_result = aperture_photometry(data-sky_counts, source_aperture, error=data_err)
        source_flux = phot_result['aperture_sum'][0]
        source_err = phot_result['aperture_sum_err'][0]

        source_mag = -2.5*math.log10(source_flux) + zero_point - apcor
        source_mag_err = ((2.5/math.log(10)*source_err/source_flux)**2 + apcor_err**2)**0.5
    except Exception as ex:
        logging.error("Failed getting magnitude: {}".format(ex))
    return positions[0]+1, positions[1]+1, source_mag, source_mag_err

