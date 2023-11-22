"""
The validate app
"""
import math
import sys
from astropy.time import Time
from mp_ephem import ObsRecord, BKOrbit
from photutils.aperture import aperture_photometry, CircularAperture, CircularAnnulus, ApertureStats
from .display import Display
import logging
from .wcs import WCS

band_filter_map = {'gri.MP9605': 'w'}


class ObservationValidator:

    def __init__(self, artifacts, apcor_downloader, comparison_finder, obs_records):
        self.comparison_finder = comparison_finder
        self.apcor_downloader = apcor_downloader
        self.obs_records = obs_records
        self.update_trial_orbit()
        self.artifacts = artifacts
        self.display = Display()
        self._cursor = None
        for i, artifact in enumerate(self.artifacts):
            artifact.previous_artifact = self.artifacts[i-1]
            artifact.next_artifact = self.artifacts[(i+1) % len(self.artifacts)]
        self.current_artifact = self.artifacts[0]
        self.display_current_artifact()

    def display_current_artifact(self):
        artifact = self.current_artifact
        hdulist = artifact.hdulist
        try:
            self.display.load(hdulist)
            self.trial_orbit.predict(Time(artifact.obs_date, format='mjd'))
            self.display.focus = self.trial_orbit.coordinate
            logging.debug(f"{artifact.observation_id} -> {artifact.obs_record}")
            colour = ( 'yellow' if artifact.obs_record is None
                       else 'green' if not artifact.obs_record.null_observation
            else 'red' )
            self.display.mark_ellipse(self.trial_orbit.coordinate,
                                      self.trial_orbit.dra,
                                      self.trial_orbit.ddec,
                                      self.trial_orbit.pa,
                                      colour=colour)
        except Exception as ex:
            logging.error(f"Failed to display {hdulist} : {ex}")

    def display_next_artifact(self):
        self.current_artifact = self.current_artifact.next_artifact
        self.display_current_artifact()

    def display_previous_artifact(self):
        self.current_artifact = self.current_artifact.previous_artifact
        self.display_current_artifact()

    def measure_current_source(self):
        artifact = self.current_artifact
        recentroid = (self.cursor['key'] == 'a')
        null_observation = (self.cursor['key'] == 'r')
        x = self.cursor['x']
        y = self.cursor['y']
        extname = self.cursor['extname']
        hdu = artifact.hdulist[extname]
        apcor = {'rin': 5, 'rout': 15, 'apcor': None, 'apcor_err': None}
        try:
            apcor = self.apcor_downloader(uri=artifact.get_apcor_uri(extname))
        except Exception as ex:
            logging.warning(f"Failed to download apcor for {artifact.observation_id}: {ex}")
        x, y, mag, merr = do_phot(hdu.data,
                                  hdu.header['PHOTZP'], x, y,
                                  recentroid=recentroid, **apcor)
        self.display_current_artifact()
        self.display.mark_annuli(x, y, annuli=[apcor['rin'], apcor['rout'] + apcor['rin'],
                                               apcor['rout'] + 2 * apcor['rin']])
        ra, dec = WCS(hdu.header).xy2sky(x, y, usepv=True)
        new_obs_record = build_obs_record(hdu.header, mag=mag, merr=merr, x=x, y=y, ra=ra, dec=dec,
                                          null_observation=null_observation)
        artifact.obs_record = new_obs_record
        self.obs_records[f"{artifact.obs_date:.6g}"] = new_obs_record
        self.update_trial_orbit()

    def update_trial_orbit(self):
        self.trial_orbit = BKOrbit(list(self.obs_records.values()))
        print(self.trial_orbit.summarize())


    def display_help(self):
        for key in self.commands:
            sys.stderr.write(f"{key}  -> {self.commands[key][0]}\n")

    @property
    def commands(self):
        return {'q': ['Quit', ],
                'a': ['Accept, recentroid', self.measure_current_source],
                'm': ['Accept, do not recentroid', self.measure_current_source],
                'n': ['Next Frame', self.display_next_artifact],
                'p': ['Previous Frame', self.display_previous_artifact],
                'c': ['Load comparison', self.display_comparison_artifact],
                'r': ['Reject observation', self.measure_current_source],
                '?': ['Help', self.display_help]
                }

    def display_comparison_artifact(self):
        logging.info(f"Getting a comparison image.")
        artifact = self.current_artifact
        if artifact.comparison is None:
            artifact.comparison = self.comparison_finder(artifact)
            if artifact.comparison is None:
                logging.warning(f"Could not find comparison for {artifact.observation_id}")
                return
            artifact.comparison.next_artifact = artifact.next_artifact
            artifact.comparison.previous_artifact = artifact
        self.current_artifact = artifact.comparison
        self.display_current_artifact()

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
        while self.cursor['key'] != 'q':
            self.commands.get(self.cursor['key'], [None, self.invalid_key])[1]()
            self.reset_cursor()
        return self.obs_records


def build_obs_record(header, mag, merr, x, y, ra, dec, null_observation=False):
    """
    Build an observation record from the given artifact and photometry results.
    """
    date = Time((header['MJDEND'] + header['MJD-OBS']) / 2.0, format='mjd', precision=6).mpc
    frame = f"{header['EXPNUM']}p{header['EXTVER']:02d}"
    plate_uncertainty = header.get('ASTERR', None)
    astrometric_level = header.get('ASTLEVEL', None)
    band = band_filter_map.get(header['FILTER'], 'r')
    return ObsRecord(provisional_name=header['OBJECT'],
                     null_observation=null_observation,
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

