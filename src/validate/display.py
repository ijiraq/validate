import contextlib
import logging
import time
import pyds9
from io import BytesIO
from astropy.coordinates import SkyCoord
from astropy import units
from astropy.io import fits
from astropy.units import Quantity
from .wcs import WCS
from .region import Region, Ellipse

CONFIG = {
    "INIT": {
        "bg": "black",
        "view_layout": "vertical",
        "view_info": "no",
        "view_magnifier": "no",
        "view_buttons": "no",
        "view_panner": "no",
        "view_filename": "no",
        "view_frame": "no",
        "view_physical": "no",
        "view_object": "no",
        "view_colorbar": "no",
        "width": 512,
        "height": 512,
        "scale_mode": "zscale",
        "zoom to": 4,
        "frame": "delete all"
    },
    "PREF": {
        "frame_lock": "wcs",
        "scale": "linear",
        "scale_mode": "zscale",
        "cmap_lock": "yes",
        "cmap_invert": "yes",
        "smooth": "no"
    }
}


class DisplayError(Exception):
    pass


class Display(object):

    def __init__(self):
        self.minimum_pan = 10 * units.arcsec
        while True:
            try:
                self.ds9 = pyds9.DS9('validate')
                self.ds9.set('frame delete all')
                self.ds9.set('frame new')
                self.ds9.set('zoom to 4')
                self.ds9.set('wcs align yes')
                break
            except ValueError as ex:
                logging.critical(f"Waiting for ds9 to start: {ex}")
                time.sleep(2)
        self.current_centre = SkyCoord(0,0, unit=units.degree)
        self.initialize()

    def next_frame(self):
        self.ds9.set('frame next')

    def previous_frame(self):
        self.ds9.set('frame previous')

    def goto_frame(self, frameno):
        self.ds9.set(f'frame {frameno}')

    @property
    def frameno(self) -> int:
        return int(self.ds9.get('frame frameno'))

    def initialize(self):
        [self.ds9.set(f"{cmd.replace('_', ' ')} {CONFIG['INIT'][cmd]}") for cmd in CONFIG['INIT']]

    def load(self, hdulist) -> int:
        """
        Display the image in ds9.
        """
        if hdulist is not None:
            # Get centre of
            with (contextlib.closing(BytesIO()) as newFitsFile):
                hdulist.writeto(newFitsFile)
                temp_fits = newFitsFile.getvalue()
                self.ds9.set('fits mosaicimage wcs', temp_fits, len(temp_fits))
        else:
            self.ds9.set('frame new')
        return int(self.ds9.get('frame frameno'))

    def load(self, hdulist) -> int:
        """
        Display the image in ds9.
        """
        if hdulist is not None:
            with contextlib.closing(BytesIO()) as newFitsFile:
                hdulist.writeto(newFitsFile)
                newfits = newFitsFile.getvalue()
                got = self.ds9.set('fits new mosaicimage wcs', newfits, len(newfits))
                self.ds9.set('zoom to 4')
                self.ds9.set('wcs align yes')
        else:
            self.ds9.set('frame new')
        return int(self.ds9.get('frame frameno'))

    def pan(self, ra, dec):
        """
        Pan the image in ds9.
        """
        return self.ds9.set(f'pan to {ra} {dec} wcs icrs')

    def aligned(self):
        if self.focus is None:
            return False
        focus = self.ds9.get("pan wcs degrees").split()
        focus = SkyCoord(focus[0], focus[1], unit=units.degree)
        return focus.separation(self.focus) < self.minimum_pan

    @property
    def focus(self):
        return self._focus

    @focus.setter
    def focus(self, focus):
        try:
            if not isinstance(focus, SkyCoord):
                focus = SkyCoord(focus[0], focus[1])
            if self.focus is None or focus.separation(self.focus) > self.minimum_pan:
                self._focus = focus
            self.ds9.set("pan to {} {} wcs fk5".format(self.focus.ra.degree,
                                                       self.focus.dec.degree))
        except Exception as ex:
            logging.error("Focus setting failed {focus}")
            self._focus = None

    def _do_move_focus(self):
        if self.focus is None:
            return
        self.ds9.set("pan to {} {} wcs fk5".format(self.focus.ra.degree,
                                                   self.focus.dec.degree))

    def pan_to(self, pos):
        self.focus = pos
        self._do_move_focus()

    def get_cursor(self) -> dict:
        """
        send back the cursor ds9 pointer location and key pressed as a dictionary:

        'key', 'frameno', 'x', 'y', 'extname', 'ra', 'dec'
        :return:
        """
        while True:
            try:
                values = self.ds9.get('iexam key coordinate wcs icrs degrees').split()
                if len(values) != 3:
                    raise ValueError("Wrong number of values returned by get.")
                break
            except Exception as ex:
                logging.debug("Error reading from ds9: {}".format(type(ex)))
        key = values[0]
        ra = Quantity(float(values[1]), unit=units.degree)
        dec = Quantity(float(values[2]), unit=units.degree)
        cursor = {'key': key, 'frameno': self.frameno, 'ra': ra, 'dec': dec}
        cursor.update(self.resolve_mosaic_extension(ra, dec))
        return cursor

    def resolve_mosaic_extension(self, ra: Quantity, dec: Quantity) -> dict:
        """determine which of the extension in a mosaic of extensions the selected source is in and return
        the extension name and x/y location of source."""
        try:
            for header in self.headers:
                this_wcs = WCS(header)
                x, y = this_wcs.sky2xy(ra, dec, usepv=True)
                contains_pixel = 0 < x < header['NAXIS1'] and 0 < y < header['NAXIS2']
                if contains_pixel:
                    return {'extname': header['EXTNAME'], 'x': x, 'y': y}
        except Exception as ex:
            logging.debug(f"Error resolving mosaic extension: {type(ex)} -> {str(ex)}")
        logging.warning("Could not match RA/DEC returned by iexem to value inside image.")
        return {'extname': None, 'x': None, 'y': None}

    @property
    def headers(self) -> list[fits.Header]:
        """Get the headers from the currently selected ds9 frame."""
        headers = []
        while True:
            try:
                header_str = self.ds9.get(f'fits header {len(headers) + 1}')
                headers.append(fits.Header.fromstring(header_str, sep='\n'))
            except TypeError as e:
                logging.debug(f"{type(e)} -> {str(e)}")
                break
        return headers

    def mark_annuli(self, x, y, annuli, colour='b'):
        r = Region((x, y), style='annulus', colour=colour, shape=annuli)
        self.ds9.set(*r)

    def mark_cirlce(self, position, radius, colour='b'):
        r = Region(position, style='circle', colour=colour, shape=radius)
        self.ds9.set(*r)

    def mark_ellipse(self, position, a, b, angle, colour='b'):
        r = Region(position, style='ellipse', colour=colour, shape=Ellipse(a, b, angle))
        self.ds9.set(*r)
