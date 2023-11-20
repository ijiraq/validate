from astropy import units
from astropy.coordinates import SkyCoord
from astropy.units import Quantity


class Region(object):
    """
    A DS9 region object.

    given a 'point' creates a circle region at that point.

    style can be one of 'cirlce', 'ellipse', 'annulus', 'point'

    if 'circle' then the radius of the circle is passed as the shape.
    if 'ellipse' then an ellipse object should be passed as the 'shape'
    if 'annulus' then a set of annulus sizes should be passed as the 'shape'
    if 'point' then the radius is ignored. point(1:37:02.476,+12:35:42.12) # point=x

    """

    def __init__(self, point, style='circle', colour='g', shape=10, option=""):
        """
        :param shape: The parameters that define the shape.
        :type shape: int, list
        """
        self._point = point
        self._coo = None
        self._colour = None
        self.colour = colour
        self.style = style
        self.shape = shape
        self.option = option

    def __str__(self):
        s = ""
        for c in self:
            s += str(c)+" "
        return s

    def __iter__(self):
        if self.style == 'point':
            return iter(('regions', '{}; {}({},{}) # color={} point=x '.format(self.coosys,
                                                                               self.style,
                                                                               self.point[0],
                                                                               self.point[1],
                                                                               self.colour)))
        try:
            if isinstance(self.shape, Ellipse):
                shape = self.shape
            elif isinstance(self.shape[0], Quantity):
                shape = ",".join(['{}"'.format(p.to(units.arcsec).value) for p in self.shape])
            else:
                shape = ",".join(['{}'.format(p) for p in self.shape])
        except Exception:
            try:
                shape = '{}"'.format(self.shape.to(units.arcsec).value)
            except:
                shape = '{}"'.format(0.185*self.shape)

        return iter(('regions', '{}; {}({},{},{}) # color={} {}'.format(self.coosys,
                                                                        self.style,
                                                                        self.point[0],
                                                                        self.point[1],
                                                                        shape,
                                                                        self.colour,
                                                                        self.option)))

    @property
    def point(self):
        if not self._coo:
            if isinstance(self._point, SkyCoord):
                self._coo = self._point.ra.degree, self._point.dec.degree
            elif isinstance(self._point[0], Quantity):
                try:
                    self._coo = self._point[0].to(units.degree).value, self._point[1].to(units.degree).value
                except:
                    self._coo = self._point[0].value, self._point[1].value
            else:
                self._coo = self._point
        return self._coo

    @property
    def coosys(self):
        if isinstance(self._point, SkyCoord):
            return "wcs"
        if isinstance(self._point[0], Quantity):
            try:
                self._point[0].to(units.degree)
                return "wcs"
            except:
                return "image"
        return "image"

    @property
    def colour(self):
        """
        The colour of the marker to create.
        """
        return self._colour

    @colour.setter
    def colour(self, colour):
        self._colour = {'r': 'red', 'b': 'blue', 'y': 'yellow', 'c': 'cyan', 'g': 'green'}.get(colour, colour)


class Ellipse(object):
    """
    An ellipse region for use in DS9.

    a = semi-major axis
    b = semi-minor axis
    pa = position angle in degrees (East is 0).

    """

    def __init__(self, a, b, pa):
        """
        :param a: semi-major axis
        :type a: Quantity
        :param b: semi-minor axis
        :type b: Quantity
        :param pa: Position Angle (East is 0)
        :type pa: Quantity
        """
        assert isinstance(a, Quantity)
        assert isinstance(b, Quantity)
        assert isinstance(pa, Quantity)
        self._a = a
        self._b = b
        self._pa = pa

    def __str__(self):
        return '{}", {}", {}'.format(self.a.to(units.arcsec).value,
                                     self.b.to(units.arcsec).value,
                                     self.pa.to(units.degree).value + 90)

    @property
    def a(self):
        """
        The semi-major axis of the ellipse
        @return: Quantity
        """
        return self._a

    @a.setter
    def a(self, value):
        self._a = value

    @property
    def b(self):
        """
        The semi-minor axis of the ellipse.
        @return: Quantity
        """
        return self._b

    @b.setter
    def b(self, value):
        self._b = value

    @property
    def pa(self):
        """
        The orientation angle of the ellipse, counter-clock wise from due west.
        @return: Quantity
        """
        return self._pa

    @pa.setter
    def pa(self, value):
        self._pa = value
