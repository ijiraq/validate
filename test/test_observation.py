import os
from unittest import TestCase

from astropy.io import fits
from astropy.table import Table
from validate.artifact import Artifact, observation_from_mocas_row


class Test(TestCase):

    def setUp(self) -> None:
        self.mocas_table = Table.read('data/mocas_test_result.csv', format='ascii.csv', delimiter='\t')

    def test_observation_from_mocas_row(self):
        for row in self.mocas_table:
            observation = observation_from_mocas_row(row)
            self.assertAlmostEqual(observation.coordinate.ra.degree, row['Object_RA'], places=4)
            self.assertAlmostEqual(observation.coordinate.dec.degree, row['Object_Dec'], places=4)

    def test_uri(self):
        observation = observation_from_mocas_row(self.mocas_table[0])
        self.assertEqual(observation.uri,
                         f"{observation.root}/{observation.observation_id}/{observation.observation_id}p.fits")

    def test_hdulist(self):
        observation = observation_from_mocas_row(self.mocas_table[0])
        self.assertIsInstance(observation.hdulist, fits.hdu.hdulist.HDUList)
        self.assertEqual(len(observation.hdulist), 1)
