import matplotlib.pyplot as plt
import numpy as np

from raytracer.elements import SphericalRefraction
from raytracer.elements import OutputPlane

from raytracer.rays import RayBundle            

class PlanoConvex:
    """Defines Plano-Convex Optical Element"""
    def __init__(self, z_0=100.0, curvature1=0.02, curvature2=0.0, n_inside=1.5168, n_outside=1, thickness=5, aperture=50):
        """Initializes the parameters defining the Plano-Convex optical element"""
        self._z_0 = z_0
        self._curvature1 = curvature1
        self._curvature2 = curvature2
        self._n_inside = n_inside
        self._n_outside = n_outside
        self._thickness = thickness
        self._aperture = aperture

        self.wavelength = 0.000588

    def propagate_ray(self, ray):
        if self._curvature1 > 0:
            plano_lens = SphericalRefraction(z_0=self._z_0, aperture=self._aperture, curvature=self._curvature1, n_1=self._n_outside, n_2=self._n_inside)
            convex_lens = SphericalRefraction(z_0=self._z_0 + self._thickness, aperture=self._aperture, curvature=0, n_1=self._n_inside, n_2=self._n_outside)
        elif self._curvature2 < 0:
            plano_lens = SphericalRefraction(z_0=self._z_0, aperture=self._aperture, curvature=0, n_1=self._n_outside, n_2=self._n_inside)
            convex_lens = SphericalRefraction(z_0=self._z_0 + self._thickness, aperture=self._aperture, curvature=self._curvature2, n_1=self._n_inside, n_2=self._n_outside)

        plano_lens.propagate_ray(ray)
        convex_lens.propagate_ray(ray)

    def focal_point(self):
        z_0 = self._z_0
        n_1 = self._n_inside
        n_2 = self._n_outside

        if self._curvature1 > 0:
            R_1 = 1 / self._curvature1
            f = R_1/(self._n_inside-1)
            BFL = f * (1 - ((self._n_inside - 1) * self._thickness) / ((1/self._curvature1) * self._n_inside))
            location_f = z_0 + self._thickness
            return location_f + BFL

        elif self._curvature2 < 0:
            R_2 = 1 / self._curvature2
            f_reciprocal = (n_1 / n_2 - 1) * (-1 / R_2)
            f = 1 / f_reciprocal
            location_f = z_0 + f + self._thickness
            return location_f

        else:
            raise ValueError("Lock in. Make sure the curvature has one 0 and the other  the correct + or -")


