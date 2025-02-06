import kete
import datetime
import astropy

import numpy as np

from astropy.io import fits
from astropy.wcs import WCS

class AsteroidChecker:

    def __init__(self):

        # Load orbit data from the MPC
        self.mpc_obs = kete.mpc.fetch_known_orbit_data()

        # Convert that data to State objects.
        self.mpc_states = kete.mpc.table_to_states(self.mpc_obs)

    def propagate_asteroids(self, time_jd):
        """
        Propagate the asteroids to the some time.
        
        :param time_jd: The current time in JD.
        """

        # Update the states based on the time. 
        self.mpc_states = kete.propagate_n_body(self.mpc_states, time_jd)

        # TODO: we need to find a way to save this data so we don't have to propagate it every time.


    def get_asteroid_list(self, frame_wcs, time_jd):

        corners = []
        dx, dy = frame_wcs.pixel_shape
        for x, y in zip([0, 0, dx, dx], [0, dy, dy, 0]):
            coord = frame_wcs.pixel_to_world(x, y).icrs
            corners.append(kete.Vector.from_ra_dec(coord.ra.deg, coord.dec.deg))

        fov = kete.fov.RectangleFOV.from_corners(corners, earth)

        # Update the states based on the time. 
        self.propagate_asteroids(time_jd)

        # Check the FOV for asteroids.
        asteroid_list = kete.fov_state_check(self.mpc_states, [fov], include_asteroids=True)[0]

        return asteroid_list

    def cross_match_sources_with_asteroids(self, asteroid_list, sources, radius_threshold):
        
        pass

if __name__=="__main__":

    # Feb 4, 2025
    jd = kete.Time.from_ymd(2025, 2, 6).jd
    earth = kete.spice.get_state("Earth", jd)

    # Load a fits
    frame = fits.open("ngc2403_V_dr4.fits")[0]
    frame_wcs = WCS(frame.header)



    ac = AsteroidChecker()
    print(ac.get_asteroid_list(frame_wcs, jd))




