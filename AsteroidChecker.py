import kete
import psycopg2

import astropy.units as u

from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy.io import fits
from astropy.wcs import WCS

from AsteroidDB import db_params

default_radius_threshold = 30 * u.arcsec

class AsteroidChecker:

    def __init__(self):
        
        # Load the database table and convert to states
        self.load_states_from_db()
    
    def load_states_from_db(self, db_params=db_params):

        try:
            conn = psycopg2.connect(**db_params)
            cursor = conn.cursor()

            # TODO See if we can apply cuts based on ra/dec and proper motion
            # Only select objects such that ra/dec + proper motion is within the field of view.
            cursor.execute("SELECT * FROM asteroids;")
            rows = cursor.fetchall()

            # Convert rows into sate objects
            self.mpc_states = []
            for row in rows:
                
                designation, jd, ra, dec, pm_ra, pm_dec, position_x, position_y, position_z, velocity_x, velocity_y, velocity_z = row
                state = kete.State(designation, jd, kete.vector.Vector([position_x, position_y, position_z]), \
                                kete.vector.Vector([velocity_x, velocity_y, velocity_z]))
                self.mpc_states.append(state)
            
            # Close the connection
            cursor.close()
            conn.close()

        except Exception as e:
            print("Error connecting to database:", e)

    def get_asteroid_list(self, frame_wcs, time_jd):

        earth = kete.spice.get_state("Earth", time_jd)

        fov = kete.fov.RectangleFOV.from_wcs(frame_wcs, earth)

        # Update the states based on the time. 
        self.mpc_states = kete.propagate_n_body(self.mpc_states, time_jd)

        # Check the FOV for asteroids.
        vis = kete.fov_state_check(self.mpc_states, [fov], dt_limit=0, include_asteroids=True)

        ra_list = []
        dec_list = []
        designation_list = []

        if len(vis) > 0:

            vis = vis[0]

            for idx in range(len(vis)):
                
                vec = vis.obs_vecs[idx]

                designation_list.append(vis[idx].desig)
                ra_list.append(vec.as_equatorial.ra)
                dec_list.append(vec.as_equatorial.dec)

            asteroid_sky_coords = SkyCoord(ra=ra_list, dec=dec_list, unit='deg')

        else:

            asteroid_sky_coords = []

        return asteroid_sky_coords, designation_list

    def cross_match_sources_with_asteroids(self, asteroid_sky_coords, asteroid_designations, sources_sky_coords, radius_threshold=default_radius_threshold):

        matched_asteroids = []
        idx, d2d, _ = sources_sky_coords.match_to_catalog_sky(asteroid_sky_coords)

        # Process and display the cross-match results
        for source, match_index, separation in zip(sources_sky_coords, idx, d2d):

            if separation < radius_threshold:
                matched_asteroids.append(asteroid_designations[match_index])
            else:
                matched_asteroids.append(None)
        
        return matched_asteroids


if __name__=="__main__":

    # Load a fits
    frame = fits.open("rings.v3.skycell.1572.032.stk.g.unconv.fits")[1]
    frame_wcs = WCS(frame.header)

    import time
    start_time = time.perf_counter()
    # Code to be timed

    ac = AsteroidChecker()

    end_time = time.perf_counter()
    elapsed_time = end_time - start_time
    print(f"Elapsed time for loading states: {elapsed_time:.4f} seconds")


    # Timing the code 
    start_time = time.perf_counter()

    jd = kete.Time.from_ymd(2025, 2, 23).jd
    asteroid_sky_coords, designations = ac.get_asteroid_list(frame_wcs, jd)
    print(asteroid_sky_coords)

    end_time = time.perf_counter()
    elapsed_time = end_time - start_time
    print(f"Elapsed time for finding asteroid in field: {elapsed_time:.4f} seconds")

    # Try to cross match with some sources
    print(ac.cross_match_sources_with_asteroids(asteroid_sky_coords, designations, asteroid_sky_coords))







