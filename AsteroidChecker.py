import os
os.environ["OMP_NUM_THREADS"] = "1"

import kete
import psycopg2

import numpy as np
import astropy.units as u

from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy.io import fits
from astropy.wcs import WCS
from memory_profiler import profile

from AsteroidDB import db_params

default_radius_threshold = 30 * u.arcsec
chip_fov_radius = 1 * u.deg


import psutil

proc = psutil.Process(os.getpid())

print(f"\nOS‑level threads: {proc.num_threads()}")

class AsteroidChecker:

    #@profile
    def __init__(self):
        
        # Load the database table and convert to states
        self.load_states_from_db()
    
    #@profile
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

    #@profile
    def get_asteroid_list(self, frame_wcs, time_jd):

        def _kete_state_ra_dec(state, sun2earth):
            """
            Private method to wrap computing the coordinates of a single state
            """

            # change the coordinate frame from heliocentric to geocentric
            # then compute the geocentric coordinate
            obj_earth_pos = state.pos - sun2earth 
            state_vec_equitorial = obj_earth_pos.change_frame(kete.vector.Frames.Equatorial)
            return state_vec_equitorial.ra, state_vec_equitorial.dec


        sun2earth = kete.spice.get_state("Earth", time_jd)

        fov = kete.fov.RectangleFOV.from_wcs(frame_wcs, sun2earth)
        fov_center_ra, fov_center_dec = fov.pointing.as_equatorial.ra, fov.pointing.as_equatorial.dec
        target_skycoord =  SkyCoord(fov_center_ra, fov_center_dec, unit="deg")

        print(fov_center_ra, fov_center_dec)

        # Update the states based on the time using the two body approximation
        new_mpc_states_2body = kete.propagate_two_body(self.mpc_states, time_jd)
        ras, decs = np.array([_kete_state_ra_dec(state, sun2earth.pos) for state in new_mpc_states_2body]).T

        # Only get MPC matches that are withing 2*survey fov radius 
        mpc_skycoords_approx = SkyCoord(ras, decs, unit="deg")
        idxs = np.where(target_skycoord.separation(mpc_skycoords_approx) < chip_fov_radius)[0]

        # Now use the full n-body solution for all of these objects
        objs_win_first_cut = np.array(self.mpc_states)[idxs]
        print(f"Found {len(objs_win_first_cut)} w/in {2 * chip_fov_radius}, running full n-body on those...")

        states = kete.propagate_n_body(objs_win_first_cut, time_jd)
        ras, decs = np.array([_kete_state_ra_dec(state, sun2earth.pos) for state in states]).T
        mpc_skycoords = SkyCoord(ras, decs, unit="deg")
        designations_list = [s.desig for s in objs_win_first_cut]

        return mpc_skycoords, designations_list

    #@profile
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

    jd = 2460730.0 #kete.Time.from_ymd(2025, 2, 23).jd
    asteroid_sky_coords, designations = ac.get_asteroid_list(frame_wcs, jd)

    print(asteroid_sky_coords)

    new_ra  = 150.123456   # degrees
    new_dec =   26.54321   # degrees

    # Update only CRVAL (leaving CRPIX unchanged)
    frame_wcs.wcs.crval = [new_ra, new_dec]
    asteroid_sky_coords, designations = ac.get_asteroid_list(frame_wcs, jd)


    print(asteroid_sky_coords)

    end_time = time.perf_counter()
    elapsed_time = end_time - start_time
    print(f"Elapsed time for finding asteroid in field: {elapsed_time:.4f} seconds")

    # Try to cross match with some sources
    print(ac.cross_match_sources_with_asteroids(asteroid_sky_coords, designations, asteroid_sky_coords))



    mem = proc.memory_info()
    print(f"RSS: {mem.rss / (1024**2):.2f} MiB")
    print(f"VMS: {mem.vms / (1024**2):.2f} MiB")

    print("CPU % (1s):", proc.cpu_percent(interval=1))  # ≤ 100%




