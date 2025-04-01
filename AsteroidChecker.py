import kete
import datetime
import astropy

import numpy as np

from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table
from memory_profiler import profile

default_db_file = "states.ecsv"


class AsteroidChecker:

    def __init__(self):
        
        self.load_states_from_db()
    
    def download_orbit_data_and_update_db(self):

        # Load orbit data from the MPC
        self.mpc_obs = kete.mpc.fetch_known_orbit_data(force_download=True)

        # Convert that data to State objects.
        self.mpc_states = kete.mpc.table_to_states(self.mpc_obs)

        # Convert the states from ecliptic to equatorial to easily access the RA and Dec.
        for i, state in enumerate(self.mpc_states):
            self.mpc_states[i] = state.as_equatorial

        # Save the states to a database.
        self.save_states_to_db()

    def save_states_to_db(self, filename=default_db_file):

        designations = []
        jd = []
        ra = []
        dec = []
        positions = []
        velocities = []

        for state in self.mpc_states:

            designations.append(state.desig)
            jd.append(state.jd)
            ra.append(state.pos.ra)
            dec.append(state.pos.dec)
            positions.append(state.pos) # Position of the object in AU with respect to the central object.
            velocities.append(state.vel) # Velocity of the object in AU/Day.

        # TODO: Replace this with storing to a database
        self.state_table = Table([designations, jd, ra, dec, positions, velocities],
           names=['designations', 'jd', 'ra', 'dec', 'position', 'velocity'])
        
        # TODO: Store this in a database instead.
        self.state_table.write(filename, format='ascii.ecsv', overwrite=True)

    def load_states_from_db(self, filename=default_db_file):

        # TODO: Load this from a database instead. Can do this for only some ra and dec.
        try:
            self.state_table = Table.read(filename)
        except FileNotFoundError:
            # If the file doesn't exist, download the orbit data and update the database.
            self.download_orbit_data_and_update_db()
            self.state_table = Table.read(filename)

        # Convert the table to a list of states.
        self.mpc_states = []
        for row in self.state_table:
            state = kete.State(row['designations'], row['jd'], kete.vector.Vector(row['position']), kete.vector.Vector(row['velocity']), frame=kete.vector.Frames.Equatorial)
            self.mpc_states.append(state)

    #@profile
    def propagate_asteroids(self, time_jd):
        """
        Propagate the asteroids to the some time.
        
        :param time_jd: The current time in JD.
        """

        # Update the states based on the time. 
        self.mpc_states = kete.propagate_n_body(self.mpc_states, time_jd)


    #@profile
    def get_asteroid_list(self, frame_wcs, time_jd):

        earth = kete.spice.get_state("Earth", time_jd)


        fov = kete.fov.RectangleFOV.from_wcs(frame_wcs, earth)

        # Update the states based on the time. 
        self.propagate_asteroids(time_jd)

        # Check the FOV for asteroids.
        asteroid_list = kete.fov_state_check(self.mpc_states, [fov])

        return asteroid_list

    def cross_match_sources_with_asteroids(self, asteroid_list, sources, radius_threshold):
        
        raise NotImplementedError

if __name__=="__main__":

    # Feb 4, 2025
    jd = kete.Time.from_ymd(2025, 2, 23).jd

    print(jd)

    # Load a fits
    frame = fits.open("ngc2403_V_dr4.fits")[0]
    frame_wcs = WCS(frame.header)

    import time
    start_time = time.perf_counter()
    # Code to be timed

    ac = AsteroidChecker()

    end_time = time.perf_counter()
    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time:.4f} seconds")


    start_time = time.perf_counter()
    # Code to be timed

    print(ac.get_asteroid_list(frame_wcs, jd))

    end_time = time.perf_counter()
    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time:.4f} seconds")

    ac.save_states_to_db()

    start_time = time.perf_counter()
    # Code to be timed

    print(ac.get_asteroid_list(frame_wcs, kete.Time.from_ymd(2025, 2, 24).jd))

    end_time = time.perf_counter()
    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time:.4f} seconds")







