import kete

from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table

from AsteroidDB import default_db_file


class AsteroidChecker:

    def __init__(self):
        
        # TODO: Establish a connection to the database.
        self.load_states_from_db()
    
    def load_states_from_db(self, filename=default_db_file):

        # TODO: Load this from a database instead. Can do this for only some ra and dec.
        try:
            self.state_table = Table.read(filename)
        except FileNotFoundError:
            print(f"File {filename} not found. Make sure to download the orbit data first.")

        # Convert the table to a list of states.
        self.mpc_states = []
        for row in self.state_table:
            state = kete.State(row['designations'], row['jd'], kete.vector.Vector(row['position']), kete.vector.Vector(row['velocity']), frame=kete.vector.Frames.Equatorial)
            self.mpc_states.append(state)


    def get_asteroid_list(self, frame_wcs, time_jd):

        earth = kete.spice.get_state("Earth", time_jd)

        fov = kete.fov.RectangleFOV.from_wcs(frame_wcs, earth)

        # Update the states based on the time. 
        self.mpc_states = kete.propagate_n_body(self.mpc_states, time_jd)

        # Check the FOV for asteroids.
        asteroid_list = kete.fov_state_check(self.mpc_states, [fov])

        return asteroid_list

    def cross_match_sources_with_asteroids(self, asteroid_list, sources, radius_threshold):
        
        raise NotImplementedError

if __name__=="__main__":

    # Load a fits
    frame = fits.open("ngc2403_V_dr4.fits")[0]
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

    jd = kete.Time.from_ymd(2025, 2, 24).jd
    print(ac.get_asteroid_list(frame_wcs, jd))

    end_time = time.perf_counter()
    elapsed_time = end_time - start_time
    print(f"Elapsed time for finding asteroid in field: {elapsed_time:.4f} seconds")









