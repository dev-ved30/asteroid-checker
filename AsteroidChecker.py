import kete
import psycopg2

from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table

from AsteroidDB import db_params


class AsteroidChecker:

    def __init__(self):
        
        # Load the database table and convert to states
        self.load_states_from_db()
    
    def load_states_from_db(self, db_params=db_params):

        try:
            conn = psycopg2.connect(**db_params)
            cursor = conn.cursor()

            # TODO See if we can apply cuts based on ra/dec and proper motion
            cursor.execute("SELECT * FROM asteroids;")
            rows = cursor.fetchall()


            # Convert rows into sate objects
            self.mpc_states = []
            for row in rows:
                
                designation, jd, ra, dec, position_x, position_y, position_z, velocity_x, velocity_y, velocity_z = row
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
        asteroid_list = kete.fov_state_check(self.mpc_states, [fov], dt_limit=0)

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









