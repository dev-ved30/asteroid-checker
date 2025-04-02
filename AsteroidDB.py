import time
import kete

from astropy.table import Table


default_db_file = "states.ecsv"


def download_and_update_db(time_jd, filename=default_db_file):

    # Load orbit data from the MPC
    mpc_obs = kete.mpc.fetch_known_orbit_data(force_download=True)

    # Convert that data to State objects.
    mpc_states = kete.mpc.table_to_states(mpc_obs)

    # Update the states based on the time. 
    mpc_states = kete.propagate_n_body(mpc_states, time_jd)

    # Convert the states from ecliptic to equatorial to easily access the RA and Dec.
    for i, state in enumerate(mpc_states):
        mpc_states[i] = state.as_equatorial

    # Save the states to a database.
    designations = []
    jd = []
    ra = []
    dec = []
    positions = []
    velocities = []

    for state in mpc_states:

        designations.append(state.desig)
        jd.append(state.jd)
        ra.append(state.pos.ra)
        dec.append(state.pos.dec)
        positions.append(state.pos) # Position of the object in AU with respect to the central object.
        velocities.append(state.vel) # Velocity of the object in AU/Day.

    # TODO: Replace this with storing to a database
    state_table = Table([designations, jd, ra, dec, positions, velocities],
        names=['designations', 'jd', 'ra', 'dec', 'position', 'velocity'])
    
    # TODO: Store this in a database instead.
    state_table.write(filename, format='ascii.ecsv', overwrite=True)

if __name__=="__main__":

    start_time = time.perf_counter()

    # Feb 4, 2025
    jd = kete.Time.from_ymd(2025, 2, 23).jd
    download_and_update_db(jd)
    
    end_time = time.perf_counter()

    print(f"Time taken: {end_time - start_time:.2f} seconds")