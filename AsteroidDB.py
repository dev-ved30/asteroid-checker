import time
import kete
import psycopg2


# Database connection parameters
db_params = {
    "dbname": "asteroid",
    "user": "vedshah",  # Change to your PostgreSQL username
    "host": "localhost",  # Change if PostgreSQL is hosted remotely
    "port": "5432"  # Default PostgreSQL port
}

# <-- Units for all the fields in the database -->
# Designation: Text
# JD: Float (days)
# RA: Float (degrees)
# Dec: Float (degrees)
# PM_RA: Float (degrees/day)
# PM_Dec: Float (degrees/day)
# Position: Float (AU)
# Velocity: Float (AU/day)

def create_database_table(db_params=db_params):
    """
    Create a database table to store the asteroid data.
    """
    try:
        # Establish connection
        conn = psycopg2.connect(**db_params)
        cursor = conn.cursor()
        
        # TODO: Crate a temporary table here.
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS asteroids (
                designation TEXT PRIMARY KEY,
                jd FLOAT NOT NULL,
                ra FLOAT NOT NULL,
                dec FLOAT NOT NULL,
                pm_ra FLOAT NOT NULL,
                pm_dec FLOAT NOT NULL,
                position_x FLOAT NOT NULL,
                position_y FLOAT NOT NULL,
                position_z FLOAT NOT NULL,
                velocity_x FLOAT NOT NULL,
                velocity_y FLOAT NOT NULL,
                velocity_z FLOAT NOT NULL
            );
        """)
        conn.commit()

        # Close the connection
        cursor.close()
        conn.close()

    except Exception as e:

        print("Error creating database table:", e)


def download_and_update_db(time_jd, db_params=db_params):

    # Establish connection
    conn = psycopg2.connect(**db_params)
    cursor = conn.cursor()

    # Load orbit data from the MPC
    mpc_obs = kete.mpc.fetch_known_orbit_data(force_download=True)

    # Convert that data to State objects.
    mpc_states = kete.mpc.table_to_states(mpc_obs)

    # Update the states based on the current time. 
    mpc_states_jd = kete.propagate_n_body(mpc_states, time_jd)

    # Update the states for the next day. This will be used to find the proper motion.
    mpc_states_jd_plus_1 = kete.propagate_n_body(mpc_states_jd, time_jd + 1)

    # Convert the states from ecliptic to equatorial to easily access the RA and Dec.
    for s_jd, s_jd_plus_1 in zip(mpc_states_jd, mpc_states_jd_plus_1):

        # Convert the state to equatorial coordinates.
        ra_jd = s_jd.as_equatorial.pos.ra 
        dec_jd = s_jd.as_equatorial.pos.dec 

        ra_jd_plus_1 = s_jd_plus_1.as_equatorial.pos.ra 
        dec_jd_plus_1 = s_jd_plus_1.as_equatorial.pos.dec 

        # Compute the average proper motion for the next day. In deg/day.
        pm_ra = ra_jd_plus_1 - ra_jd
        pm_dec = dec_jd_plus_1 - dec_jd
        
        # Add the state to the database. 
        # TODO: Add them to a temporary table first and then replace the main table.
        cursor.execute("INSERT INTO asteroids (designation, jd, ra, dec, pm_ra, pm_dec, position_x, position_y, position_z, velocity_x, velocity_y, velocity_z) VALUES \
                       (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)", \
                       (s_jd.desig, s_jd.jd, ra_jd, dec_jd, pm_ra, pm_dec, s_jd.pos.x, s_jd.pos.y, s_jd.pos.z, s_jd.vel.x, s_jd.vel.y, s_jd.vel.z))

    conn.commit()

    # Close the connection
    cursor.close()
    conn.close()

if __name__=="__main__":

    start_time = time.perf_counter()

    # Feb 4, 2025
    jd = kete.Time.from_ymd(2025, 2, 23).jd
    
    create_database_table()
    download_and_update_db(jd)

    
    end_time = time.perf_counter()

    print(f"Time taken: {end_time - start_time:.2f} seconds")