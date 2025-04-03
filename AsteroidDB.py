import time
import kete
import psycopg2

from astropy.table import Table

# Database connection parameters
db_params = {
    "dbname": "asteroid",
    "user": "vedshah",  # Change to your PostgreSQL username
    "host": "localhost",  # Change if PostgreSQL is hosted remotely
    "port": "5432"  # Default PostgreSQL port
}

def create_database_table(db_params=db_params):
    """
    Create a database table to store the asteroid data.
    """
    try:
        # Establish connection
        conn = psycopg2.connect(**db_params)
        cursor = conn.cursor()
        
        # Create table if it doesn't exist
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS asteroids (
                designation TEXT PRIMARY KEY,
                jd FLOAT NOT NULL,
                ra FLOAT NOT NULL,
                dec FLOAT NOT NULL,
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

    # Update the states based on the time. 
    mpc_states = kete.propagate_n_body(mpc_states, time_jd)

    # Convert the states from ecliptic to equatorial to easily access the RA and Dec.
    for i, s in enumerate(mpc_states):

        # Convert the state to equatorial coordinates.
        ra = s.as_equatorial.pos.ra
        dec = s.as_equatorial.pos.dec

        # Add the state to the database.
        cursor.execute("INSERT INTO asteroids (designation, jd, ra, dec, position_x, position_y, position_z, velocity_x, velocity_y, velocity_z) VALUES \
                       (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)", (s.desig, s.jd, ra, dec, s.pos.x, s.pos.y, s.pos.z, s.vel.x, s.vel.y, s.vel.z))


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