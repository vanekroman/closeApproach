from math import radians

import numpy as np
from astropy import time
from astropy import units as u
from flask import Flask, render_template, request
from numpy import arctan, cos, sin, sqrt, tan
from poliastro.bodies import Sun
from poliastro.plotting import OrbitPlotter
from poliastro.plotting.orbit.backends import Plotly3D
from poliastro.twobody import Orbit

from heliopy import (
    Body,
    polar_to_cartesian_position,
    polar_to_cartesian_velocity,
    transformation_matrix,
)

# from typing import List


app = Flask(__name__)

# Gravitational parameters
uSun = 2.959122083e-4  # AU^3 day^-2
au = 149597870.700  # km

# Example predefined bodies
predefined_bodies = [
    Body("Earth", 1, 0.0167, 0, -11.26064, 114.20783, 2451547.5191),
    Body("LV2021", 1.3117, 0.4316, 16.4732, 246.4812, 276.6864, 2459305.2444),
]

# Main body list for rendering
bodies = []

frame = OrbitPlotter(backend=Plotly3D())


@app.route("/", methods=["GET", "POST"])
def index():
    global bodies
    if request.method == "POST":
        if "add_predefined" in request.form:
            # Add a predefined body
            selected_name = request.form.get("predefined_body")
            selected_body = next((b for b in predefined_bodies if b.name == selected_name), None)
            if selected_body:
                print(selected_body, flush=True)

        elif "add_custom" in request.form:
            # Add a custom-defined body
            name = request.form.get("name")
            a = float(request.form.get("a"))
            ecc = float(request.form.get("ecc"))
            inc = float(request.form.get("inc"))
            raan = float(request.form.get("raan"))
            argp = float(request.form.get("argp"))
            tpassp = float(request.form.get("tpassp"))
            bodies.append(Body(name, a, ecc, inc, raan, argp, tpassp))

        # Render orbits after adding bodies
        distances, velocities = calculate_orbits_and_plot(selected_body)
        return render_template("index.html", predefined_bodies=predefined_bodies, distances=distances, velocities=velocities)

    return render_template("index.html", predefined_bodies=predefined_bodies)


ctime = 2459362.0555


def calculate_orbits_and_plot(body):
    global ctime
    global bodies
    ## Main Calculations

    # Mean motion
    body.n = np.sqrt(uSun / (body.a * body.a * body.a))

    #  Unit conversion
    body.a = body.a * au
    body.inc = radians(body.inc)
    body.raan = radians(body.raan)
    body.argp = radians(body.argp)

    # Type of orbit
    if body.ecc < 1:
        # elliptic
        # Orbital parameter
        body.p = body.a * (1 - body.ecc * body.ecc)

        # Mean anomaly
        body.meananomaly = body.n * (ctime - body.tpassp)

        # Eccentric anomaly
        body.eccentricanomaly[0] = body.meananomaly
        # Ecaluate eccentric anomaly
        for i in range(0, body.eccentricanomaly.size - 1):
            body.eccentricanomaly[i + 1] = body.meananomaly + body.ecc * sin(body.eccentricanomaly[i])

        # True anomaly
        body.trueanomaly = 2 * arctan(sqrt((1 + body.ecc) / (1 - body.ecc)) * tan(body.eccentricanomaly[-1] / 2))

    elif body.ecc == 1:
        # parabolic
        body.p = 2 * body.a
    else:
        # hyperbolic
        body.p = body.a * (pow(body.ecc, 2) - 1)

    # 2D State position vector
    position_perifocal = polar_to_cartesian_position(body.p, body.ecc, body.trueanomaly)
    # 2D State velocity vector
    uSun_km = uSun * (au * au * au) / ((24 * 60 * 60) * (24 * 60 * 60))
    velocity_perifocal = polar_to_cartesian_velocity(body.p, body.ecc, body.trueanomaly, uSun_km)

    # Transformation matrix
    T = transformation_matrix(body.argp, body.raan, body.inc)

    # 3D State position vector
    body.pos = T @ position_perifocal
    # 3D State velocity vector
    body.velocity = T @ velocity_perifocal

    distances = {}
    velocities = {}
    # Calculate distances between each pair of bodies
    for i in range(len(bodies)):
        for j in range(i + 1, len(bodies)):
            body1 = bodies[i]
            body2 = bodies[j]

            # Calculate the Euclidean distance
            distance = np.linalg.norm(body1.pos - body2.pos)
            # Calculate the relative velocity (Euclidean norm of velocity difference)
            relative_velocity = np.linalg.norm(body1.velocity - body2.velocity)

            # Store distance and relative velocity in dictionaries
            distances[(body1.name, body2.name)] = distance
            velocities[(body1.name, body2.name)] = relative_velocity

    # Output distances and relative speed
    for (body1_name, body2_name), distance in distances.items():
        relative_velocity = velocities[(body1_name, body2_name)]
        print(f"Distance between {body1_name} and {body2_name}: {(distance / au):.2f} AU or {distance:.2f} km")
        print(f"Relative velocity between {body1_name} and {body2_name}: {(relative_velocity / au):.2f} AU or {relative_velocity:.2f} km/s")

    global frame

    # Generate plot for each body
    a = body.a / au << u.AU
    ecc = body.ecc << u.one
    inc = body.inc << u.rad
    raan = body.raan << u.rad
    argp = body.argp << u.rad
    nu = body.trueanomaly << u.rad

    jd = time.Time(ctime, format="jd")

    orb = Orbit.from_classical(Sun, a, ecc, inc, raan, argp, nu, epoch=jd)

    frame.plot(orb)

    bodies.append(body)
    # Retrieve the HTML div with Plotly figure
    frame.show()

    return distances, velocities


if __name__ == "__main__":
    app.run(debug=True)
