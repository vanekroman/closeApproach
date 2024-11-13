import json

from astropy import units as u
from astropy.time import Time
from flask import Flask, render_template, request
from plotly.graph_objects import Figure
from plotly.utils import PlotlyJSONEncoder
from poliastro.bodies import Sun
from poliastro.plotting import OrbitPlotter
from poliastro.plotting.orbit.backends import Plotly3D
from poliastro.twobody import Orbit

# My lib
from heliopy import Body, au, calculate_orbits, gregorian_to_julian

app = Flask(__name__)


ctime = 2459362.0555

# Predefined bodies
predefined_bodies = [
    Body("Earth", 1, 0.0167, 0, -11.26064, 114.20783, 2451547.5191),
    Body("LV2021", 1.3117, 0.4316, 16.4732, 246.4812, 276.6864, 2459305.2444),
    Body("Oumuamua", 1.272345, 1.20113, 122.7417, 24.59691, 241.81054, 2458006.007321),
]

# Main body list for rendering
bodies = []

fig = Figure()
frame = OrbitPlotter(backend=Plotly3D(figure=fig, use_dark_theme=True))


@app.route("/", methods=["GET", "POST"])
def index():
    global bodies, ctime

    if request.method == "POST":
        # Determine if adding predefined or custom body
        if request.form.get("body_type") == "predefined":
            selected_name = request.form.get("predefined_body")
            selected_body = next((b for b in predefined_bodies if b.name == selected_name), None)
            if selected_body:
                existing_body = next((b for b in bodies if b.name == selected_name), None)
                if not existing_body:
                    bodies.append(selected_body)

        elif request.form.get("body_type") == "custom":
            # Collect custom body parameters
            name = request.form.get("name")
            a = float(request.form.get("a"))
            ecc = float(request.form.get("ecc"))
            inc = float(request.form.get("inc"))
            raan = float(request.form.get("raan"))
            argp = float(request.form.get("argp"))
            tpassp = float(request.form.get("tpassp"))
            existing_body = next((b for b in bodies if b.name == name), None)
            if existing_body:
                existing_body.a = a
                existing_body.ecc = ecc
                existing_body.inc = inc
                existing_body.raan = raan
                existing_body.argp = argp
                existing_body.tpassp = tpassp
            else:
                new_body = Body(name, a, ecc, inc, raan, argp, tpassp)
                bodies.append(new_body)

        # Handle time input
        time_format = request.form.get("time_format")
        if time_format == "gregorian":
            year = int(request.form.get("year"))
            month = int(request.form.get("month"))
            day = int(request.form.get("day"))
            hour = float(request.form.get("hour", 0))
            ctime = gregorian_to_julian(year, month, day, hour)
        elif time_format == "julian":
            ctime = float(request.form.get("julian_date"))

        # Calculate orbits with updated bodies and time
        distances, velocities = plot_bodies(bodies, frame, time=ctime)
        plot_data = json.dumps(fig, cls=PlotlyJSONEncoder)

        return render_template("index.html", predefined_bodies=predefined_bodies, distances=distances, velocities=velocities, plot_data=plot_data, current_time=ctime)

    return render_template("index.html", predefined_bodies=predefined_bodies)


# Just for creation of the plot
def plot_bodies(bodies, frame, time):
    # Main calcution calculated here
    distances, velocities = calculate_orbits(bodies, time)

    for body in bodies:
        if body.calculated:
            continue

        ecc = body.ecc << u.one
        inc = body.inc << u.rad
        raan = body.raan << u.rad
        argp = body.argp << u.rad
        nu = body.trueanomaly << u.rad
        jd = Time(time, format="jd")

        if body.ecc > 1:
            # Hyperbolic orbit
            a = body.a / au * (-1) << u.AU

            # nu_limit = hyp_nu_limit(ecc)
            # ephem = orb1.to_ephem(strategy=TrueAnomalyBounds())
            # orb = Orbit.from_ephem(Sun, ephem, epoch=jd)

        else:
            a = body.a / au << u.AU

        orb = Orbit.from_classical(Sun, a, ecc, inc, raan, argp, nu, epoch=jd)
        frame.plot(orb)
        body.calculated = True

    return distances, velocities


if __name__ == "__main__":
    # app.run(debug=True)
    import os

    port = int(os.environ.get("PORT", 5000))
    app.run(host="0.0.0.0", port=port)
