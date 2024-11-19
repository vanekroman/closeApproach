from math import cosh, floor, radians, sinh, sqrt, tanh

from numpy import arctan, array, cos, linalg, log, sin, tan, zeros

# Gravitational parameters
uSun = 2.959122083e-4  # AU^3 day^-2
au = 149597870.700  # km


class Body:
    def __init__(
        self,
        name: str,
        a: float,
        ecc: float,
        inc: float,
        raan: float,
        argp: float,
        tpassp: float,
    ):
        # Instance variables
        self.name = name
        self.a = a
        self.ecc = ecc
        self.inc = inc
        self.raan = raan
        self.argp = argp
        self.tpassp = tpassp

        # These should be instance variables, so each object has its own version
        self.n = 0.0
        self.trueanomaly = 0.0
        self.meananomaly = 0.0
        self.eccentricanomaly = zeros(30)
        self.pos = zeros(3)
        self.velocity = zeros(3)
        self.p = 0.0
        self.calculated = False


def gregorian_to_julian(year, month, day, ut=0.0):
    """
    Convert a Gregorian date to Julian Date (JD).
    year: int, Gregorian year
    month: int, Gregorian month
    day: int, Gregorian day
    ut: float, Universal Time in hours (default is 0 for midnight)

    Returns the Julian Date.

    # Example usage:
    gregorian_date = (2021, 10, 30, 0)  # Jan 15, 2000, noon UT
    JD = gregorian_to_julian(*gregorian_date)
    print("Julian Date:", JD)
    """
    # Adjust month and year if month is January or February
    JD0 = 367 * year - floor(7 * (year + floor((month + 9) / 12)) / 4) + floor(275 * month / 9) + day + 1721013.5

    # Add the fractional day based on UT
    JD = JD0 + ut / 24.0
    return JD


def julian_to_gregorian(JD):
    """
    Convert a Julian Date (JD) to a Gregorian date.
    JD: float, Julian Date (e.g., 2451545.0 for J2000)

    Returns (year, month, day, ut), where:
    year: int, Gregorian year
    month: int, Gregorian month
    day: int, Gregorian day
    ut: float, Universal Time in hours

    # Example usage:
    gregorian_converted = julian_to_gregorian(JD)
    print("Gregorian Date:", gregorian_converted)
    """
    JD += 0.5
    Z = floor(JD)
    F = JD - Z

    if Z >= 2299161:
        alpha = floor((Z - 1867216.25) / 36524.25)
        A = Z + 1 + alpha - floor(alpha / 4)
    else:
        A = Z

    B = A + 1524
    C = floor((B - 122.1) / 365.25)
    D = floor(365.25 * C)
    E = floor((B - D) / 30.6001)

    day = B - D - floor(30.6001 * E) + F
    month = E - 1 if E < 14 else E - 13
    year = C - 4716 if month > 2 else C - 4715

    # Extract the integer day and the fractional part for UT
    day_int = floor(day)
    ut = (day - day_int) * 24.0  # convert fraction of day to hours

    return int(year), int(month), int(day_int), ut


def transformation_matrix(omega, Omega, i):
    """
    Calculates the transformation matrix [T_{P?GE}].

    Parameters:
        omega (float): Argument of periapsis in radians.
        Omega (float): Longitude of the ascending node in radians.
        i (float): Inclination in radians.

    Returns:
        numpy.ndarray: The 3x3 transformation matrix.
    """
    # Define the transformation matrix elements
    T = array(
        [
            [
                cos(omega) * cos(Omega) - sin(omega) * sin(Omega) * cos(i),
                -sin(omega) * cos(Omega) - cos(omega) * sin(Omega) * cos(i),
                sin(Omega) * sin(i),
            ],
            [
                cos(omega) * sin(Omega) + sin(omega) * cos(Omega) * cos(i),
                -sin(omega) * sin(Omega) + cos(omega) * cos(Omega) * cos(i),
                -cos(Omega) * sin(i),
            ],
            [
                sin(omega) * sin(i),
                cos(omega) * sin(i),
                cos(i),
            ],
        ]
    )

    return T


def polar_to_cartesian_position(p, e, theta):
    """
    Converts polar coordinates to Cartesian coordinates for a 2D orbit.

    Parameters:
        p (float): Semi-latus rectum of the orbit.
        e (float): Eccentricity of the orbit.
        theta (float): True anomaly in radians.

    Returns:
        numpy.ndarray: The 3D Cartesian coordinates [x, y, z].
    """
    # Calculate radial distance r
    r = p / (1 + e * cos(theta))

    # Cartesian coordinates in the plane
    x = r * cos(theta)
    y = r * sin(theta)
    z = 0  # Since it's a 2D orbit in the xy-plane

    return array([x, y, z])


def polar_to_cartesian_velocity(p, e, theta, mu):
    """
    Converts the velocity vector from polar coordinates to Cartesian coordinates for a 2D orbit.

    Parameters:
        mu (float): Standard gravitational parameter in km^3 s-2.
        theta (float): True anomaly in radians.
        p (float): Semi-latus rectum of the orbit.
        e (float): Eccentricity of the orbit.
    Returns:
        numpy.ndarray: The 3D velocity vector [Vx, Vy, 0].
    """
    up = sqrt(mu / p)

    # Calculate the velocity components
    vx = up * (-sin(theta))
    vy = up * (e + cos(theta))
    vz = 0  # Since it's a 2D orbit in the xy-plane

    return array([vx, vy, vz])


def calculate_orbits(bodies, time):
    ## Main Calculations

    for body in bodies:
        if body.calculated:
            continue

        # Mean motion
        body.n = sqrt(uSun / (body.a * body.a * body.a))

        #  Unit conversion
        body.a = body.a * au
        body.inc = radians(body.inc)
        body.raan = radians(body.raan)
        body.argp = radians(body.argp)

        # Type of orbit
        if body.ecc < 1:
            # ELLIPSE
            # Orbital parameter
            body.p = body.a * (1 - body.ecc * body.ecc)

            # Mean anomaly
            body.meananomaly = body.n * (time - body.tpassp)

            # Eccentric anomaly
            body.eccentricanomaly[0] = body.meananomaly
            # Ecaluate eccentric anomaly
            for i in range(0, body.eccentricanomaly.size - 1):
                body.eccentricanomaly[i + 1] = body.meananomaly + body.ecc * sin(body.eccentricanomaly[i])

            # True anomaly
            body.trueanomaly = 2 * arctan(sqrt((1 + body.ecc) / (1 - body.ecc)) * tan(body.eccentricanomaly[-1] / 2))

        elif body.ecc == 1:
            # PARABOLA
            # Orbital parameter
            body.p = 2 * body.a

            # Mean anomaly
            body.meananomaly = body.n * (time - body.tpassp)

            g = 3 * sqrt(2) * body.meananomaly / 4
            y = (g + sqrt(g * g + 1)) ** (1 / 3)

            # True anomaly
            body.trueanomaly = 2 * arctan(y - 1 / y)

        else:
            # HYPERBOLA
            # Orbital paramter
            body.p = body.a * (body.ecc * body.ecc - 1)

            # Mean anomaly
            body.meananomaly = body.n * (time - body.tpassp)

            # Eccentric anomaly
            body.eccentricanomaly[0] = log(1.8 + 2 * (body.meananomaly / body.ecc))
            # Ecaluate eccentric anomaly
            for i in range(0, body.eccentricanomaly.size - 1):
                body.eccentricanomaly[i + 1] = body.eccentricanomaly[i] - ((body.ecc * sinh(body.eccentricanomaly[i]) - body.eccentricanomaly[i] - body.meananomaly) / (body.ecc * cosh(body.eccentricanomaly[i]) - 1))

            # True anomaly
            body.trueanomaly = 2 * arctan(sqrt((body.ecc + 1) / (body.ecc - 1)) * tanh(body.eccentricanomaly[-1] / 2))

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
            distance = linalg.norm(body1.pos - body2.pos)
            # Calculate the relative velocity (Euclidean norm of velocity difference)
            relative_velocity = linalg.norm(body1.velocity - body2.velocity)

            # Store distance and relative velocity in dictionaries
            distances[(body1.name, body2.name)] = distance
            velocities[(body1.name, body2.name)] = relative_velocity

    return distances, velocities
