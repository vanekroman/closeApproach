from math import floor

import numpy as np

# For plotting
from astropy import time
from astropy import units as u
from numpy import arctan, cos, sin, sqrt, tan
from poliastro.bodies import Moon, Sun
from poliastro.constants import J2000, GM_earth, R_pluto
from poliastro.examples import churi, iss
from poliastro.plotting import OrbitPlotter
from poliastro.plotting.orbit.backends import Plotly3D
from poliastro.plotting.orbit.plotter import Trajectory
from poliastro.twobody import Orbit


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
        self.eccentricanomaly = np.zeros(30)
        self.pos = np.zeros(3)
        self.velocity = np.zeros(3)
        self.p = 0.0


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
    T = np.array(
        [
            [
                np.cos(omega) * np.cos(Omega) - np.sin(omega) * np.sin(Omega) * np.cos(i),
                -np.sin(omega) * np.cos(Omega) - np.cos(omega) * np.sin(Omega) * np.cos(i),
                np.sin(Omega) * np.sin(i),
            ],
            [
                np.cos(omega) * np.sin(Omega) + np.sin(omega) * np.cos(Omega) * np.cos(i),
                -np.sin(omega) * np.sin(Omega) + np.cos(omega) * np.cos(Omega) * np.cos(i),
                -np.cos(Omega) * np.sin(i),
            ],
            [
                np.sin(omega) * np.sin(i),
                np.cos(omega) * np.sin(i),
                np.cos(i),
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
    r = p / (1 + e * np.cos(theta))

    # Cartesian coordinates in the plane
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    z = 0  # Since it's a 2D orbit in the xy-plane

    return np.array([x, y, z])


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
    vx = up * (-np.sin(theta))
    vy = up * (e + np.cos(theta))
    vz = 0  # Since it's a 2D orbit in the xy-plane

    return np.array([vx, vy, vz])


# frame = OrbitPlotter(backend=Plotly3D(use_dark_theme=True))

# for body in bodies:
#     a = body.a / au << u.AU
#     ecc = body.ecc << u.one
#     inc = body.inc << u.rad
#     raan = body.raan << u.rad
#     argp = body.argp << u.rad
#     nu = body.trueanomaly << u.rad
#
#     jd = time.Time(ctime, format="jd")
#
#     orb = Orbit.from_classical(Sun, a, ecc, inc, raan, argp, nu, epoch=jd)
#     frame.plot(orb)
#
# frame.show()
