import numpy as np
import ussa1976

R_EARTH = 6371 # [km]
H_TOP = 86     # [km]

############### ATMOSPHERE ##########################

z = np.arange(0, H_TOP * 1000.0 + 1, 100) # [m]
atm = ussa1976.compute(z=z, variables=["t", "p", "rho"])

T=atm["t"].values     # [K]
P=atm["p"].values     # [Pa]
rho=atm["rho"].values # [kg/m³]

def density(z):
    z_arr = np.atleast_1d(np.asarray(z, dtype=float)) * 1000.0

    atm = ussa1976.compute(z=z_arr, variables=["rho"]) # [kg/m³]
    rho = atm["rho"].values * 1e-3  # kg/m³ -> g/cm³

    return rho[0] if np.ndim(z) == 0 else rho


def temperature(z):
    z_arr = np.atleast_1d(np.asarray(z, dtype=float)) * 1000.0

    atm = ussa1976.compute(z=z_arr, variables=["t"])
    T = atm["t"].values  # [K]

    return T[0] if np.ndim(z) == 0 else T


def pressure(z):
    z_arr = np.atleast_1d(np.asarray(z, dtype=float)) * 1000.0

    atm = ussa1976.compute(z=z_arr, variables=["p"])
    P = atm["p"].values  # [Pa]

    return P[0] if np.ndim(z) == 0 else P


############### GEOMETRY ##########################
class Ray:
    def __init__(
        self,
        theta_deg,
        h_obs=0.0,
        h_prod=15.0,
        curved=True,
        n=40000
    ):
        theta = np.radians(theta_deg)
        cos_theta = np.cos(theta)

        # radius at detector
        R_obs = R_EARTH + h_obs

        # radius at production point
        R_prod = R_EARTH + h_prod


        if curved:

            # check whether the ray reaches the requested altitude
            discriminant = ((R_obs * cos_theta)**2+ R_prod**2 - R_obs**2)
            Lmax = (-R_obs * cos_theta+ np.sqrt(discriminant))
        else:
            Lmax = (h_prod - h_obs) / max(cos_theta, 1e-12)
            Lmax = min(Lmax, 1e7)

        #integration grid
        l = np.concatenate([[0.0],np.geomspace(1e-4, Lmax, n)])

        if curved:
            r = np.sqrt( R_obs**2 + l**2 + 2.0 * R_obs * l * cos_theta)
            h = ((l**2 + 2.0 * R_obs * l * cos_theta) / (r + R_obs)+ h_obs )

        else:
            h = h_obs + l * cos_theta

        rho = density(h)
        dl = np.diff(l)                 # [km]
        
        #integrate rho*dl
        rho_avg = 0.5 * (rho[1:] + rho[:-1])   # [g/cm^3]
        dx = (rho_avg* dl* 1e5)         # [g/cm^2]

        #cumulative depth from detector
        x = np.concatenate([ [0.0],np.cumsum(dx)])

        self.l = l
        self.h = h
        self.rho = rho
        self.x = x

        self.X_total = x[-1]

        self.theta = theta_deg
        self.h_obs = h_obs
        self.h_prod = h_prod

    def l_of_x(self, x):
            #convert slant depth [g/cm^2] to distance along ray [km]
            return np.interp(x, self.x, self.l)

    def h_of_x(self, x):
            #convert slant depth [g/cm^2] to altitude [km]
            return np.interp(x, self.x, self.h)

############### TOTAL SLANT DEPTH [g/cm^2] ##########################
def slant_depth(theta_deg, h_obs=0.0, curved=True):
    ray = Ray(theta_deg=theta_deg, h_obs=h_obs, curved=curved)
    return ray.X_total

############### CHIRKIN EFFECTIVE ZENITH ANGLE ##########################
def cos_theta_star(theta_deg):

    x = np.cos(np.radians(theta_deg))
    p1 = 0.102573
    p2 = -0.068287
    p3 = 0.958633
    p4 = 0.0407253
    p5 = 0.817285
    return np.sqrt( ( x ** 2 + p1 ** 2 + p2 * x ** p3 + p4 * x ** p5 ) / ( 1 + p1 ** 2 + p2 + p4 ) )

############### CHIRKIN SLANT DEPTH FOR 86km ATMOSPHERE ##########################
def slant_depth_chirkin(theta_deg):
    x = np.cos(np.radians(theta_deg))

    p1 = -0.017326
    p2 = 0.114236
    p3 = 1.15043
    p4 = 0.0200854
    p5 = 1.16714

    X_mwe = 1.0 / ( p1 + p2 * x**p3 + p4 * (1.0 - x**2)**p5 )  #[mwe]

    return X_mwe * 100.0 # 1 mwe = 100 g/cm^2
############### LAYER BY LAYER METHOD ##########################
def user_layer_method( theta_deg, zmax=15):
    # Atmospheric layers
    z = np.arange(0.0, zmax + 1.0, 0.1)  # [km]
    rho = density(z)  # [g/cm^3]

    theta = np.radians(theta_deg)
    cos_theta = np.cos(theta)

    # Detector radius
    R_obs = R_EARTH

    total = 0.0

    for i in range(len(z) - 1):
        h1 = z[i]
        h2 = z[i + 1]

        # radii of bottom and top of the layer
        r1 = R_EARTH + h1
        r2 = R_EARTH + h2

        # distance along the ray at which it reaches radius r
        l1 = (-R_obs * cos_theta+ np.sqrt(r1 ** 2 - R_obs ** 2 * np.sin(theta) ** 2))
        l2 = (-R_obs * cos_theta + np.sqrt( r2 ** 2 - R_obs ** 2 * np.sin(theta) ** 2))

        # path length through this layer
        dl = l2 - l1  # [km]

        # average density in the layer
        rho_avg = 0.5 * (rho[i] + rho[i + 1])

        # contribution to slant depth
        total += rho_avg * dl * 1e5

    return total

############### TEST CODE ##########################
if __name__ == "__main__":
# Atmospheric parameters at 10 km
    print("Atmosphere at 10 km:")
    print( "rho =", density(10.0), "g/cm^3" )
    print( "T =", temperature(10.0), "K" )
    print( "P =", pressure(10.0), "Pa" )
    # Slant depth print()
    print("Slant depth:")
    for theta in [0, 30, 60, 85,89, 90]:
        X = slant_depth(theta)
        X_chirkin = slant_depth_chirkin(theta)
        print("Ray method:")
        print( f"theta = {theta:2d} deg : " f"X = {X:.3f} g/cm^2" )
        print("Layer method:")
        print(f"theta = {theta:2d} deg : " f"X = {user_layer_method(theta)} g/cm^2")
        print("Chirkin method:")
        print(f"theta = {theta:2d} deg : " f"X = {X_chirkin} g/cm^2")

    # Full ray
    ray = Ray(90.0)
    print()
    print("Ray at 60 degrees:")
    print("X_total =", ray.X_total, "g/cm^2")
    print("Lmax =", ray.l[-1], "km")
    print("h_top =", ray.h[-1], "km")
