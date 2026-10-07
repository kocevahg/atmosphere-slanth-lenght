# Atmospheric Slant Depth

Python code for calculating atmospheric properties and the slant depth traversed by a particle through the Earth's atmosphere.

The atmosphere is modeled using the **US Standard Atmosphere 1976** via the [`ussa1976`](https://pypi.org/project/ussa1976/) package.

The code supports both **curved-Earth** and **flat-Earth** geometries.

## Installation

Install the required packages with:

```bash
pip install numpy ussa1976
```

## Basic Usage

Calculate the total atmospheric slant depth for a particle arriving at a given zenith angle:

```python
from atmosphere import slant_depth

X = slant_depth(60.0)

print(X, "g/cm²")
```

For more information about the trajectory, use the `Ray` class:

```python
from atmosphere import Ray

ray = Ray(60.0)

print("Total depth:", ray.X_total, "g/cm²")
print("Path length:", ray.l[-1], "km")
print("Top altitude:", ray.h[-1], "km")
```

The main arrays are:

```python
ray.l      # distance along ray [km]
ray.h      # altitude [km]
ray.rho    # density [g/cm³]
ray.x      # cumulative slant depth [g/cm²]
```

## Atmospheric Model

The atmosphere extends from sea level to:

```python
H_TOP = 86  # km
```

The Earth is modeled with:

```python
R_EARTH = 6371  # km
```

Atmospheric quantities are obtained from the US Standard Atmosphere 1976.

```python
density(10.0)       # g/cm³
temperature(10.0)   # K
pressure(10.0)      # Pa
```

## Atmospheric Slant Depth

For a given zenith angle, the particle trajectory is represented by a **ray** originating at the detector and extending through the atmosphere. The variable `l` denotes the distance traveled along this ray, with

$$
l=0
$$

at the detector.

The atmospheric density changes along the trajectory because the particle moves through regions of different altitude. The altitude is therefore written as a function of the distance along the ray,

$$
h=h(l),
$$

and the corresponding atmospheric density is

$$
\rho(l)=\rho\left(h(l)\right).
$$

The **atmospheric slant depth** encountered by the particle after traveling a distance $l$ is

$$
X(l)=\int_0^l \rho\left(h(l')\right)\,dl',
$$

The **total atmospheric slant depth** is obtained by integrating from the detector to the point where the ray reaches the upper boundary of the atmospheric model:

$$
X_{\mathrm{total}}=\int_0^{L_{\mathrm{max}}}\rho\left(h(l)\right)\,dl,
$$

where $L_{\mathrm{max}}$ is the distance along the ray to the upper atmospheric boundary.

The slant depth is a **column density**, with units

$$
\mathrm{g/cm^2}.
$$

In the numerical calculation, the distance along the ray is stored in kilometres, while the atmospheric density is given in $\mathrm{g/cm^3}$. Therefore, the path length is converted from kilometres to centimetres using

$$
1~\mathrm{km}=10^5~\mathrm{cm}.
$$

### Numerical Integration

The integral is evaluated numerically using the **trapezoidal rule**. For two neighboring points along the ray, $l_i$ and $l_{i+1}$, the contribution to the slant depth is approximated by

$$
\Delta X_i\approx\frac{\rho_i+\rho_{i+1}}{2}\left(l_{i+1}-l_i\right),
$$

with the path length converted from kilometres to centimetres.

The cumulative slant depth is then obtained by summing these contributions:

$$
X(l_i)\approx\sum_{j=0}^{i-1}\frac{\rho_j+\rho_{j+1}}{2}\left(l_{j+1}-l_j\right).
$$

Thus, the `Ray` class stores both the geometric trajectory and the atmospheric quantities evaluated along it:

```python
ray.l      # distance along particle trajectory [km]
ray.h      # altitude along trajectory [km]
ray.rho    # atmospheric density along trajectory [g/cm³]
ray.x      # cumulative atmospheric slant depth [g/cm²]
```

The final value,

```python
ray.X_total
```

is the total atmospheric slant depth encountered by the particle between the detector and the upper atmospheric boundary.

### Curved-Earth Geometry

For curved-Earth geometry, the detector is located at

$$
R_{\mathrm{obs}}=R_{\mathrm{Earth}}+h_{\mathrm{obs}},
$$

where $h_{\mathrm{obs}}$ is the detector altitude.

For a ray with zenith angle $\theta$, the distance from the Earth's center after traveling a distance $l$ along the trajectory is

$$
r(l)=\sqrt{R_{\mathrm{obs}}^2+l^2+2R_{\mathrm{obs}}l\cos\theta}.
$$

The altitude above the Earth's surface is therefore

$$
h(l)=r(l)-R_{\mathrm{Earth}}.
$$

The atmospheric density along the trajectory is then evaluated as

$$
\rho(l)=\rho\left(h(l)\right).
$$

The maximum integration distance $L_{\mathrm{max}}$ is determined by the point where the ray reaches the upper atmospheric boundary,

$$
R_{\mathrm{top}}=R_{\mathrm{Earth}}+H_{\mathrm{top}},
$$

with

$$
H_{\mathrm{top}}=86~\mathrm{km}.
$$

The curved-Earth treatment is particularly important for large zenith angles. For nearly horizontal trajectories, the Earth's curvature significantly changes the altitude of the particle along its path, making the flat-atmosphere approximation increasingly inaccurate.

## Layer-by-Layer Cross-Check

An independent layer-based calculation is provided by:

```python
user_layer_method(theta_deg)
```

For example:

```python
X_ray = slant_depth(60.0)
X_layer = user_layer_method(60.0)

print(X_ray)
print(X_layer)
```

The two methods can be compared to check the numerical integration.

## Chirkin Effective Zenith Angle

The function

```python
cos_theta_star(theta_deg)
```

implements the effective zenith-angle parametrization introduced by Dmitry Chirkin.

It is given by

$$
\cos\theta^*=\sqrt{\frac{x^2+p_1^2+p_2x^{p_3}+p_4x^{p_5}}{1+p_1^2+p_2+p_4}},
$$

where

$$
x=\cos\theta.
$$

The parameters are

```text
p1 = 0.102573
p2 = -0.068287
p3 = 0.958633
p4 = 0.0407253
p5 = 0.817285
```

This parametrization accounts for the effect of atmospheric curvature on inclined atmospheric particle trajectories. 

## Chirkin Parametrization

As a reference for the numerical calculation, the code also includes the five-parameter parametrization of the total atmospheric overburden from Chirkin [1]:

$$
X_{\mathrm{tot}}(\theta)=\frac{1\ \mathrm{mwe}}{p_1+p_2\cos^{p_3}(\theta)+p_4\left(1-\cos^2(\theta)\right)^{p_5}}.
$$

For the `Xtot` fit, Chirkin gives

```text
p1 = -0.017326
p2 =  0.114236
p3 =  1.15043
p4 =  0.0200854
p5 =  1.16714
```

The implementation is available as

```python
slant_depth_chirkin(theta_deg)
```

and returns the total atmospheric overburden in **g/cm²**.

This parametrization provides an independent reference for the numerical curved-Earth integration:

```python
X_ray = slant_depth(theta)
X_chirkin = slant_depth_chirkin(theta)
```

The two results can therefore be compared as a function of zenith angle to assess the agreement between the numerical atmospheric integration and the Chirkin parametrization.

**Reference**

[1] D. Chirkin, *Fluxes of Atmospheric Leptons at 600 GeV - 60 TeV*, [arXiv:hep-ph/0407078](https://arxiv.org/abs/hep-ph/0407078). 


## Units

| Quantity           | Unit  |
| ------------------ | ----- |
| Altitude           | km    |
| Distance along ray | km    |
| Density            | g/cm³ |
| Temperature        | K     |
| Pressure           | Pa    |
| Slant depth        | g/cm² |

## Notes

This calculation assumes:

* US Standard Atmosphere 1976
* spherical Earth with radius 6371 km
* atmospheric boundary at 86 km
* straight-line particle trajectories
* no magnetic deflection or scattering
