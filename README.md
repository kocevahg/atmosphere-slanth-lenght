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

## Slant Depth

The atmospheric slant depth is calculated as

$$
X = \int \rho(h(l))\,dl,
$$

where `ρ` is the atmospheric density and `l` is the distance along the particle trajectory.

The result is expressed in:

$$
\mathrm{g/cm^2}.
$$

For curved-Earth geometry, the altitude along the ray is calculated from

$$
r(l)=
\sqrt{
R_{\rm obs}^2+l^2+
2R_{\rm obs}l\cos\theta
},
$$

with

$$
h(l)=r(l)-R_{\rm Earth}.
$$

This becomes important for large zenith angles, where the flat-atmosphere approximation breaks down.

## Curved vs. Flat Geometry

Curved geometry is the default:

```python
ray = Ray(85.0, curved=True)
```

For comparison, a flat atmosphere can be used with:

```python
ray = Ray(85.0, curved=False)
```

The curved calculation is especially important for trajectories close to the horizon.

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
\cos\theta^*
=
\sqrt{
\frac{
x^2+p_1^2+p_2x^{p_3}+p_4x^{p_5}
}{
1+p_1^2+p_2+p_4
}
},
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

### Reference

D. Chirkin, **“Fluxes of Atmospheric Leptons at 600 GeV - 60 TeV”**, 2004.

[arXiv:hep-ph/0407078](https://arxiv.org/abs/hep-ph/0407078)

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
