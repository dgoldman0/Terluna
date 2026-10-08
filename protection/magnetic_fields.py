"""Finite circular current paths in vacuum; no plasma or material response.

Elliptic-integral field: Simpson et al., NASA 20010038494, eqs. 24--32.
Current means ampere-turns. Coordinates and all returned quantities are SI.
"""
import numpy as np
from scipy.special import ellipk, ellipe
from shared import constants as K


def basis(normal):
    n = np.asarray(normal, dtype=float)
    n = n / np.linalg.norm(n)
    axis = np.eye(3)[np.argmin(abs(n))]
    x = np.cross(axis, n); x /= np.linalg.norm(x)
    return np.column_stack([x, np.cross(n, x), n])


def loop_field(points, radius, current=1., centre=(0., 0., 0.), normal=(0., 0., 1.)):
    """Exact filament field off the wire, with a near-axis expansion."""
    if radius <= 0:
        raise ValueError('Positive loop radius required')
    original = np.asarray(points, dtype=float)
    frame = basis(normal)
    p = (original.reshape(-1, 3)-centre) @ frame
    rho = np.linalg.norm(p[:, :2], axis=1); z = p[:, 2]
    den = (radius-rho)**2+z*z
    if np.any(den < radius**2*1e-24):
        raise ValueError('Field evaluated on filament')
    b = np.zeros_like(p)
    axial = rho < 1e-6*np.sqrt(radius**2+z*z)
    zz = z[axial]; rr = radius**2+zz**2
    b[axial, 2] = K.VACUUM_PERMEABILITY*current*radius**2/(2*rr**1.5)
    b[axial, :2] = (3*K.VACUUM_PERMEABILITY*current*radius**2*zz/(4*rr**2.5))[:, None]*p[axial, :2]
    r = rho[~axial]; zz = z[~axial]; d = den[~axial]
    root = np.sqrt((radius+r)**2+zz*zz)
    parameter = 4*radius*r/root**2
    ek, ee = ellipk(parameter), ellipe(parameter)
    prefactor = K.VACUUM_PERMEABILITY*current/(2*np.pi*root)
    br = prefactor*zz/r*(-ek+(radius**2+r*r+zz*zz)/d*ee)
    b[~axial, :2] = br[:, None]*p[~axial, :2]/r[:, None]
    b[~axial, 2] = prefactor*(ek+(radius**2-r*r-zz*zz)/d*ee)
    return (b @ frame.T).reshape(original.shape)


def field(points, loops):
    return sum((loop_field(points, **loop) for loop in loops), np.zeros_like(np.asarray(points, dtype=float)))


def loop_samples(loop, count=512):
    phi = np.arange(count)*2*np.pi/count
    frame = basis(loop['normal']); a = loop['radius']
    offset = np.c_[a*np.cos(phi), a*np.sin(phi), np.zeros(count)] @ frame.T
    dl = np.c_[-a*np.sin(phi), a*np.cos(phi), np.zeros(count)] @ frame.T*(2*np.pi/count)
    return offset+loop['centre'], dl


def mutual_load(source, target, count=512):
    p, dl = loop_samples(target, count)
    df = target['current']*np.cross(dl, loop_field(p, **source))
    return dict(force_N=df.sum(axis=0),
                torque_N_m=np.cross(p-target['centre'], df).sum(axis=0),
                maximum_line_load_N_m=float(np.max(np.linalg.norm(df, axis=1)/np.linalg.norm(dl, axis=1))))


def mutual_inductance(source, target, count=256):
    """Neumann integral of two distinct filaments; signed by their normals."""
    p, dp = loop_samples(source, count)
    q, dq = loop_samples(target, count)
    distance = np.linalg.norm(p[:, None, :]-q[None, :, :], axis=2)
    if distance.min() == 0:
        raise ValueError('Intersecting current paths')
    return float(K.VACUUM_PERMEABILITY/(4*np.pi)*np.sum((dp @ dq.T)/distance))


def sphere(radius, latitudes=61, longitudes=120):
    th = np.linspace(0, np.pi, latitudes)
    ph = np.arange(longitudes)*2*np.pi/longitudes
    t, p = np.meshgrid(th, ph, indexing='ij')
    return radius*np.c_[np.sin(t).ravel()*np.cos(p).ravel(),
                       np.sin(t).ravel()*np.sin(p).ravel(), np.cos(t).ravel()]


def regional_loops(radius, moment, count=4, historical=False):
    """Four disjoint spherical caps, or two polar caps. No planetary ring.

    New paths lie on the reference sphere; relief is not represented.
    Historical geometry retains its parallel displaced planar loops.
    """
    if not 0 < radius < K.MOON_RADIUS or count not in (2, 4):
        raise ValueError('Regional radius and two/four circuits required')
    loops = []
    if historical:
        if count != 4:
            raise ValueError('Historical source has four circuits')
        h = np.sqrt(K.MOON_RADIUS**2-(200e3)**2)
        for x, z in [(200e3, h), (-200e3, h), (200e3, -h), (-200e3, -h)]:
            loops.append(dict(radius=radius, centre=np.array([x, 0., z]),
                              normal=np.array([0., 0., 1.]), current=moment/(4*np.pi*radius**2)))
    else:
        beta = np.arcsin(radius/K.MOON_RADIUS)
        alpha = beta+np.deg2rad(5.) if count == 4 else 0.
        h = np.sqrt(K.MOON_RADIUS**2-radius**2)
        for pole, phi in ([(1, 0), (1, np.pi), (-1, np.pi/2), (-1, 3*np.pi/2)] if count == 4 else [(1, 0), (-1, 0)]):
            radial = np.array([np.sin(alpha)*np.cos(phi), np.sin(alpha)*np.sin(phi), pole*np.cos(alpha)])
            loops.append(dict(radius=radius, centre=h*radial, normal=pole*radial,
                              current=moment/(count*np.cos(alpha)*np.pi*radius**2)))
        axes = np.array([l['centre']/np.linalg.norm(l['centre']) for l in loops])
        nearest = np.arccos(np.clip((axes @ axes.T)[np.triu_indices(count, 1)], -1, 1)).min()
        if nearest <= 2*beta:
            raise ValueError('Regional caps overlap')
    return loops
