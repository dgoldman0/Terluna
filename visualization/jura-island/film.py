"""The Jura Island flyover's timeline: frame size and rate, the camera path, the chapter panels and the frame's texts.

Nothing here needs OpenGL, so the tests can check the film without rendering it. The camera's position follows
clamped cubic splines in longitude and latitude and in the logarithm of its altitude. Its heading and pitch follow
their own clamped cubic splines through the look directions at the keyframes, the heading unwrapped, so a target
passing under the camera cannot flip it. Five chapters each hold one fixed panel of two named reference points, and
a panel changes with a smoothstep fade out and in. Travel time is compressed; heights and distances are true.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.interpolate import CubicSpline

from shared.constants import MOON_RADIUS

WIDTH, HEIGHT, FPS, DURATION = 1600, 900, 60, 44
FIELD_OF_VIEW_DEG = 58.0                  # vertical
R_KM = MOON_RADIUS / 1000.0               # the reference sphere, on which the common water surface is drawn
LON0, LAT0 = -36.0, 45.0                  # origin of the local frame (degrees east, north)
MESH = dict(south=34.0, north=56.0, west=-53.0, east=-19.0)   # the regional mesh at 16 px/deg; a coarser one beyond

_la0, _lo0 = np.radians([LAT0, LON0])
BASIS = np.array([[-np.sin(_lo0), np.cos(_lo0), 0.0],
                  [-np.sin(_la0) * np.cos(_lo0), -np.sin(_la0) * np.sin(_lo0), np.cos(_la0)],
                  [np.cos(_la0) * np.cos(_lo0), np.cos(_la0) * np.sin(_lo0), np.sin(_la0)]])   # east, north, up


def xyz(lon, lat, h_km=0.0):
    """Positions (km) in the local frame, east, north and up from the origin on the reference sphere, of points at
    longitude and latitude (degrees) and height h_km above the sphere."""
    la, lo = np.radians(lat), np.radians(lon)
    r = R_KM + np.asarray(h_km)
    raw = np.stack(np.broadcast_arrays(r * np.cos(la) * np.cos(lo), r * np.cos(la) * np.sin(lo), r * np.sin(la)),
                   axis=-1)
    return (raw - R_KM * BASIS[2]) @ BASIS.T


def norm(v):
    return v / np.maximum(np.linalg.norm(v, axis=-1, keepdims=True), 1e-9)


def local_axes(lon, lat):
    """East, north and up at a point, as vectors in the local frame."""
    lo, la = np.radians([lon, lat])
    east = np.array([-np.sin(lo), np.cos(lo), 0.0]) @ BASIS.T
    north = np.array([-np.sin(la) * np.cos(lo), -np.sin(la) * np.sin(lo), np.cos(la)]) @ BASIS.T
    up = np.array([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)]) @ BASIS.T
    return east, north, up


def look_at(eye, target, up):
    f = norm(target - eye)
    s = norm(np.cross(f, up))
    u = np.cross(s, f)
    m = np.eye(4, dtype=np.float32)
    m[0, :3], m[1, :3], m[2, :3] = s, u, -f
    m[:3, 3] = -m[:3, :3] @ eye
    return m


def perspective(fov_deg, aspect, near, far):
    t = 1 / math.tan(math.radians(fov_deg) / 2)
    return np.array([[t / aspect, 0, 0, 0], [0, t, 0, 0],
                     [0, 0, (far + near) / (near - far), 2 * far * near / (near - far)], [0, 0, -1, 0]], np.float32)


# Time (s); the camera's longitude and latitude (degrees) and altitude above the water surface (km); and the
# longitude, latitude and height (km) of the point it looks towards at that key.
KEYS = np.array([
    [0, -28.0, 38.6, 250.0, -36.5, 45.0, 0.0],
    [9, -28.7, 42.1, 62.0, -33.5, 46.1, 0.6],
    [18, -32.8, 46.2, 22.0, -35.4, 47.2, 0.8],
    [25, -38.6, 48.0, 18.0, -40.1, 46.7, 0.8],
    [32, -43.1, 44.1, 27.0, -41.6, 42.2, 0.6],
    [40, -35.0, 40.2, 425.0, -36.5, 44.8, 0.0],
    [44, -34.5, 41.5, 480.0, -36.5, 44.8, 0.0],
])
_POSITION = CubicSpline(KEYS[:, 0], KEYS[:, 1:], axis=0, bc_type='clamped')
_ALTITUDE = CubicSpline(KEYS[:, 0], np.log(KEYS[:, 3]), bc_type='clamped')


def path(t):
    """The camera's longitude, latitude and altitude (km), and the interpolated target's longitude, latitude and
    height, at time t (s, a scalar or an array)."""
    value = _POSITION(t)
    value[..., 2] = np.exp(_ALTITUDE(t))
    return value


def _key_angles():
    angles = []
    for _, lon, lat, alt, target_lon, target_lat, target_h in KEYS:
        east, north, up = local_axes(lon, lat)
        d = norm(xyz(target_lon, target_lat, target_h) - xyz(lon, lat, alt))
        angles.append([np.arctan2(d @ east, d @ north), np.arcsin(d @ up)])
    angles = np.array(angles)
    angles[:, 0] = np.unwrap(angles[:, 0])
    return angles


KEY_ANGLES = _key_angles()   # heading (clockwise from north) and pitch (radians) at each key, the heading unwrapped
_ORIENTATION = CubicSpline(KEYS[:, 0], KEY_ANGLES, axis=0, bc_type='clamped')


def orientation(t):
    """Heading and pitch (radians) at time t."""
    return _ORIENTATION(t)


def view_direction(t):
    """The unit vector the camera looks along at time t, in the local frame."""
    lon, lat = path(t)[:2]
    east, north, up = local_axes(lon, lat)
    heading, pitch = orientation(t)
    return np.cos(pitch) * (np.sin(heading) * east + np.cos(heading) * north) + np.sin(pitch) * up


def camera_pose(t):
    """The eye's position (km, local frame) and the view matrix at time t."""
    lon, lat, alt = path(t)[:3]
    eye = xyz(lon, lat, alt)
    up = local_axes(lon, lat)[2]
    return eye, look_at(eye, eye + view_direction(t), up)


def frame_times():
    return np.arange(FPS * DURATION) / FPS


# Each chapter's start and end (s), its title and its two named reference points, shown in a fixed panel.
CHAPTERS = (
    (0, 9, '01  Across the Sea of Rains', ('Mare Imbrium', 'Bay mouth')),
    (9, 19, '02  Into the Bay of Rainbows', ('Inner bay', 'Bay floor')),
    (19, 27, '03  Along the Jura mountain front', ('Jura summit', 'Northern water')),
    (27, 35, '04  Over the western uplands', ('Western water', 'Flooded crater')),
    (35, 44, '05  The island and its open bay', ('Jura summit', 'Bay floor')),
)
FADE_SECONDS = 0.75


def chapter_state(t):
    """The chapter showing at time t and its panel's opacity: the panel fades out over the last FADE_SECONDS of its
    chapter and the next fades in over the first, each by smoothstep."""
    index = min(sum(t >= c[0] for c in CHAPTERS) - 1, len(CHAPTERS) - 1)
    start, end = CHAPTERS[index][:2]
    opacity = 1.0
    if index > 0:
        opacity = min(opacity, np.clip((t - start) / FADE_SECONDS, 0, 1))
    if index < len(CHAPTERS) - 1:
        opacity = min(opacity, np.clip((end - t) / FADE_SECONDS, 0, 1))
    return index, float(opacity * opacity * (3 - 2 * opacity))


def camera_frame(t):
    """The camera's orientation at time t as a rotation matrix (rows: right, up, backwards), in double precision."""
    f = view_direction(t)
    s = norm(np.cross(f, local_axes(*path(t)[:2])[2]))
    return np.array([s, np.cross(s, f), -f])


def camera_checks(times=None):
    """What bounds the camera over the frame times: its lowest altitude, its pitch range, and how far the whole camera
    and its heading turn from one frame to the next."""
    times = frame_times() if times is None else np.asarray(times)
    values = path(times)
    heading, pitch = orientation(times).T
    frames = np.array([camera_frame(t) for t in times])
    # rotation angle between successive frames from the chord, |R1 - R2| = 2 sqrt(2) sin(angle / 2)
    chord = np.linalg.norm(frames[1:] - frames[:-1], axis=(1, 2))
    turn = np.degrees(2 * np.arcsin(np.clip(chord / (2 * np.sqrt(2)), 0.0, 1.0)))
    dt = float(np.min(np.diff(times)))
    return dict(frames=int(times.size), finite=bool(np.isfinite(values).all() and np.isfinite(frames).all()),
                minimum_altitude_km=float(values[:, 2].min()),
                pitch_range_deg=[float(np.degrees(pitch.min())), float(np.degrees(pitch.max()))],
                maximum_rotation_deg_per_frame=float(turn.max()), maximum_rotation_deg_per_s=float(turn.max() / dt),
                maximum_heading_step_deg=float(np.degrees(np.abs(np.diff(heading)).max())))


def captions(sea_level_m, water_share, pixels_per_degree):
    """The frame's fixed texts, filled from the island study's product and the atlas it was measured against."""
    level = f'{sea_level_m:,.1f}'.replace('-', '−')
    return dict(share=f'{water_share:.0%} of the Moon under water',
                level=f'Common water level {level} m relative to the geoid',
                datum='Heights and depths from the water surface',
                view='Common-water-level view  •  illustrative lighting and haze',
                terrain=f'LOLA terrain at {pixels_per_degree} px/deg  •  true relief  •  depth bands '
                        'overlaid  •  camera travel compressed')
