"""The depth scale the map and the flyover share, and the named points whose depths they show.

Bands of 250 m to 1,500 m keep the bay's detail; the last band holds everything from 2,000 m down.
"""
import numpy as np

DEPTH_LEVELS_M = np.array([0, 250, 500, 750, 1000, 1250, 1500, 2000, 3000], float)
DEPTH_COLOURS = ['#d7eee4', '#a4dace', '#74bfbf', '#4a9eaf', '#337d9b', '#285b83', '#25416b', '#172b4b']
DEPTH_RGB = np.array([[int(c[i:i + 2], 16) / 255 for i in (1, 3, 5)] for c in DEPTH_COLOURS])
CONTOUR_INTERVAL_M = 250.0
# Named reference points: longitude (degrees east, negative west) and latitude
SOUNDINGS = (
    ('Inner bay', -32.0, 46.5),
    ('Bay floor', -31.5, 45.0),
    ('Bay mouth', -28.0, 42.0),
    ('Mare Imbrium', -25.0, 38.0),
    ('Northern water', -35.0, 53.0),
    ('Western water', -48.0, 44.0),
    ('Flooded crater', -43.5, 41.6),
)


def colours(depth_m):
    """The band colour (RGB, 0-1) of each depth in metres."""
    i = np.searchsorted(DEPTH_LEVELS_M[1:], np.asarray(depth_m), side='right')
    return DEPTH_RGB[np.minimum(i, len(DEPTH_RGB) - 1)]
