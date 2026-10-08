"""Ring 0 of the summit port: its design as explicit data (artistic: a design proposal, drawn by hand; the structure it
sits on is informed by research/studies/summit_tower/form.py).

Ring 0 stands on the transfer ring, 3.06 km above the summit, 150 m wide and 20.5 km round at its mid-width. The
air there averages 8.9 degC and breathes like 2,450 m on Earth, so it is open to everyone. It holds the base zone's
programme (6 km2): markets, food, big halls, the arrivals from the six legs, and some hotels. Two storeys (9 m and
6 m) under a park.

Positions: s is the distance along the ring's mid-width circle (radius R_MID) counter-clockwise from a gate's centre
line; d is the distance in from the ring's outer edge (0 at the frame, 150 at the edge over the hollow).
"""

R_OUTER = 3334.0                  # m, from the form product (ring 0's outer edge)
WIDTH = 150.0
R_MID = R_OUTER - WIDTH / 2
LEVELS = {'L1': 0.0, 'L2': 9.0, 'roof': 15.0}   # m above ring 0's floor (3,063 m above the summit)
QUARTER_LENGTH = 3412.8           # m along R_MID between gate centres (60 degrees)
NODE_SPACING = 853.2              # m along R_MID between frame nodes (15 degrees)
LEG0_DEG = 22.0                   # the legs' meridians, 22 + 60k degrees from east (the map's terrain fit)

# The cross-section, from the outer edge in. Column lines stand on the radial trusses' panel points.
ZONES = [
    ('outer gallery', 0.0, 8.0),   # double height, glazed to the view out through the frame
    ('outer band', 8.0, 62.0),     # halls, hotels, shops
    ('promenade', 62.0, 88.0),     # the ring's street, full height, lit by a skylight in the park
    ('inner band', 88.0, 142.0),   # food, restaurants, hotels facing the hollow
    ('inner gallery', 142.0, 150.0),   # double height, balconies over the hollow
]
COLUMN_LINES = [8.0, 35.0, 62.0, 88.0, 115.0, 142.0]
SKYLIGHT = (69.0, 81.0)            # the promenade's glazed strip in the roof park
UNDERCROFT_LINE = (71.0, 79.0)     # the people mover, in the trusses' depth under the promenade

GATES = [
    # the six gates, where each leg's inclined lifts arrive and its lift spine rises
    dict(n=1, deg=22.0, faces='ENE'), dict(n=2, deg=82.0, faces='N'), dict(n=3, deg=142.0, faces='NW'),
    dict(n=4, deg=202.0, faces='WSW'), dict(n=5, deg=262.0, faces='S'), dict(n=6, deg=322.0, faces='SE'),
]
GATE_HALF = 110.0                  # each gate hall runs 110 m either side of the leg's meridian

QUARTERS = [
    dict(n=1, name='Market Quarter', anchors='Grand Market, craft and flower markets, food hall, bazaar',
         roof='orchards and kitchen gardens'),
    dict(n=2, name='Theatre Quarter', anchors='the Great Hall (20,000 seats), two theatres, cinemas',
         roof='meadow and an open-air amphitheatre'),
    dict(n=3, name='Baths Quarter', anchors='the Baths, courtyard hotels', roof='water gardens'),
    dict(n=4, name='Museum Quarter', anchors='Museum of the Open Moon, library, exhibition halls',
         roof='sculpture and rock gardens'),
    dict(n=5, name='Food Quarter', anchors='food halls, night market, restaurants over the hollow',
         roof='dining terraces and pergolas'),
    dict(n=6, name='Sports Quarter', anchors='sports halls, pools', roof='playing fields and a running track'),
]

# Quarter 1, hall level (L1, with L2 above as noted): places as (name, kind, s0, s1, d0, d1, note).
# kinds: gate, square, court (open to the sky through the park), node, market, crafts, hotel, food, services,
# bazaar, promenade, gallery.
Q1_L1 = [
    ('Gate 1 hall', 'gate', 0, 110, 0, 150, 'arrivals from leg 1; the lift spine; a glass vault rises through the park'),
    ('Gate square', 'square', 110, 190, 8, 142, 'Undercroft Line station below; L2 gallery round'),
    ('Grand Market: stalls', 'market', 190, 640, 8, 62, 'four aisles of stalls; L2 market restaurants'),
    ('Grand Market: nave', 'promenade', 190, 640, 62, 88, 'the promenade as the market street'),
    ('Grand Market: stalls', 'market', 190, 640, 88, 142, 'four aisles of stalls; L2 restaurants over the nave'),
    ('Lane 1', 'court', 640, 700, 8, 142, 'light court open to the park; stairs up'),
    ('Market Hotel lobby, shops', 'hotel', 700, 800, 8, 62, 'rooms on L2'),
    ('Restaurants', 'food', 700, 800, 88, 142, 'facing the hollow'),
    ('Node A room', 'node', 800, 906, 0, 35, 'the outer gallery widens round frame node A'),
    ('Node A court', 'court', 800, 906, 35, 142, 'Undercroft Line station below'),
    ('Craft workshops and shops', 'crafts', 906, 1300, 8, 62, 'makers at L1, studios on L2'),
    ('Flower market', 'crafts', 906, 1300, 88, 142, 'top-lit through the meadow above'),
    ('Lane 2', 'court', 1300, 1360, 8, 142, 'light court'),
    ('Terrace Hotel', 'hotel', 1360, 1650, 8, 62, 'rooms on L2 look out through the frame'),
    ('Terrace Hotel: restaurants, suites', 'hotel', 1360, 1650, 88, 142, 'suites on L2 over the hollow'),
    ('Node B room', 'node', 1650, 1760, 0, 35, ''),
    ('Node B court', 'court', 1650, 1760, 35, 142, 'station below'),
    ('Bakeries, grocers', 'food', 1760, 2100, 8, 62, ''),
    ('Food Hall', 'food', 1760, 2100, 88, 142, 'a 340 m hall of kitchens and tables'),
    ('Lane 3', 'court', 2100, 2160, 8, 142, 'light court'),
    ("Traders' House", 'hotel', 2160, 2500, 8, 62, "traders' hotel and offices"),
    ('Family restaurants', 'food', 2160, 2500, 88, 142, ''),
    ('Node C room', 'node', 2500, 2610, 0, 35, ''),
    ('Node C court', 'court', 2500, 2610, 35, 142, 'station below'),
    ('Tea and Spice Bazaar', 'bazaar', 2610, 3000, 8, 62, 'narrow lanes of small shops'),
    ('Tea and Spice Bazaar', 'bazaar', 2610, 3000, 88, 142, ''),
    ('Lane 4', 'court', 3000, 3060, 8, 142, 'light court'),
    ('Travel services, luggage, clinic', 'services', 3060, 3303, 8, 62, ''),
    ('Waiting lounge, cafes', 'food', 3060, 3303, 88, 142, 'facing the hollow'),
    ('Gate 2 hall', 'gate', 3303, 3413, 0, 150, 'arrivals from leg 2'),
]
Q1_STATIONS = [150, 853, 1706, 2560, 3413]     # Undercroft Line stations (s)

# Quarter 1, the roof park: (name, kind, s0, s1, d0, d1). The outer walk (d 0-10), the skylight walks (d 62-88)
# and the inner walk (d 130-150) run the whole quarter, bridging the lanes' voids.
Q1_ROOF = [
    ('Gate 1 plaza', 'plaza', 0, 190, 10, 130),
    ('Orchard', 'orchard', 190, 640, 10, 62),
    ('Kitchen gardens', 'beds', 190, 640, 88, 130),
    ('Lane 1 void', 'void', 640, 700, 10, 130),
    ('Tea garden', 'lawn', 700, 800, 10, 130),
    ('Node A terrace', 'terrace', 800, 906, 10, 130),
    ('Meadow', 'meadow', 906, 1300, 10, 62),
    ('Meadow', 'meadow', 906, 1300, 88, 130),
    ('Lane 2 void', 'void', 1300, 1360, 10, 130),
    ('Hotel gardens', 'lawn', 1360, 1650, 10, 62),
    ('Reflecting pool', 'water', 1420, 1590, 96, 122),
    ('Hotel gardens', 'lawn', 1360, 1650, 88, 130),
    ('Node B terrace', 'terrace', 1650, 1760, 10, 130),
    ('Pergola terraces', 'pergola', 1760, 2100, 10, 62),
    ('Dining terraces', 'plaza', 1760, 2100, 88, 130),
    ('Lane 3 void', 'void', 2100, 2160, 10, 130),
    ('Play garden', 'play', 2160, 2500, 10, 62),
    ('Splash garden', 'water', 2250, 2400, 96, 122),
    ('Play lawn', 'lawn', 2160, 2500, 88, 130),
    ('Node C terrace', 'terrace', 2500, 2610, 10, 130),
    ('Herb garden', 'beds', 2610, 3000, 10, 62),
    ('Glasshouse', 'pavilion', 2740, 2860, 96, 124),
    ('Herb garden', 'beds', 2610, 3000, 88, 130),
    ('Lane 4 void', 'void', 3000, 3060, 10, 130),
    ('Event lawn', 'lawn', 3060, 3303, 10, 130),
    ('Gate 2 plaza', 'plaza', 3303, 3413, 10, 130),
]
# Kiosk groups (the open-air shops) and pavilions on the roof: (name, s, d, count).
Q1_KIOSKS = [('kiosks', 150, 40, 8), ('kiosks', 150, 110, 6), ('tea house', 750, 70, 1), ('cafe', 828, 60, 1),
             ('kiosks', 1730, 50, 6), ('kiosks', 1930, 110, 10), ('cafe', 2585, 60, 1), ('kiosks', 3180, 70, 12)]
