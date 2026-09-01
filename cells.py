"""Where the measured-radar cells are.

Lifted out of `observed.py` so the lightning renderer can import the layout
without importing a GRIB reader with it. One definition, so a cell can never
exist for radar and not for lightning — `test_cells.py` guards it for gaps.
"""

CELL_SIZE = (2400, 1920)
CELL_SPAN = (5.0, 4.0)

# Where the cells are. Hand-placed over the water this app is for — the coasts,
# the Great Lakes and the Gulf bays — rather than tiled across the whole
# country, which would be 160 cells to cover a great deal of Nebraska.
#
# Generous overlap on purpose: a spot near a seam should find a cell holding
# its whole view rather than one that clips it at the edge.
CELL_ORIGINS = [
    # Atlantic, north to south
    ("downeast",      -71.0, 42.0), ("southern-ne",   -74.0, 39.5),
    ("mid-atlantic",  -77.0, 36.5), ("carolinas",     -80.0, 33.0),
    ("sc-georgia",    -82.0, 30.0), ("florida-ne",    -82.0, 27.0),
    ("florida-se",    -83.0, 23.5),
    # Gulf, east to west
    ("florida-sw",    -84.0, 25.5), ("florida-nw",    -87.0, 28.0),
    ("alabama-miss",  -91.0, 28.0), ("louisiana",     -94.0, 27.5),
    ("houston-galv",  -96.0, 28.0), ("texas-central", -98.0, 25.5),
    ("texas-south",  -100.0, 24.0),
    # Pacific
    ("socal",        -121.0, 31.5), ("central-ca",   -124.0, 34.0),
    ("norcal",       -125.0, 37.0), ("oregon",       -126.0, 40.0),
    ("washington",   -126.0, 43.0), ("puget",        -125.0, 46.0),
    # Great Lakes. Superior takes two cells and Erie starts a degree further
    # west than it first did: at -92/45 and -83/41 the layout missed Duluth,
    # Whitefish Bay, Toledo and the whole western basin of Erie — four pieces
    # of heavily boated water, each of which fell back to NOAA's tiles
    # without anything on screen saying why. `test_cells.py` is the guard.
    ("lake-michigan",   -89.0, 41.0),
    ("lake-superior-w", -93.0, 45.0), ("lake-superior-e", -88.0, 45.5),
    ("lake-huron",      -85.0, 43.0), ("lake-erie",       -84.0, 41.0),
    ("lake-ontario",    -79.0, 42.5),
]
