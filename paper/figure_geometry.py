"""Coordinate-only helpers shared by figure composition and lightweight tests."""


def window_bounds(data, win, extent=318, origin=(30, 65)):
    members = data.iloc[win["indices"]]
    xmin, xmax = members.grid_x.min(), members.grid_x.max()
    ymin, ymax = members.grid_y.min(), members.grid_y.max()
    # Occupied columns need not reach the edge of the candidate lattice.
    xstep = (xmax - xmin) / 2
    ystep = (ymax - ymin) / 2
    return (origin[0] + (xmin - xstep / 2) * extent,
            origin[1] + (1 - ymax - ystep / 2) * extent,
            3 * xstep * extent, 3 * ystep * extent)
