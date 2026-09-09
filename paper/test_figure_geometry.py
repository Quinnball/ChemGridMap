import numpy as np
import pandas as pd

from figure_geometry import window_bounds


def test_window_bounds_follow_coordinates_when_outer_columns_are_empty():
    # A 48-column lattice with occupied columns only up to 42 reproduced the bug.
    cells = [(col, row) for col in range(4, 43) for row in range(48)]
    data = pd.DataFrame(cells, columns=["grid_col", "grid_row"])
    data["grid_x"] = data.grid_col / 47
    data["grid_y"] = data.grid_row / 47
    for col, row in [(20, 34), (13, 12)]:
        selected = data.grid_col.between(col, col + 2) & data.grid_row.between(row, row + 2)
        win = {"indices": data.index[selected].tolist()}
        x, y, width, height = window_bounds(data, win)
        assert np.allclose((x, y, width, height),
                           (30 + (col - 0.5) * 318 / 47,
                            65 + (47 - row - 2.5) * 318 / 47,
                            3 * 318 / 47, 3 * 318 / 47))
        px, py = 30 + data.grid_x * 318, 65 + (1 - data.grid_y) * 318
        enclosed = px.between(x, x + width) & py.between(y, y + height)
        assert data.index[enclosed].tolist() == win["indices"]
