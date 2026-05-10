import numpy as np


def generate_hexagonal_grid(num_layers, filter_threshold=2):
    """
    Generates a hexagonal grid up to num_layers.
    Filters out points in negative physical quadrants for layers > filter_threshold.

    Returns:
        layer_points: dict mapping layer integer to a list of (q, r, i, j) tuples.
    """
    sqrt3_2 = np.sqrt(3) / 2.0
    layer_points = {layer: [] for layer in range(num_layers + 1)}

    for q in range(-num_layers, num_layers + 1):
        for r in range(-num_layers, num_layers + 1):
            # Hexagon boundary condition
            if -q - r in range(-num_layers, num_layers + 1):
                layer = max(abs(q), abs(r), abs(-q - r))
                i = q + (r / 2.0)
                j = r * sqrt3_2

                layer_points[layer].append((q, r, i, j))

    return layer_points
