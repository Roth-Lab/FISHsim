import numpy as np
import generate_emitters


def test_emitter_position_bounds():
    x_dim = (0, 100)
    y_dim = (0, 100)
    z_dim = (0, 100)

    # generate_emitters.cell_emitter_position(x_dim, y_dim, z_dim, 1000, 10, , True)
    emitter_positions = generate_emitters.random_emitter_position(
        x_dim, y_dim, z_dim, 10000
    )

    print(emitter_positions)
    is_negative = False
    for pos in emitter_positions:
        if pos[0] < 0 or pos[1] < 0 or pos[2] < 0:
            is_negative = True
            break

    assert is_negative == False


def test_cell_emitter_position_bounds():
    x_dim = (0, 100)
    y_dim = (0, 100)
    z_dim = (0, 100)

    cell_axes_bounds = {"a": [30, 40], "b": [20, 30], "c": [20, 30]}

    emitter_positions = generate_emitters.cell_emitter_position(
        x_dim, y_dim, z_dim, 10000, 10, cell_axes_bounds, True
    )
    print(emitter_positions)
    is_negative = False
    for pos in emitter_positions:
        if pos[0] < 0 or pos[1] < 0 or pos[2] < 0:
            is_negative = True
            break

    assert is_negative == False


if __name__ == "__main__":
    test_emitter_position_bounds()
    test_cell_emitter_position_bounds()
