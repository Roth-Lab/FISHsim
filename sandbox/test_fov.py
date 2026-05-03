import numpy as np
import pandas as pd

from fishsim.codebook import Codebook
from fishsim.data_organisation import DataOrganisation
from fishsim.fov import FieldOfView, ProbedFieldOfView

boundary_box = {"x": (0, 1024), "y": (0, 1024), "z": (-10, 10)}

cell_axes_bounds = {"a": (30, 40), "b": (60, 70), "c": (50, 60)}

fov = FieldOfView(boundary_box, cell_axes_bounds, 20)

print(fov.x_size, fov.y_size, fov.z_size)

codebook_df = pd.DataFrame(
    [
        {"target": "gene_1", "bit_1": 1, "bit_2": 0, "bit_3": 0, "bit_4": 0},
        {"target": "gene_2", "bit_1": 0, "bit_2": 1, "bit_3": 0, "bit_4": 0},
        {"target": "gene_3", "bit_1": 0, "bit_2": 0, "bit_3": 1, "bit_4": 0},
        {"target": "gene_4", "bit_1": 0, "bit_2": 0, "bit_3": 0, "bit_4": 1},
    ]
)

codebook = Codebook(codebook_df)

print(codebook.targets)

print(codebook.target_dist)

print(codebook.get_barcode(1))

data_org_df = pd.DataFrame(
    [
        {"bit_id": "bit_1", "channel": "650", "round": 0},
        {"bit_id": "bit_2", "channel": "750", "round": 0},
        {"bit_id": "bit_3", "channel": "650", "round": 1},
        {"bit_id": "bit_4", "channel": "750", "round": 1},
    ]
)

probed_fov = ProbedFieldOfView(codebook, fov, 1000, bit_add_prob=0, bit_drop_prob=0, sim_nucleus=True)

print(probed_fov.emitter_positions)

print(probed_fov.emitter_targets)

print(probed_fov.emitter_cell_ids)
