import numpy as np

from fishsim.cells import EllipsoidCell


class FieldOfView(object):
    def __init__(self, boundary_box, cell_axes_bounds, num_cells, rng):
        self.rng = rng

        self.boundary_box = boundary_box

        self._init_cells(cell_axes_bounds, num_cells)

    @property
    def x_size(self):
        return self.boundary_box["x"][1] - self.boundary_box["x"][0]

    @property
    def y_size(self):
        return self.boundary_box["y"][1] - self.boundary_box["y"][0]

    @property
    def z_size(self):
        return self.boundary_box["z"][1] - self.boundary_box["z"][0]

    def _init_cells(self, cell_axes_bounds, num_cells, max_attempt_multiple=10):
        self.cells = []

        attempt = 0

        max_attempts = max_attempt_multiple * num_cells

        while (len(self.cells) < num_cells) and (attempt < max_attempts):
            attempt += 1

            # Generate center position for cells
            cell_pos = [
                self.rng.uniform(*self.boundary_box["x"]),
                self.rng.uniform(*self.boundary_box["y"]),
                self.rng.uniform(*self.boundary_box["z"]),
            ]

            current_cell = EllipsoidCell(cell_axes_bounds, cell_pos)

            overlap = False

            for cell in self.cells:
                if current_cell.check_overlap(cell):
                    overlap = True
                    break

            if not overlap:
                self.cells.append(current_cell)


class ProbedFieldOfView(object):
    def __init__(
        self, codebook, fov, num_emitters, rng, bit_add_prob=0, bit_drop_prob=0, sim_nucleus=True, subpixel=True
    ):
        self.codebook = codebook

        self.fov = fov

        self.rng = rng

        self._init_emitters(num_emitters, sim_nucleus, subpixel)

        self._init_emitter_barcodes(bit_add_prob, bit_drop_prob)

    @property
    def cells(self):
        return self.fov.cells

    @property
    def emitter_cell_ids(self):
        cell_ids = []

        for i, c in enumerate(self.cells):
            for _ in c.emitter_target_ids:
                cell_ids.append(i)

        return np.array(cell_ids)

    @property
    def emitter_positions(self):
        return np.concat([c.emitter_positions for c in self.cells if len(c.emitter_positions) > 0])

    @property
    def emitter_target_ids(self):
        return np.concat([c.emitter_target_ids for c in self.cells if len(c.emitter_target_ids) > 0])

    @property
    def emitter_targets(self):
        targets = []

        for c in self.cells:
            for t in c.emitter_target_ids:
                targets.append(self.codebook.targets[t])

        return np.array(targets)

    def _init_emitters(self, num_emitters, sim_nucleus, subpixel):
        p = self.rng.dirichlet([c.volume for c in self.cells])

        cell_emitter_count = self.rng.multinomial(num_emitters, p)

        for i, cell in enumerate(self.cells):
            emitter_positions = []

            emitter_target_ids = []

            p = self.rng.dirichlet(self.codebook.target_dist + 1e-6)

            target_emitter_counts = self.rng.multinomial(cell_emitter_count[i], p)

            for t in range(self.codebook.num_targets):
                emitter_positions.extend(
                    cell.generate_emitters(
                        self.fov.boundary_box,
                        target_emitter_counts[t],
                        sim_nucleus=sim_nucleus,
                    )
                )

                emitter_target_ids.extend(t * np.ones(target_emitter_counts[t]))

            emitter_positions = np.array(emitter_positions)

            if not subpixel:
                emitter_positions = np.floor(emitter_positions).astype(int)

            cell.emitter_positions = emitter_positions

            cell.emitter_target_ids = np.array(emitter_target_ids, dtype=int)

    def _init_emitter_barcodes(self, bit_add_prob, bit_drop_prob):
        barcodes = []

        for i in self.emitter_target_ids:
            b = self.codebook.get_barcode(i).copy()

            if self.rng.random() < bit_add_prob:
                idx = self.rng.choice(np.where(b == 0)[0])

                b[idx] = 1

            elif self.rng.random() < bit_drop_prob:
                idx = self.rng.choice(np.where(b == 1)[0])

                b[idx] = 0

            barcodes.append(b)

        self.emitter_barcodes = np.array(barcodes)
