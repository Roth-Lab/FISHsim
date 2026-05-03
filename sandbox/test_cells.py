from fishsim.cells import EllipsoidCell

axes_bounds = {"a": (30, 40), "b": (60, 70), "c": (50, 60)}

center = (0, 0, 0)

cell = EllipsoidCell(axes_bounds, center)

print(cell.volume)

print(cell.background(10).shape)

print(cell.check_overlap(cell))

boundary_box = {"x": (0, 1024), "y": (0, 1024), "z": (-10, 10)}

print(cell.generate_emitters(boundary_box, 100))
