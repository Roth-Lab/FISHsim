from fishsim.ellipsoid import Ellipsoid

axes_bounds = {"a": (30, 40), "b": (60, 70), "c": (50, 60)}

center = (0, 0, 0)

axes = Ellipsoid.random_axes(axes_bounds)

print(axes)

e = Ellipsoid(axes, center)

print(e.volume)

print(e.check_overlap(e))

boundary_box_dim = {"x": (0, 1024), "y": (0, 1024), "z": (-10, 10)}

print(e.generate_random_points(boundary_box_dim, (0, 1), num_points=1000))
