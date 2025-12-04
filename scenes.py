from manim import *
import numpy as np


class TrajectoryAnimation(ThreeDScene):
    def construct(self):
        # Reading data
        data_e = np.loadtxt("dipole_electron_trajectory.csv", delimiter=",", skiprows=1)
        x_e = data_e[:, 1]
        y_e = data_e[:, 2]
        z_e = data_e[:, 3]

        data_p = np.loadtxt("dipole_proton_trajectory.csv", delimiter=",", skiprows=1)
        x_p = data_p[:, 1]
        y_p = data_p[:, 2]
        z_p = data_p[:, 3]

        SAMPLING_STEP = 1
        x_e = x_e[::SAMPLING_STEP]
        y_e = y_e[::SAMPLING_STEP]
        z_e = z_e[::SAMPLING_STEP]

        x_p = x_p[::SAMPLING_STEP]
        y_p = y_p[::SAMPLING_STEP]
        z_p = z_p[::SAMPLING_STEP]

        all_coords_abs_e = np.abs(np.concatenate([x_e, y_e, z_e]))
        all_coords_abs_p = np.abs(np.concatenate([x_p, y_p, z_p]))
        SCALING_FACTOR_E = 3 / np.max(all_coords_abs_e)
        SCALING_FACTOR_P = 7 / np.max(all_coords_abs_p)

        x_e = x_e * SCALING_FACTOR_E
        y_e = y_e * SCALING_FACTOR_E
        z_e = z_e * SCALING_FACTOR_E

        x_p = x_p * SCALING_FACTOR_P
        y_p = y_p * SCALING_FACTOR_P
        z_p = z_p * SCALING_FACTOR_P

        coords_3d_e = np.vstack([x_e, y_e, z_e]).T
        coords_3d_p = np.vstack([x_p, y_p, z_p]).T

        # Camera setup
        self.set_camera_orientation(phi=75 * DEGREES, theta=40 * DEGREES)

        # Create 3D Axes
        axes = ThreeDAxes(
            x_range=[-5, 5, 1],
            y_range=[-5, 5, 1],
            z_range=[-5, 5, 1],
            x_length=7,
            y_length=7,
            z_length=7,
        )
        # axes.add(axes.get_axis_labels(x_label="X", y_label="Y", z_label="Z"))

        self.add(axes)

        # Path

        # Map real coordinates to Manim coordinates
        def coords_to_point(coords):
            return axes.c2p(*coords)

        manim_points_e = [coords_to_point(coords) for coords in coords_3d_e]
        manim_points_p = [coords_to_point(coords) for coords in coords_3d_p]
        manim_points = [manim_points_e, manim_points_p]

        colors = [BLUE, RED]
        paths = VGroup()
        for i, points in enumerate(manim_points):
            path = VMobject()
            path.set_points_as_corners(points)
            path.make_smooth()
            path.set_color(colors[i])
            paths.add(path)

        # Create the particle
        # particle = Dot3D(manim_points[0], radius=0.1, color=RED)
        # trail = TracedPath(
        #     particle.get_center,
        #     stroke_width=4,
        #     stroke_color=BLUE,
        #     dissipating_time=3,
        #     stroke_opacity=[0, 1],
        # )

        # Animation

        self.play(FadeIn(axes, scale=0.5), run_time=1)
        self.begin_ambient_camera_rotation(rate=0.2)
        create_list = [Create(path) for path in paths]
        self.play(*create_list, run_time=10, rate_func=linear)
        # self.add(trail)
        #
        # # Animate the particle moving along the path
        # # The total run_time of this animation determines the speed
        # self.play(
        #     MoveAlongPath(particle, trajectory_path),
        #     run_time=50,
        #     rate_func=linear,  # Use linear for constant speed
        # )

        self.wait(3)
