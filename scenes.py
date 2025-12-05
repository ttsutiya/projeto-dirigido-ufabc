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


class LorentzDerivation(Scene):
    def construct(self):
        full_lorentz = MathTex(
            r"\mathbf{F} = q(\mathbf{E} + \mathbf{v} \times \mathbf{B})"
        )
        full_lorentz.to_edge(UP).shift(DOWN * 0.5)
        self.play(Write(full_lorentz))

        magnetic_force = MathTex(r"\mathbf{F} = q( \mathbf{v} \times \mathbf{B})")
        magnetic_force.to_edge(UP).shift(DOWN * 0.5)
        self.play(ReplacementTransform(full_lorentz, magnetic_force))

        newton_second_law = MathTex(r"\mathbf{F}" r"=m\mathbf{a}")
        newton_second_law.to_edge(UP).shift(DOWN * 1.5)
        self.play(Write(newton_second_law))

        merged_equation = MathTex(
            r"m\mathbf{a}", r"=", r"q(\mathbf{v} \times \mathbf{B})"
        )

        self.play(
            ReplacementTransform(newton_second_law, merged_equation),
            ReplacementTransform(magnetic_force, merged_equation),
        )

        acc = MathTex(r"\mathbf{a}", r"=", r"\frac{q}{m}(\mathbf{v} \times \mathbf{B})")
        self.play(ReplacementTransform(merged_equation, acc))

        dxdt = MathTex(
            r"\frac{d\mathbf{x}}{dt}",
            r"=",
            r"\mathbf{v}",
        ).shift(UP * 1.5)
        dvdt = MathTex(
            r"\frac{d\mathbf{v}}{dt}",
            r"=",
            r"\mathbf{a}",
        )

        self.play(
            Write(dxdt),
            Write(dvdt),
            acc.animate.shift(DOWN * 1.5),
        )

        merge_dvdt = MathTex(
            r"\frac{d\mathbf{v}}{dt}",
            r"=",
            r"\frac{q}{m}(\mathbf{v} \times \mathbf{B})",
        )
        merge_dvdt
        self.play(ReplacementTransform(dvdt, merge_dvdt), FadeOut(acc))


class NumericalMethod1(Scene):
    def construct(self):
        axes = Axes(
            x_range=[0, 10, 1],
            y_range=[0, 10, 1],
            x_length=7,
            y_length=5,
        )

        labels = axes.get_axis_labels(x_label=r"t", y_label=r"y")
        self.wait()
        self.play(Create(axes), Create(labels))

        x_unit = axes.get_x_axis().get_unit_size()
        y_unit = axes.get_y_axis().get_unit_size()

        def func(t):
            return 0.1 * t**2 + 1

        def derivative_func(t):
            return 0.2 * t

        exact_curve = axes.plot(func, color=BLUE)
        self.play(Create(exact_curve))

        # Euler method
        euler_title = Text("Método de Euler", color=RED, font_size=40).to_edge(UP)
        self.play(Write(euler_title))
        self.wait()

        t0 = 2.0
        y0 = func(t0)
        h = 5.0

        start_coords = axes.coords_to_point(t0, y0)
        start_dot = Dot(start_coords, color=YELLOW)
        start_label = MathTex(r"(t_i, y_i)", font_size=30).next_to(
            start_dot, UP + LEFT, buff=0.1
        )
        self.wait()

        h_line = Line(
            axes.coords_to_point(t0, 0.5),
            axes.coords_to_point(t0 + h, 0.5),
            color=RED,
            stroke_width=2,
        )
        h_label = MathTex(r"h", font_size=35).next_to(h_line, UP, buff=0.1)

        self.add(start_dot)
        self.play(Write(start_label))
        self.play(Create(h_line), Write(h_label))
        self.wait()

        k1 = derivative_func(t0)
        self.wait()

        k1_angle = np.arctan(k1 * y_unit / x_unit)
        k1_line = Line(
            start_dot.get_center()
            + 0.5 * np.array([-np.cos(k1_angle), -np.sin(k1_angle), 0]),
            start_dot.get_center()
            + 0.5 * np.array([np.cos(k1_angle), np.sin(k1_angle), 0]),
            color=PURPLE_A,
        )
        k1_label = MathTex(r"k1", font_size=20).next_to(k1_line, DOWN, buff=0.1)

        self.play(Create(k1_line), Write(k1_label))
        self.wait()

        t1 = t0 + h
        y1_euler = y0 + h * k1

        euler_step_line = Line(
            start_dot.get_center(), axes.coords_to_point(t1, y1_euler), color=RED
        )
        euler_dot = Dot(axes.coords_to_point(t1, y1_euler), color=RED)
        euler_label = MathTex(r"(t_{i+1}, y_{i+1})", font_size=30).next_to(
            euler_dot, RIGHT, buff=0.1
        )

        self.play(Create(euler_step_line), Create(euler_dot), Write(euler_label))
        self.wait()

        y1_exact = func(t1)
        exact_dot_t1 = Dot(axes.coords_to_point(t1, y1_exact), color=YELLOW)

        error_line = DashedLine(
            euler_dot.get_center(), exact_dot_t1.get_center(), color=GRAY
        )
        error_label = MathTex(
            r"Erro",  # Usa \text{} para o texto não-matemático
            color=RED_E,
            font_size=30,
        ).next_to(error_line, RIGHT, buff=0.2)

        self.play(Create(exact_dot_t1))
        self.play(Create(error_line), Write(error_label))
        self.wait()

        self.play(
            FadeOut(euler_step_line),
            FadeOut(euler_dot),
            FadeOut(euler_label),
            FadeOut(error_line),
            FadeOut(error_label),
            FadeOut(exact_dot_t1),
            FadeOut(euler_title),
            run_time=0.8,
        )


class NumericalMethod2(Scene):
    def construct(self):
        axes = Axes(
            x_range=[0, 10, 1],
            y_range=[0, 10, 1],
            x_length=7,
            y_length=5,
        )

        labels = axes.get_axis_labels(x_label=r"t", y_label=r"y")
        self.add(axes, labels)

        x_unit = axes.get_x_axis().get_unit_size()
        y_unit = axes.get_y_axis().get_unit_size()

        def func(t):
            return 0.1 * t**2 + 1

        def derivative_func(t):
            return 0.2 * t

        exact_curve = axes.plot(func, color=BLUE)
        self.add(exact_curve)

        t0 = 2.0
        y0 = func(t0)
        h = 5.0
        start_coords = axes.coords_to_point(t0, y0)
        start_dot = Dot(start_coords, color=YELLOW)
        start_label = MathTex(r"(t_i, y_i)", font_size=30).next_to(
            start_dot, UP + LEFT, buff=0.1
        )

        h_line = Line(
            axes.coords_to_point(t0, 0.5),
            axes.coords_to_point(t0 + h, 0.5),
            color=RED,
            stroke_width=2,
        )
        h_label = MathTex(r"h", font_size=35).next_to(h_line, UP, buff=0.1)

        self.add(start_dot)

        k1 = derivative_func(t0)

        k1_angle = np.arctan(k1 * y_unit / x_unit)
        k1_line = Line(
            start_dot.get_center()
            + 0.5 * np.array([-np.cos(k1_angle), -np.sin(k1_angle), 0]),
            start_dot.get_center()
            + 0.5 * np.array([np.cos(k1_angle), np.sin(k1_angle), 0]),
            color=PURPLE_A,
        )
        k1_label = MathTex(r"k1", font_size=20).next_to(k1_line, DOWN, buff=0.1)
        self.add(k1_line, k1_label)

        t1 = t0 + h

        # RK2
        self.wait()
        rk2_title = Text("Runge-Kutta 2", color=RED, font_size=40).to_edge(UP)
        self.play(Write(rk2_title))
        self.wait()

        t_mid = t0 + h / 2
        y_k1_mid = y0 + (h / 2) * k1

        mid_dot = Dot(axes.coords_to_point(t_mid, y_k1_mid), color=PURPLE)
        mid_estimate_line = DashedLine(
            start_dot.get_center(), mid_dot.get_center(), color=PURPLE
        )

        self.play(Create(mid_estimate_line), Create(mid_dot))

        k2 = derivative_func(t_mid)

        k2_angle = np.arctan(k2 * y_unit / x_unit)
        k2_line = Line(
            mid_dot.get_center()
            + 0.5 * np.array([-np.cos(k2_angle), -np.sin(k2_angle), 0]),
            mid_dot.get_center()
            + 0.5 * np.array([np.cos(k2_angle), np.sin(k2_angle), 0]),
            color=PURPLE_A,
        )
        k2_label = MathTex(r"k2", font_size=20).next_to(k2_line, DOWN, buff=0.1)

        self.play(Create(k2_line), Write(k2_label))
        self.wait()

        k_avg = 0.5 * (k1 + k2)

        y1_rk2 = y0 + h * k_avg
        rk2_dot = Dot(axes.coords_to_point(t1, y1_rk2), color=GREEN)
        rk2_label = MathTex(r"(t_{i+1}, y_{i+1})", font_size=30).next_to(
            rk2_dot, RIGHT, buff=0.1
        )

        rk2_final_line = Line(start_dot.get_center(), rk2_dot.get_center(), color=GREEN)

        self.play(
            Create(rk2_final_line),
            Create(rk2_dot),
            Create(rk2_label),
        )
        self.wait()

        y1_exact = func(t1)
        exact_dot_t1 = Dot(axes.coords_to_point(t1, y1_exact), color=YELLOW)

        rk2_error_line = DashedLine(
            rk2_dot.get_center(), exact_dot_t1.get_center(), color=GRAY_A
        )
        k2_error_label = MathTex(
            r"Erro",  # Usa \text{} para o texto não-matemático
            color=RED_E,
            font_size=30,
        ).next_to(rk2_error_line, RIGHT, buff=0.3)

        self.play(Create(exact_dot_t1))
        self.play(Create(rk2_error_line), Write(k2_error_label))
        self.wait()

        self.play(FadeOut(Group(*self.mobjects)))
        self.wait()


class Intro(Scene):
    def construct(self):
        main_title = Text("Trajetória de Partículas", font_size=50, color=BLUE_A).shift(
            UP * 1.0
        )
        main_title2 = Text("em Campos Magnéticos", font_size=50, color=BLUE_A).shift(
            UP * 0.25
        )

        self.play(Write(main_title))
        self.play(Write(main_title2))
        self.wait(1)

        subtitle = Text("Explicação Qualitativa", font_size=30, color=GREEN_C)
        subtitle.next_to(main_title2, DOWN, buff=0.5)

        self.play(Write(subtitle), run_time=1.5)
        self.wait(1)

        self.play(
            FadeOut(main_title), FadeOut(main_title2), FadeOut(subtitle), run_time=1
        )
        self.wait()


class Conteudo(Scene):
    def construct(self):
        list_title = Text("Tópicos", font_size=50, color=YELLOW_A)
        list_title.to_edge(UP)
        self.wait()
        self.play(Write(list_title), run_time=1)
        self.wait(0.5)

        bullet_points = VGroup(
            Tex(r"\text{Força de Lorentz}", font_size=40),
            Tex(r"\text{Método de Euler}", font_size=40),
            Tex(r"\text{Runge-Kutta 2}", font_size=40),
            Tex(r"\text{Simulação Runge-Kutta 4}", font_size=40),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.8)

        bullet_points.next_to(list_title, DOWN, buff=2)

        for i, point in enumerate(bullet_points):
            self.play(FadeIn(point, shift=LEFT), run_time=0.8)
            self.wait(1.5)

        self.wait(2)
        self.play(FadeOut(VGroup(list_title, bullet_points)), run_time=1)
