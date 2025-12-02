from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt

from particle import Particle
from magnetic_field import Magfield
from magnetic_field import FieldType


def plot(x, y, z):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    ax.plot(
        x,
        y,
        z,
        label="Particle Trajectory",
        linestyle="-",
        markersize=2,
        linewidth=1,
    )

    ax.set_xlabel("X Position")
    ax.set_ylabel("Y Position")
    ax.set_zlabel("Z Position")
    ax.set_title("Trajectory Plot")

    ax.scatter(x[0], y[0], z[0], color="green", s=50, label="Start")
    ax.scatter(x[-1], y[-1], z[-1], color="red", s=50, label="End")
    plt.gca().set_aspect("equal")

    plt.legend()
    plt.show()


def main():
    pos = np.array(
        [
            1e5,
            1e5,
            1e4,
        ]
    )
    vel = np.array(
        [
            1e5,
            1e4,
            1e5,
        ]
    )
    mass = 1.670e-27
    charge = 1.600e-19

    part = Particle(pos, vel, mass, charge)
    mag = Magfield(mode=FieldType.DIPOLE, field_constant=7e8)

    def func(t, y):
        # dx/dt = v
        # dv/dt = q/m (v X B)
        pos = y[0:3]
        vel = y[3:6]

        acc = part.charge / part.mass * np.cross(vel, mag.magfield(pos))

        dydt = np.concatenate((vel, acc))

        return dydt

    t_span = [0, 1e2]
    t_points = np.linspace(t_span[0], t_span[1], int(t_span[1] * 1e2))
    initial = np.concatenate((part.pos, part.vel))
    sol = solve_ivp(func, t_span, initial, t_eval=t_points, atol=1e-6, rtol=1e-3)
    x = sol.y[0, :]
    y = sol.y[1, :]
    z = sol.y[2, :]

    plot(x, y, z)


if __name__ == "__main__":
    main()
