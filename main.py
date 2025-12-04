from enum import Enum

from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt

from particle import Particle
from magnetic_field import Magfield
from magnetic_field import FieldType


class Presets(Enum):
    DIPOLE_ELECTRON = 1
    DIPOLE_PROTON = 2


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


def preset(mode: Presets = Presets.DIPOLE_ELECTRON) -> tuple[Particle, Magfield]:
    match mode:
        case Presets.DIPOLE_ELECTRON:
            pos = np.array(
                [
                    2.5e6,
                    2.5e6,
                    2.5e6,
                ]
            )
            vel = np.array(
                [
                    1e5,
                    1e4,
                    1e5,
                ]
            )
            mass = 9.11e-31
            charge = -1.600e-19

            part = Particle(pos, vel, mass, charge)
            mag = Magfield(mode=FieldType.DIPOLE, field_constant=7e8)

            return part, mag

        case Presets.DIPOLE_PROTON:
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

            return part, mag


def solve(part, mag, t_span):
    def func(t, y):
        # dx/dt = v
        # dv/dt = q/m (v X B)
        pos = y[0:3]
        vel = y[3:6]

        acc = part.charge / part.mass * np.cross(vel, mag.magfield(pos))

        dydt = np.concatenate((vel, acc))

        return dydt

    t_points = np.linspace(t_span[0], t_span[1], int(t_span[1] * 1e2))
    initial = np.concatenate((part.pos, part.vel))
    sol = solve_ivp(func, t_span, initial, t_eval=t_points, atol=1e-6, rtol=1e-3)
    return sol


def main():
    mode = Presets.DIPOLE_ELECTRON
    part, mag = preset(mode)

    t_span = [0, 1e4]
    sol = solve(part, mag, t_span)

    t = sol.t
    x, y, z = sol.y[:3]

    file_name = f"{str.lower(mode.name)}_trajectory.csv"
    data = np.vstack([t, x, y, z]).T

    np.savetxt(
        file_name,
        data,
        delimiter=",",
        fmt="%.8f",
        header="time,x_position,y_position,z_position",
    )

    plot(x, y, z)


if __name__ == "__main__":
    main()
