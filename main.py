from scipy.integrate import solve_ivp
import numpy as np

from particle import Particle
from magnetic_field import Magfield
from magnetic_field import FieldType


def main():
    pos = np.array([1, 0, 0])
    vel = np.array([1, 0, 1])
    mass = 1
    charge = 1

    part = Particle(pos, vel, mass, charge)
    mag = Magfield()

    def func(t, y):
        # dx/dt = v
        # dv/dt = q/m (v X B)
        pos = y[0:3]
        vel = y[3:6]

        acc = part.charge / part.mass * np.cross(vel, mag.magfield(pos))

        dydt = np.concatenate((vel, acc))

        return dydt

    t_span = [0, 100]
    t_points = np.linspace(t_span[0], t_span[1], 1000)
    sol = solve_ivp(func, t_span, part.position, t_eval=t_points)
    x = sol.y[0, :]
    y = sol.y[1, :]
    z = sol.y[2, :]

    print(x)


if __name__ == "__main__":
    main()
