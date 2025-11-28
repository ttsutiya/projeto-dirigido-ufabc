import numpy as np


class Particle:
    def __init__(
        self,
        pos: np.array = np.array([0, 0, 0]),
        vel: np.array = np.array([0, 0, 1]),
        mass: float = 1,
        charge: float = 1,
    ):
        self.intial_pos = pos
        self.intial_vel = vel

        self.pos = self.initial_pos
        self.vel = self.initial_vel

        self.mass = mass
        self.charge = charge
