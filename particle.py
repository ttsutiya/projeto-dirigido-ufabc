import numpy as np


class Particle:
    def __init__(
        self,
        position: np.array = np.array([0, 0, 0]),
        velocity: np.array = np.array([0, 0, 1]),
        mass: float = 1,
        charge: float = 1,
    ):
        self.position = position
        self.velocity = velocity

        self.mass = mass
        self.charge = charge
