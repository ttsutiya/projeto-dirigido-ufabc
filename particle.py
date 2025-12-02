import numpy as np

# create settings for electron and proton


class Particle:
    def __init__(
        self,
        pos=np.array([0, 0, 0]),
        vel=np.array([0, 0, 1]),
        mass: float = 1,
        charge: float = 1,
    ):
        self.pos = pos
        self.vel = vel

        self.mass = mass
        self.charge = charge
