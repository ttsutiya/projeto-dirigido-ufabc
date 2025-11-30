from enum import Enum
import numpy as np


class FieldType(Enum):
    CONSTANT = 1  # B0 z_hat
    INV_R = 2  # B0 * 1/r z_hat (r is the distance from z plane)
    INV_R2 = 3  # B0 * 1/r^2 z_hat (r is the distance from z plane)
    DIPOLE = 4


class Magfield:
    def __init__(
        self,
        mode: FieldType = FieldType.CONSTANT,
        field_constant: float = 1,
        origin: np.array = np.array([0, 0, 0]),
    ):
        self.mode = mode
        self.origin = origin
        self.field_constant = field_constant

        field_calculators = {
            FieldType.CONSTANT: self._calcualtor_constant,
            FieldType.INV_R: self._calcualtor_inv_r,
            FieldType.INV_R2: self._calcualtor_inv_r2,
            FieldType.DIPOLE: self._calcualtor_dipole,
        }

        self._field_calculator = field_calculators[self.mode]

    def magfield(self, position: np.array):
        magfield = self._field_calculator(position)
        return magfield

    def _distance(self, position: np.array) -> float:
        diff = position - self.origin
        r = np.linalg.norm(diff)
        return r

    def _calcualtor_constant(self, position: np.array):
        return self.field_constant

    def _calcualtor_inv_r(self, position: np.array):
        r = position[3]  # distance form the plane z=0

        if not r:
            return None

        return self.field_cosntant * 1 / r

    def _calcualtor_inv_r2(self, position: np.array):
        r = position[3]  # distance form the plane z=0

        if not r:
            return None

        return self.field_cosntant * 1 / pow(r, 2)

    def _calcualtor_dipole(self, position: np.array):
        # B[0] =  (3*initialB * pos[0] * pos[2]) / pow(r,5);
        # B[1] =  (3*initialB * pos[1] * pos[2]) / pow(r,5);
        # B[2] =  (initialB * (3*pow(pos[2],2) - pow(r,2))) / pow(r,5);
        pass
