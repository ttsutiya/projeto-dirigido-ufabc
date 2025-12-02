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
        origin=np.array([0, 0, 0]),
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

    def magfield(self, pos):
        magfield = self._field_calculator(pos)
        return magfield

    def _distance(self, pos):
        diff = pos - self.origin
        r = np.linalg.norm(diff)
        return r

    def _calcualtor_constant(self, pos):
        magfield = np.array([0, 0, self.field_constant])
        return magfield

    def _calcualtor_inv_r(self, pos):
        r = pos[2]  # distance form the plane z=0

        if not r:
            return None

        magfield = np.array([0, 0, self.field_constant * 1 / r])
        return magfield

    def _calcualtor_inv_r2(self, pos):
        r = pos[2]  # distance form the plane z=0

        if not r:
            return None

        magfield = np.array([0, 0, self.field_constant * 1 / pow(r, 2)])
        return magfield

    def _calcualtor_dipole(self, pos):
        r = self._distance(pos)

        if not r:
            return None

        Bx = (3 * self.field_constant * pos[0] * pos[2]) / pow(r, 5)
        By = (3 * self.field_constant * pos[1] * pos[2]) / pow(r, 5)
        Bz = (self.field_constant * (3 * pow(pos[2], 2) - pow(r, 2))) / pow(r, 5)

        magfield = np.array([Bx, By, Bz])
        return magfield
