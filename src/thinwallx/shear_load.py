"""Shear load representation for ThinWallX v0.2."""

from __future__ import annotations

from dataclasses import dataclass
import math
import numpy as np

from thinwallx.exceptions import GeometryError


@dataclass(frozen=True)
class ShearLoad:
    """Applied transverse shear force vector V = [Vx, Vy]^T.

    Attributes:
        vx: Transverse shear force in global +x direction.
        vy: Transverse shear force in global +y direction.
    """

    vx: float
    vy: float

    def __post_init__(self) -> None:
        """Validate shear load components for finiteness."""
        if not math.isfinite(self.vx) or not math.isfinite(self.vy):
            raise GeometryError(
                f"Shear load components must be finite. Got vx={self.vx}, vy={self.vy}."
            )

    @property
    def vector(self) -> np.ndarray:
        """Return [vx, vy] as a 1D numpy array."""
        return np.array([self.vx, self.vy], dtype=float)

    def __add__(self, other: ShearLoad) -> ShearLoad:
        if not isinstance(other, ShearLoad):
            return NotImplemented
        return ShearLoad(self.vx + other.vx, self.vy + other.vy)

    def __sub__(self, other: ShearLoad) -> ShearLoad:
        if not isinstance(other, ShearLoad):
            return NotImplemented
        return ShearLoad(self.vx - other.vx, self.vy - other.vy)

    def __mul__(self, scalar: float) -> ShearLoad:
        if not isinstance(scalar, (int, float)):
            return NotImplemented
        return ShearLoad(self.vx * scalar, self.vy * scalar)

    def __rmul__(self, scalar: float) -> ShearLoad:
        return self.__mul__(scalar)

    def __neg__(self) -> ShearLoad:
        return ShearLoad(-self.vx, -self.vy)
