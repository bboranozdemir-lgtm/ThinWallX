"""Custom exceptions for Sectalix."""

from __future__ import annotations


class SectalixError(Exception):
    """Base exception class for all Sectalix errors."""


class ValidationError(SectalixError):
    """Raised when geometric or topological validation fails."""


class GeometryError(ValidationError):
    """Raised when geometric properties are invalid (e.g. non-finite coordinates, zero length, non-positive thickness)."""


class TopologyError(ValidationError):
    """Raised when section topology is invalid (e.g. empty section, disconnected segments, closed loops, duplicate segments)."""


class SingularSectionError(ValidationError):
    """Raised when the section coordinate second-moment matrix C is singular or rank-deficient."""

