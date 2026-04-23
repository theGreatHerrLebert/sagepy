import warnings
from typing import List

import sagepy_connector
psc = sagepy_connector.py_mass


def monoisotopic(aa: str) -> float:
    return psc.py_monoisotopic(aa)


class Composition:
    def __init__(self, carbon, sulfur):
        """Composition class

        Args:
            carbon (int): The number of carbon atoms
            sulfur (int): The number of sulfur atoms
        """
        self.__composition_ptr = psc.PyComposition(carbon, sulfur)

    @classmethod
    def from_py_composition(cls, composition: psc.PyComposition):
        instance = cls.__new__(cls)
        instance.__composition_ptr = composition
        return instance

    @property
    def carbon(self):
        return self.__composition_ptr.carbon

    @property
    def sulfur(self):
        return self.__composition_ptr.sulfur

    def __repr__(self):
        return f"Composition(carbon: {self.carbon}, sulfur: {self.sulfur})"

    def get_py_ptr(self):
        return self.__composition_ptr

    @staticmethod
    def sum(composition_list: List['Composition']) -> 'Composition':
        return psc.PyComposition.sum([c.get_py_ptr() for c in composition_list])

    @staticmethod
    def aa_composition(aa: str) -> 'Composition':
        return Composition.from_py_composition(psc.PyComposition.py_composition(aa))


class CONSTANTS:
    """Mass constants for common molecules and particles."""
    NEUTRON: float = psc.neutron()
    PROTON: float = psc.proton()
    H2O: float = psc.h2o()
    NH3: float = psc.nh3()


class Tolerance:
    def __init__(self, da: (float, float) = None, ppm: (float, float) = None):
        """Tolerance class

        Args:
            da (float, optional): Signed (low, high) tolerance in Da. The
                conventional symmetric form is e.g. ``(-0.5, 0.5)``. Defaults
                to None.
            ppm (float, optional): Signed (low, high) tolerance in ppm. The
                conventional symmetric form is e.g. ``(-10.0, 10.0)``. Note
                that this is a SIGNED interval — passing ``(10.0, 10.0)``
                yields a zero-width window that never matches anything.
                Defaults to None.

        Raises:
            ValueError: if both / neither of da and ppm are provided, or if
                low > high (reversed interval).
        """
        if da is not None and ppm is not None:
            raise ValueError("Only one of da or ppm can be set")
        elif da is None and ppm is None:
            raise ValueError("One of da or ppm must be set")

        # Validate the (low, high) interval for whichever unit was supplied.
        # A reversed interval is always a bug; a zero-width interval almost
        # always is too — historically this footgun manifested as silent
        # 0-PSM searches when callers wrote `Tolerance(ppm=(20, 20))` meaning
        # "±20 ppm". Loud-warn on equal bounds, hard-error on swapped bounds.
        bounds = ppm if ppm is not None else da
        unit   = "ppm" if ppm is not None else "da"
        try:
            lo, hi = float(bounds[0]), float(bounds[1])
        except (TypeError, IndexError, ValueError) as e:
            raise ValueError(
                f"{unit} must be a (low, high) tuple of two numbers"
            ) from e
        if lo > hi:
            raise ValueError(
                f"Tolerance({unit}=({lo}, {hi})): low > high — the interval "
                f"is reversed; expected low <= high "
                f"(e.g. {unit}=({-abs(hi)}, {abs(hi)}) for ±{abs(hi)})"
            )
        if lo == hi:
            warnings.warn(
                f"Tolerance({unit}=({lo}, {hi})): zero-width interval — "
                f"this matches nothing. Did you mean "
                f"{unit}=({-abs(hi)}, {abs(hi)})?",
                UserWarning,
                stacklevel=2,
            )

        self.__tolerance_ptr = psc.PyTolerance(da, ppm)

    def get_py_ptr(self):
        return self.__tolerance_ptr

    @classmethod
    def from_py_tolerance(cls, tolerance: psc.PyTolerance) -> 'Tolerance':
        instance = cls.__new__(cls)
        instance.__tolerance_ptr = tolerance
        return instance

    @property
    def da(self) -> (float, float):
        return self.__tolerance_ptr.da

    @property
    def ppm(self) -> (float, float):
        return self.__tolerance_ptr.ppm

    def bounds(self, center: float) -> (float, float):
        return self.__tolerance_ptr.bounds(center)

    def contains(self, center: float, target: float) -> bool:
        return self.__tolerance_ptr.contains(center, target)

    def __repr__(self) -> str:
        if self.da is not None:
            return f"Tolerance(da={self.da})"
        else:
            return f"Tolerance(ppm={self.ppm})"

    def __mul__(self, other) -> 'Tolerance':
        if isinstance(other, float):
            return Tolerance.from_py_tolerance(self.__tolerance_ptr * other)

        elif isinstance(other, int):
            return Tolerance.from_py_tolerance(self.__tolerance_ptr * float(other))

        else:
            raise ValueError("Tolerance can only be multiplied by a float or an int")
