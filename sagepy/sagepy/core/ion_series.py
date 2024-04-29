import sagepy_connector

from sagepy.core.peptide import Peptide

psc = sagepy_connector.py_ion_series


class IonType:
    def __init__(self, ion_type: str):
        """IonType class

        Args:
            ion_type (str): The ion type, allowed values are: a, b, c, x, y, z
        """
        try:
            self.__ion_type_ptr = psc.PyKind(ion_type)
        except ValueError:
            raise ValueError("Invalid ion type, allowed values are: a, b, c, x, y, z")

    @classmethod
    def from_py_kind(cls, kind: psc.PyKind):
        instance = cls.__new__(cls)
        instance.__ion_type_ptr = kind
        return instance

    @classmethod
    def y(cls):
        return cls.from_py_kind(psc.PyKind('y'))

    @classmethod
    def b(cls):
        return cls.from_py_kind(psc.PyKind('b'))

    def __repr__(self):
        return f"IonType({self.__ion_type_ptr.kind_as_string()})"

    def __hash__(self):
        return hash(self.__ion_type_ptr.kind_as_string())

    def __eq__(self, other):
        if not isinstance(other, IonType):
            return False
        return self.__ion_type_ptr.kind_as_string() == other.__ion_type_ptr.kind_as_string()

    def get_py_ptr(self):
        return self.__ion_type_ptr

    def to_string(self) -> str:
        return self.__ion_type_ptr.kind_as_string()


class Ion:
    """Ion class

    Args:
        ion_type (IonType): The ion type, e.g. b, y
        mass (float): The mass of the ion
    """
    def __init__(self, ion_type: IonType, mass: float):
        self.__ion_ptr = psc.PyIon(ion_type.get_py_ptr(), mass)

    @classmethod
    def from_py_ion(cls, ion: psc.PyIon):
        instance = cls.__new__(cls)
        instance.__ion_ptr = ion
        return instance

    @property
    def ion_type(self):
        return IonType.from_py_kind(self.__ion_ptr.kind)

    @property
    def mono_isotopic_mass(self):
        return self.__ion_ptr.monoisotopic_mass

    def __repr__(self):
        return f"Ion({self.ion_type}, {self.mono_isotopic_mass})"


class IonSeries:
    """IonSeries class

    Args:
        peptide (Peptide): The peptide
        ion_type (IonType): The ion type, e.g. b, y
    """
    def __init__(self, peptide: Peptide, ion_type: IonType):
        self.__ion_series_ptr = psc.PyIonSeries(peptide.get_py_ptr(), ion_type.get_py_ptr())

    @classmethod
    def from_py_ion_series(cls, ion_series: psc.PyIonSeries):
        instance = cls.__new__(cls)
        instance.__ion_series_ptr = ion_series
        return instance

    @property
    def ion_type(self):
        return IonType.from_py_kind(self.__ion_series_ptr.kind)

    @property
    def cumulative_mass(self):
        return self.__ion_series_ptr.cumulative_mass

    @property
    def peptide(self):
        return Peptide.from_py_peptide(self.__ion_series_ptr.peptide)

    def __repr__(self):
        return f"IonSeries({self.ion_type}, {self.cumulative_mass}, {self.peptide})"

    def get_py_ptr(self):
        return self.__ion_series_ptr

    def get_ion_series(self):
        return [Ion.from_py_ion(i) for i in self.__ion_series_ptr.get_ion_series()]