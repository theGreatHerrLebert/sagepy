import sagepy_connector
psc = sagepy_connector.py_tmt


class Isobaric:
    def __init__(self, type_name: str):
        types = ["tmt6", "tmt10", "tmt11", "tmt16", "tmt18"]
        if type_name in types:
            self.__isobaric_ptr = psc.PyIsobaric(type_name)
        else:
            raise ValueError(f"Invalid isobaric type, allowed values are: {types}")

    @classmethod
    def from_py_isobaric(cls, isobaric: psc.PyIsobaric):
        instance = cls.__new__(cls)
        instance.__isobaric_ptr = isobaric
        return instance

    @property
    def type_name(self):
        return self.__isobaric_ptr.type_name

    def __repr__(self):
        return f"Isobaric({self.__isobaric_ptr.type_name})"

    def get_py_ptr(self):
        return self.__isobaric_ptr
