import attrs

@attrs.define
class MissingValuesError(TypeError):
    message : str
    a: set