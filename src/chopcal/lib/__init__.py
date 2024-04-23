from chopcal._chopper_lib_impl import Chopper, inverse_velocity_windows, inverse_velocity_limits, wavelength_limits


def _chopper_str(chopper):
    return f"Chopper[{chopper.speed} Hz, {chopper.phase} deg, {chopper.angle} deg, {chopper.distance} m]"


Chopper.__str__ = _chopper_str


all = [
    'Chopper',
    'inverse_velocity_windows',
    'inverse_velocity_limits',
    'wavelength_limits',
]