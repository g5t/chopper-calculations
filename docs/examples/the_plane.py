"""Why the plane is drawn in inverse velocity rather than wavelength."""
# --8<-- [start:proportional]
from chopcal import constants
from chopcal.lib import inverse_velocity_to_wavelength

# Wavelength is proportional to inverse velocity: lambda = (h/m) * a, with no other terms.
for a in (1e-4, 5e-4, 1e-2):
    print(f'{a:.0e} s/m -> {inverse_velocity_to_wavelength(a):.6f} AA'
          f'   ratio {inverse_velocity_to_wavelength(a) / a:.4f}')

print(f'H_OVER_M = {constants.H_OVER_M:.4f} AA m/s')
# --8<-- [end:proportional]

# --8<-- [start:arrival]
# A disc 78 m downstream, and a neutron emitted at 1 ms with inverse velocity 5e-4 s/m.
path, emission, a = 78.0, 1e-3, 5e-4

# The scale is h/m, whichever way it is spelled; ask the library for its own.
h_over_m = inverse_velocity_to_wavelength(1.0)

in_inverse_velocity = emission + path * a
in_wavelength = emission + path / h_over_m * inverse_velocity_to_wavelength(a)

print(f'{in_inverse_velocity * 1e3:.6f} ms   {in_wavelength * 1e3:.6f} ms')
# --8<-- [end:arrival]


def main():
    # A pure scaling: the ratio does not depend on where you evaluate it.
    ratios = [inverse_velocity_to_wavelength(a) / a
              for a in (1e-4, 5e-4, 1.26e-3, 1e-2)]
    assert all(abs(r - ratios[0]) < 1e-6 for r in ratios), 'not a scaling'

    # And the scale is h/m, to the 1.1e-9 the library's V2K and K2V differ from the
    # SI-derived pair by -- see the warning on the constants page.
    assert abs(ratios[0] - constants.H_OVER_M) / constants.H_OVER_M < 2e-9

    # So both coordinates give the same arrival time; only the coefficient differs, and
    # with the library's own h/m on both sides there is nothing left over at all.
    assert in_inverse_velocity == in_wavelength

    # Reaching for chopcal.constants.H_OVER_M instead leaves the 1.4e-9 those two
    # definitions of h/m differ by -- 55 ps on a 40 ms flight. Small, and it is the whole
    # reason the conversions are exposed rather than left to the caller.
    mixed = emission + path / constants.H_OVER_M * inverse_velocity_to_wavelength(a)
    assert 0 < abs(mixed - in_inverse_velocity) < 1e-10


main()
