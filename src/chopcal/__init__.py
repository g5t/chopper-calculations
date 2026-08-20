from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("chopcal")
except PackageNotFoundError:
    # package is not installed
    pass

del version
del PackageNotFoundError

from chopcal._chopcal_impl import bifrost as _bifrost
import chopcal.lib as lib

# The band the BIFROST chopper train passes, fixed by the 14 Hz source period and the
# 155.7 m between the pulse shaping choppers and the detectors.
BIFROST_BANDWIDTH = 1.7715  # angstrom


def bifrost(energy_min=None, wavelength_max=None, shaping_time=2e-4):
    """Chopper settings for BIFROST, from the slow edge of the band you want.

    BIFROST passes a fixed bandwidth of about 1.77 angstrom, so one number decides where
    that band sits. Give it as either the **longest wavelength** or, equivalently, the
    **lowest energy** you want to reach the sample; the band runs from there to about
    1.77 angstrom shorter::

        wavelength_max=3.0  ->  1.21 to  2.98 angstrom  (   9.2 to  55.7 meV)
        energy_min=7.0      ->  1.63 to  3.40 angstrom  (   7.1 to  30.8 meV)

    Note the sense of it: the number you give is the *slow* end of the band. Asking for
    ``energy_min=7`` does not cap the incident energy at 7 meV -- it puts 7 meV at the
    bottom and lets everything faster through, up to about 31 meV.

    Parameters
    ----------
    energy_min : float, optional
        The lowest incident energy to pass, in meV.
    wavelength_max : float, optional
        The longest incident wavelength to pass, in angstrom. Overrides ``energy_min``.
    shaping_time : float, optional
        How long the pulse shaping pair should burst for, in seconds. They co-rotate at
        a harmonic of 14 Hz, one lagging the other, so the burst is one slit crossing
        less the lag; asking for longer than a crossing at 196 Hz drops them to a lower
        harmonic and says so. Defaults to 0.2 ms.

    Returns
    -------
    chopcal.lib.ChopperSet
        The six choppers by name, in beam order -- ``ps1`` and ``ps2`` shaping the pulse,
        ``fo1`` and ``fo2`` suppressing frame overlap, ``bw1`` and ``bw2`` cutting the
        band -- each a :class:`~chopcal.lib.Chopper` of ``speed``, ``delay``, ``angle``
        and ``path``. Feed ``.values()`` to :func:`chopcal.lib.wavelength_limits` to see
        what the train actually passes.

    Examples
    --------
    >>> settings = bifrost(wavelength_max=3.0)                  # doctest: +SKIP
    >>> settings                                                # doctest: +SKIP
    name  speed [Hz]  delay [ms]  opening [deg]  open [ms]  path [m]
     ps1         196     3.87771            170      2.409     6.342
     ps2         196       6.087            170      2.409     6.362
     fo1          14     6.13892          38.26      7.591      8.53
     fo2          14     9.54467          52.01      10.32    14.973
     bw1          14     42.8605            161      31.94        78
     bw2         -14     42.8605            161      31.94     78.02
    >>> settings['bw1'].delay                                   # doctest: +SKIP
    0.042860481...
    """
    if not energy_min and not wavelength_max:
        raise ValueError(
            'Give the slow edge of the band you want, as either energy_min in meV or '
            'wavelength_max in angstrom'
        )
    # The compiled function takes both and prefers a positive wavelength, so a falsy
    # value for either means 'not given' -- which is what makes the positional call
    # bifrost(7.0, 0) work.
    return lib.ChopperSet(
        sorted(_bifrost(energy_min or 0.0, wavelength_max or 0.0, shaping_time).items(),
               key=lambda pair: pair[1].path)
    )


__all__ = [
    "__version__",
    "BIFROST_BANDWIDTH",
    "bifrost",
    "lib"
]
