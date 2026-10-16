"""BIFROST chopper settings as NICOS sets them.

An ESS chopper controller does not take chopcal's ``delay``. It locks its top-dead-centre
sensor to a total, which the EPICS module ``evr_seq_calc`` forms from three parts, in ns::

    TotDly = MechDly / 360 / |Spd_S| * 1e9 + BeamPosDly + ChopDly

``MechDly-S`` is measured for each disc, so that with the other two at zero its first
opening edge -- NeXus's ``slit_edges[0]`` -- is at the beam at the timing reference.
``BeamPosDly-S`` (the pulse delay) moves every disc together and is fixed by the
instrument team. ``ChopDly-S`` (the experiment delay) is what an experiment sets, and the
only delay NICOS writes, through ``<disc>_chopper_delay``. So the first opening edge is at
the beam ``pulse_delay + experiment_delay`` after the reference, and nothing about the real
disc's mounting enters what an experiment has to set.

chopcal's ``delay`` is when the disc point at ``beam`` is on the beam, which for every
chopper ``chopcal.bifrost`` describes is the centre of its one opening. Both are measured
from the proton pulse, so the experiment delay follows from chopcal's own geometry alone:
the time ``slit_edges[0]`` reaches the beam, less the pulse delay.
"""
from . import bifrost as _bifrost

#: chopcal's names for BIFROST's discs, and the prefix of each one's NICOS devices.
NICOS_NAMES = {'ps1': 'psc1', 'ps2': 'psc2', 'fo1': 'foc1', 'fo2': 'foc2',
               'bw1': 'bwc1', 'bw2': 'bwc2'}


def experiment_delay(chopper, pulse_delay=0.0):
    """ESS ``ChopDly-S`` for one chopper, in seconds, in ``[0, 1/|speed|)``.

    ``pulse_delay`` is the instrument's ``BeamPosDly-S``, in seconds. An edge at angle
    ``a`` is on the beam at ``delay + (beam - a) / (360 speed)``, as chopper-lib defines
    it, and the openings recur every turn, so the result is folded into one turn of the
    disc. A disc that is not turning has no delay; zero stands in for it.

    Only a disc of one opening has a first opening edge that does not depend on where its
    top-dead-centre mark is, which chopcal does not know; any other is refused.
    """
    if len(chopper.edges) != 2:
        raise ValueError(
            f'{len(chopper.edges) // 2} openings: which is slit_edges[0] depends on where '
            f'the disc\'s top-dead-centre mark is, which chopcal does not know')
    if not chopper.speed:
        return 0.0
    first = chopper.edges[0]
    at_beam = chopper.delay + (chopper.beam - first) / (360.0 * chopper.speed)
    return (at_beam - pulse_delay) % (1.0 / abs(chopper.speed))


def bifrost(pulse_delay, energy_min=None, wavelength_max=None, shaping_time=2e-4,
            apertures=True):
    """BIFROST's chopper settings, keyed by the NICOS device that takes each one.

    The arguments after ``pulse_delay`` are :func:`chopcal.bifrost`'s. ``pulse_delay`` is
    the ``BeamPosDly-S`` the chopper controllers hold, in **ns** as they hold it; it has
    no default because no value is right for every installation. Each disc gets its
    ``<disc>_chopper_speed`` in Hz, signed -- ``bwc2`` counter-rotates -- and its
    ``<disc>_chopper_delay``, the experiment delay ``ChopDly-S``, in ns::

        >>> bifrost(pulse_delay=10_000, energy_min=7.0)     # doctest: +SKIP
        {'psc1_chopper_speed': 196.0, 'psc1_chopper_delay': ..., ...}

    This takes chopcal's timing reference, the proton pulse, to be the one the choppers
    and the detectors are timed from (NiSync), and ``slit_edges[0]`` to be the edge
    ``MechDly-S`` is measured to.
    """
    settings = _bifrost(energy_min=energy_min, wavelength_max=wavelength_max,
                        shaping_time=shaping_time, apertures=apertures)
    values = {}
    for key, name in NICOS_NAMES.items():
        chopper = settings[key]
        values[f'{name}_chopper_speed'] = chopper.speed
        values[f'{name}_chopper_delay'] = 1e9 * experiment_delay(chopper,
                                                                 1e-9 * pulse_delay)
    return values


__all__ = ['NICOS_NAMES', 'bifrost', 'experiment_delay']
