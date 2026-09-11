from enum import IntEnum

from chopcal._chopper_lib_impl import (
    Chopper,
    inverse_velocity_windows,
    inverse_velocity_limits,
    wavelength_limits,
    MASK_EXCLUDED,
    MASK_INCLUDED,
    MASK_GROWN,
)
from chopcal._chopper_lib_impl import MaskSampler as _MaskSampler
from chopcal._chopper_lib_impl import (
    inverse_velocity_time_mask as _inverse_velocity_time_mask,
    unmasked_probability as _unmasked_probability,
)


class MaskValue(IntEnum):
    """What a cell of an `inverse_velocity_time_mask` holds.

    A finished mask only ever holds EXCLUDED or INCLUDED -- GROWN is the marker the
    growing pass leaves behind and folds back into INCLUDED before it returns -- but it
    is written to file as it stands, so a mask read back from one can carry it.
    """

    EXCLUDED = MASK_EXCLUDED
    INCLUDED = MASK_INCLUDED
    GROWN = MASK_GROWN


def _chopper_repr(chopper):
    """Unambiguous, and in the units the fields are actually stored in."""
    return (f"Chopper(speed={chopper.speed!r}, delay={chopper.delay!r}, "
            f"beam={chopper.beam!r}, edges={list(chopper.edges)!r}, "
            f"path={chopper.path!r}, aperture={chopper.aperture!r})")


def _openings(chopper):
    """The (opening, closing) edge pairs, which is how `edges` is meant to be read."""
    return list(zip(chopper.edges[::2], chopper.edges[1::2]))


def _chopper_str(chopper):
    """The same disk, said out loud."""
    count = len(chopper.edges) // 2
    if count == 1:
        low, high = _openings(chopper)[0]
        where = (f"{high - low:g} deg opening centred on the beam"
                 if low == -high else
                 f"{high - low:g} deg opening at {(low + high) / 2:g} deg")
    else:
        where = f"{count} openings totalling {chopper.opening:g} deg"
    wide = f", {chopper.aperture:g} deg of beam" if chopper.aperture else ""
    return (f"{abs(chopper.speed):g} Hz "
            f"{'clockwise' if chopper.speed < 0 else 'anticlockwise'}, "
            f"{where} on the beam at {chopper.delay * 1e3:.4g} ms, "
            f"{chopper.path:g} m from the source{wide}")


Chopper.openings = property(
    _openings,
    doc="The (opening, closing) edge pairs, in degrees from the disk's zero mark.")


def _scipp():
    """scipp, if it is installed. It is not a dependency, and does not become one."""
    try:
        import scipp
    except ImportError as error:
        raise ImportError(
            'chopper quantities need scipp: pip install scipp. The plain attributes '
            'need nothing -- speed is Hz, delay seconds, beam, edges and aperture '
            'degrees, path metres.'
        ) from error
    return scipp


def _numpy():
    """numpy, likewise. Only the mask functions need it, and only when called."""
    try:
        import numpy
    except ImportError as error:
        raise ImportError(
            'chopper masks are arrays, and need numpy: pip install numpy. Nothing else '
            'in chopcal does -- the window and limit functions return plain tuples.'
        ) from error
    return numpy


_SCALAR_UNITS = (('speed', 'Hz'), ('delay', 's'), ('beam', 'deg'),
                 ('path', 'm'), ('aperture', 'deg'))


def _chopper_quantities(chopper):
    """The fields as scipp variables, each carrying its own unit.

    The attributes are plain numbers in fixed units -- Hz, seconds, degrees, metres --
    and this is the same values with the units attached, so that whatever reads them
    can convert rather than assume. Worth reaching for whenever the numbers are going
    somewhere else: the table above prints delays in milliseconds because that is the
    scale they live on, while ``chopper.delay`` is in seconds, and a scalar cannot be
    read the wrong way round.

    ``edges`` is a list rather than a single number, so it comes back as an array
    variable over an ``edge`` dimension.
    """
    sc = _scipp()
    out = {field: sc.scalar(float(getattr(chopper, field)), unit=unit)
           for field, unit in _SCALAR_UNITS}
    out['edges'] = sc.array(dims=['edge'], values=[float(e) for e in chopper.edges],
                            unit='deg')
    return out


Chopper.__repr__ = _chopper_repr
Chopper.__str__ = _chopper_str
Chopper.quantities = property(_chopper_quantities)


def _opening_column(chopper):
    """One opening reads as its width; several read as a count and a total."""
    count = len(chopper.edges) // 2
    if count == 1:
        return f'{chopper.opening:.6g}'
    return f'{chopper.opening:.6g} ({count})'


def _open_ms(chopper):
    """How long the beam spends inside the openings each turn, in milliseconds.

    The total opening angle over the turn rate. It was a single ``angle`` before
    chopper-lib 4.0.0, which could not describe a disk of more than one opening.
    """
    if not chopper.speed:
        return float('inf')
    return chopper.opening / 360 / abs(chopper.speed) * 1e3


_COLUMNS = (
    ('name', '{}', lambda name, c: name),
    ('speed [Hz]', '{:.6g}', lambda name, c: c.speed),
    ('delay [ms]', '{:.6g}', lambda name, c: c.delay * 1e3),
    ('beam [deg]', '{:.6g}', lambda name, c: c.beam),
    ('opening [deg]', '{}', lambda name, c: _opening_column(c)),
    ('open [ms]', '{:.4g}', lambda name, c: _open_ms(c)),
    ('path [m]', '{:.6g}', lambda name, c: c.path),
)


class ChopperSet(dict):
    """The choppers of an instrument, by name, in beam order.

    A ``dict`` in every respect -- ``set['ps1']``, ``set.values()``, ``**set`` all behave
    as usual -- that prints itself as a table rather than as six memory addresses.

    ``delay`` is when the disk point at ``beam`` is on the beam path and ``open`` how
    long the openings spend there each turn, so a single-opening disk centred on the beam
    passes neutrons from ``delay - open/2`` to ``delay + open/2``, and again every
    ``1/speed`` after that. Both are shown in milliseconds because that is the scale they
    live on; the attributes themselves are in seconds.

    Use :attr:`quantities` to get them as scipp variables instead, which is the safer
    thing to hand to anything else -- the units come along and cannot be misread.
    """

    @property
    def quantities(self):
        """Every chopper's fields as scipp variables, by name.

            >>> settings.quantities['ps1']['delay']     # doctest: +SKIP
            <scipp.Variable> ()  float64  [s]  0.00595313
        """
        return {name: chopper.quantities for name, chopper in self.items()}

    def _rows(self):
        return [[fmt.format(get(name, chopper)) for _, fmt, get in _COLUMNS]
                for name, chopper in self.items()]

    def __repr__(self):
        if not self:
            return f'{type(self).__name__}()'
        headers = [head for head, _, _ in _COLUMNS]
        rows = self._rows()
        widths = [max(len(head), *(len(row[i]) for row in rows))
                  for i, head in enumerate(headers)]
        lines = ['  '.join(h.rjust(w) for h, w in zip(headers, widths))]
        lines += ['  '.join(v.rjust(w) for v, w in zip(row, widths)) for row in rows]
        return '\n'.join(lines)

    def _repr_pretty_(self, printer, cycle):
        """IPython pretty-prints dict subclasses as dicts unless asked otherwise."""
        printer.text('...' if cycle else repr(self))

    def _repr_html_(self):
        headers = ''.join(f'<th style="text-align:right">{h}</th>'
                          for h, _, _ in _COLUMNS)
        body = ''.join(
            '<tr>' + ''.join(f'<td style="text-align:right">{v}</td>' for v in row) + '</tr>'
            for row in self._rows()
        )
        return f'<table><thead><tr>{headers}</tr></thead><tbody>{body}</tbody></table>'


def _as_edges(values):
    """Bin edges as the compiled side wants them: contiguous C doubles.

    The bindings take typed arrays rather than sequences, so that the C reads the
    caller's buffer instead of a copy. Coercing here means a list works anyway, and that
    a caller passing float32 or a sliced view gets a conversion rather than a TypeError
    naming a signature they did not write.
    """
    return _numpy().ascontiguousarray(values, dtype=float)


def _as_mask(values):
    """A mask as the compiled side wants it: contiguous C ints.

    `inverse_velocity_time_mask` returns exactly this, so the common path is a no-op;
    the coercion is for a mask read back from a file, or one built by hand, which on
    most platforms defaults to a width the binding will not accept.
    """
    np = _numpy()
    return np.ascontiguousarray(values, dtype=np.int32)


def inverse_velocity_time_mask(choppers, inverse_velocities, times, grow=0):
    """Which (inverse velocity, time) bins a chopper train passes.

    `inverse_velocities` (s/m) and `times` (s, at the source) are bin *edges*, so the
    mask is one smaller in each direction. It comes back shaped
    **(time, inverse velocity)** -- time is the slow axis, which is how the C stores it
    and the opposite of the argument order.

    `grow` expands each allowed region by that many bins in every direction. It is a
    blunt instrument: it admits inverse velocities no disk ever passes. Prefer
    :attr:`Chopper.aperture` for a beam of finite width, which opens the windows in time
    only, where the width actually acts.

    Returns ``(mask, allowed_bin_count)``. Needs numpy.
    """
    return _inverse_velocity_time_mask(choppers, _as_edges(inverse_velocities),
                                       _as_edges(times), grow)


def unmasked_probability(signal, mask):
    """The fraction of `signal` lying in the allowed bins of `mask`.

    Both are shaped (time, inverse velocity), as `inverse_velocity_time_mask` returns.
    This is the transmission an instrument sees -- a *weighted* fraction -- and is not
    the weight correction a sampler needs; that is :attr:`MaskSampler.acceptance`, which
    counts draws rather than intensity.
    """
    np = _numpy()
    return _unmasked_probability(np.ascontiguousarray(signal, dtype=float),
                                 _as_mask(mask))


class MaskSampler(_MaskSampler):
    """Draws (inverse velocity, time) pairs from the allowed cells of a finished mask.

    Sampling the whole region and discarding what the mask excludes spends the ray budget
    to keep the fraction the choppers pass, which for a real chopper train is a fraction
    of a percent. Drawing from the allowed cells keeps all of it and is the same
    distribution, provided every drawn weight is multiplied by :attr:`acceptance`.

    `mask` is shaped (time, inverse velocity) and `inverse_velocities`/`times` are its
    bin edges. The four bounds describe the region the caller samples uniformly, which is
    not in general the region the grid covers -- a grid sized with ``ceil`` runs past it
    in the last row and column, and those cells are clipped before they are weighted so
    that they do not inflate :attr:`acceptance`.
    """

    def __init__(self, mask, inverse_velocities, times, inverse_velocity_minimum,
                 inverse_velocity_range, time_minimum, time_range):
        super().__init__(_as_mask(mask), _as_edges(inverse_velocities),
                         _as_edges(times), inverse_velocity_minimum,
                         inverse_velocity_range, time_minimum, time_range)


def _sample(sampler, count, generator=None):
    """``count`` (inverse_velocity, time) pairs, as two numpy arrays.

    The compiled ``draw`` and ``draw_many`` take their uniform deviates as arguments, so
    that chopper-lib needs no generator of its own and a McStas TRACE can hand over
    ``rand01()``. This is the convenience for everything else: it draws the deviates from
    ``generator`` -- ``numpy.random.default_rng()`` if none is given -- and calls
    ``draw_many``.

    The caller still owes every drawn ray's weight a factor of ``sampler.acceptance``.
    """
    np = _numpy()
    if generator is None:
        generator = np.random.default_rng()
    deviates = generator.random((3, count))
    inverse_velocity, time = sampler.draw_many(
        np.ascontiguousarray(deviates[0]), np.ascontiguousarray(deviates[1]),
        np.ascontiguousarray(deviates[2]))
    return np.asarray(inverse_velocity), np.asarray(time)


_MaskSampler.sample = _sample


__all__ = [
    'Chopper',
    'ChopperSet',
    'MaskSampler',
    'MaskValue',
    'inverse_velocity_windows',
    'inverse_velocity_limits',
    'wavelength_limits',
    'inverse_velocity_time_mask',
    'unmasked_probability',
]
