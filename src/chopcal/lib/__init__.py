from chopcal._chopper_lib_impl import Chopper, inverse_velocity_windows, inverse_velocity_limits, wavelength_limits


def _chopper_repr(chopper):
    """Unambiguous, and in the units the fields are actually stored in."""
    return (f"Chopper(speed={chopper.speed!r}, delay={chopper.delay!r}, "
            f"angle={chopper.angle!r}, path={chopper.path!r})")


def _chopper_str(chopper):
    """The same four numbers, said out loud."""
    return (f"{abs(chopper.speed):g} Hz "
            f"{'clockwise' if chopper.speed < 0 else 'anticlockwise'}, "
            f"{chopper.angle:g} deg opening on the beam at {chopper.delay * 1e3:.4g} ms, "
            f"{chopper.path:g} m from the source")


def _scipp():
    """scipp, if it is installed. It is not a dependency, and does not become one."""
    try:
        import scipp
    except ImportError as error:
        raise ImportError(
            'chopper quantities need scipp: pip install scipp. The plain attributes '
            'need nothing -- speed is Hz, delay seconds, angle degrees, path metres.'
        ) from error
    return scipp


_UNITS = (('speed', 'Hz'), ('delay', 's'), ('angle', 'deg'), ('path', 'm'))


def _chopper_quantities(chopper):
    """The four fields as scipp scalars, each carrying its own unit.

    The attributes are plain numbers in fixed units -- Hz, seconds, degrees, metres --
    and this is the same four values with the units attached, so that whatever reads them
    can convert rather than assume. Worth reaching for whenever the numbers are going
    somewhere else: the table above prints delays in milliseconds because that is the
    scale they live on, while ``chopper.delay`` is in seconds, and a scalar cannot be
    read the wrong way round.
    """
    sc = _scipp()
    return {field: sc.scalar(float(getattr(chopper, field)), unit=unit)
            for field, unit in _UNITS}


Chopper.__repr__ = _chopper_repr
Chopper.__str__ = _chopper_str
Chopper.quantities = property(_chopper_quantities)


_COLUMNS = (
    ('name', '{}', lambda name, c: name),
    ('speed [Hz]', '{:.6g}', lambda name, c: c.speed),
    ('delay [ms]', '{:.6g}', lambda name, c: c.delay * 1e3),
    ('opening [deg]', '{:.6g}', lambda name, c: c.angle),
    ('open [ms]', '{:.4g}', lambda name, c: c.angle / 360 / abs(c.speed) * 1e3
                            if c.speed else float('inf')),
    ('path [m]', '{:.6g}', lambda name, c: c.path),
)


class ChopperSet(dict):
    """The choppers of an instrument, by name, in beam order.

    A ``dict`` in every respect -- ``set['ps1']``, ``set.values()``, ``**set`` all behave
    as usual -- that prints itself as a table rather than as six memory addresses.

    ``delay`` is when an opening is on the beam and ``open`` how long it stays there, so
    a chopper passes neutrons from ``delay - open/2`` to ``delay + open/2``, and again
    every ``1/speed`` after that. Both are shown in milliseconds because that is the
    scale they live on; the attributes themselves are in seconds.

    Use :attr:`quantities` to get them as scipp scalars instead, which is the safer thing
    to hand to anything else -- the units come along and cannot be misread.
    """

    @property
    def quantities(self):
        """Every chopper's fields as scipp scalars, by name.

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


__all__ = [
    'Chopper',
    'ChopperSet',
    'inverse_velocity_windows',
    'inverse_velocity_limits',
    'wavelength_limits',
]
