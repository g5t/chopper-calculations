# Build the transmitted region

`Region` is the region a chopper train transmits, exactly — a union of convex polygons in
(inverse velocity, emission time).

## Build one

Start from the rectangle your source draws from, and transmit it through the train:

```python
--8<-- "transmitted_region.py:build"
```

```
1 polygon(s)
[(1.1458362867034736, 3.0515432206219723)]
acceptance 3.929119e-03
```

`Region.from_wavelengths(low, high, time_minimum, time_range)` takes angstrom and seconds;
`Region.rectangle` takes inverse velocity directly if you have it.

For the six BIFROST discs the answer is **one polygon of five vertices**. That is the whole
transmitted phase space of the instrument.

## Why it is exact

A neutron emitted at inverse velocity `a` and time `t` reaches path `L` at `t + L*a`, so a
disc open on `[lower, upper]` accepts

```
lower <= t + L*a <= upper
```

a slab between two parallel lines. A train's acceptance is an intersection of unions of
such slabs, and intersection distributes over union — so the exact region **is** a union of
convex pieces, one per choice of which opening and which turn of each disc a neutron goes
through. Nothing is approximated and nothing is quantised.

## Read the pieces

```python
--8<-- "transmitted_region.py:polygons"
```

Each `Polygon` has `vertices`, `area` and `contains`. The vertices are the disc edges, at
the precision of a double.

Because the pieces are **disjoint** — two turns of one disc cannot both pass the same
neutron — `region.area` is a plain sum, with no inclusion-exclusion.

## Ask about one neutron

```python
--8<-- "transmitted_region.py:contains"
```

The boundary counts as inside.

## Reuse the source rectangle

`transmit` answers a question; it does not consume the region it was asked of.

```python
--8<-- "transmitted_region.py:reuse"
```

```
1 1 1
```

So one rectangle can be transmitted through several trains, or through prefixes of one, to
see which discs are doing the cutting.

## Write it out

```python
--8<-- "transmitted_region.py:write"
```

`region.to_dict(sampled=source)` gives plain data; `region.write_json(path, sampled=source)`
writes the same thing to a file — byte for byte what the `Polygon_ESS_butterfly` McStas
component writes beside its data. Every number carries enough digits to read back
bit-exact, and the file is a few hundred bytes where a mask of any useful resolution is
megabytes.

`sampled` is what sets the acceptance; without it both it and `acceptance` are `null`.

## Next

* [Sample without wasting rays](sample-the-region.md) — what the region is *for*.
* [Allow for a guide's path spread](path-spread.md).
