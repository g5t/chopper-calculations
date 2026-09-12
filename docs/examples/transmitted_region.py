"""Building the region a train transmits, and reading it."""
import chopcal

train = list(chopcal.bifrost(wavelength_max=3.0).values())

# --8<-- [start:build]
from chopcal.lib import Region

# the rectangle the source draws from: wavelengths in AA, then emission time in s
source = Region.from_wavelengths(0.4, 45.0, 0.0, 3e-3)
region = source.transmit(train)

print(len(region), 'polygon(s)')
print(region.wavelength_ranges())
print(f'acceptance {region.area / source.area:.6e}')
# --8<-- [end:build]

# --8<-- [start:polygons]
for polygon in region.polygons:
    print(len(polygon), 'vertices, area', polygon.area)
    for inverse_velocity, time in polygon.vertices:
        print(f'  {inverse_velocity:.6e} s/m at {time * 1e3:.4f} ms')
# --8<-- [end:polygons]

# --8<-- [start:contains]
centre = region.polygons[0].vertices[0]
print(region.contains(*centre))            # a vertex is on the boundary, so inside
print(region.contains(5e-3, 1e-3))         # far outside
# --8<-- [end:contains]

# --8<-- [start:reuse]
# `transmit` answers a question; it does not consume the rectangle
without_bandwidth = source.transmit(train[:4])
print(len(source), len(region), len(without_bandwidth))
# --8<-- [end:reuse]

# --8<-- [start:write]
data = region.to_dict(sampled=source)
print(sorted(data))
print(data['acceptance'])
# region.write_json('region.json', sampled=source)   # the same thing, on disk
# --8<-- [end:write]


def main():
    assert len(region) == 1
    assert region.area > 0
    assert region.area < source.area
    # disjoint pieces, so the areas add
    assert abs(sum(p.area for p in region.polygons) - region.area) < 1e-24
    assert region.contains(*centre)
    assert not region.contains(5e-3, 1e-3)

    # the source is unchanged, and fewer disks pass more
    assert len(source) == 1
    assert without_bandwidth.area > region.area

    assert data['acceptance'] == region.area / source.area
    assert data['chopper_lib_version'].startswith('4.')


main()
