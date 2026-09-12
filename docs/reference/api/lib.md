# API: `chopcal.lib`

Everything chopper-lib exposes: discs, the bands they pass, the grid mask, and the exact
transmitted region.

The classes here are mostly nanobind bindings onto the C, so they are the same objects a
McStas instrument works with rather than a Python reimplementation of them.

## Describing a disc

::: chopcal.lib.Chopper

::: chopcal.lib.beam_aperture

::: chopcal.lib.ChopperSet

## What a train passes

::: chopcal.lib.wavelength_windows

::: chopcal.lib.wavelength_limits

::: chopcal.lib.inverse_velocity_windows

::: chopcal.lib.inverse_velocity_limits

## The exact transmitted region

::: chopcal.lib.Region

::: chopcal.lib.Polygon

::: chopcal.lib.RegionSampler

## On a grid

::: chopcal.lib.inverse_velocity_time_mask

::: chopcal.lib.unmasked_probability

::: chopcal.lib.MaskSampler

::: chopcal.lib.MaskValue

## Conversions

::: chopcal.lib.wavelength_to_inverse_velocity

::: chopcal.lib.inverse_velocity_to_wavelength
