# Chopcal

## Description
Exposes chopper calculations to Python which were otherwise hidden in McStas instruments.


## Supported components

| source          | name                                                            | component                                 | type      |
|-----------------|-----------------------------------------------------------------|-------------------------------------------|-----------|
| Instr           | BIFROST                                                         | `chopcal.bifrost`                         | function  |
| runtime library | [chopper-lib](https://github.com/mcdotstar/mcstas-chopper-lib/) | `chopcal.mcstas`                          | submodule |
|                 | chopper parameters                                              | `chopper.mcstas.Chopper`                  | class     |
|                 | transmitted inverse-velocity phase space                        | `chopper.mcstas.inverse_velocity_windows` | function  |
|                 | transmitted inverse-velocity extremea                           | `chopper.mcstas.inverse_velocity_limits`  | function  |
|                 | transmitted wavelength extremea                                 | `chopper.mcstas.wavelength_limits`        | function  |



## Installation

```bash
pip install chopcal
```

Or, from the source repository to get the latest development version
```bash
pip install git+https://github.com/g5t/chopper-calculations.git
```
