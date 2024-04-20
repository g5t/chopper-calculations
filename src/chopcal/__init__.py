from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("chopcal")
except PackageNotFoundError:
    # package is not installed
    pass

from chopcal._chopcal_impl import bifrost
import chopcal._chopper_lib_impl as mcstas

__all__ = [
    "bifrost",
    "mcstas"
]
