from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("chopcal")
except PackageNotFoundError:
    # package is not installed
    pass

from chopcal._chopcal_impl import bifrost
import chopcal.lib as lib

__all__ = [
    "bifrost",
    "lib"
]
