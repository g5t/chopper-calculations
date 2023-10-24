from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("chopcal")
except PackageNotFoundError:
    # package is not installed
    pass

from ._chopcal_impl import bifrost

__all__ = ["bifrost"]
