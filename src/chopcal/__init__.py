from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("chopcal")
except PackageNotFoundError:
    # package is not installed
    pass

del version
del PackageNotFoundError

from chopcal._chopcal_impl import bifrost
import chopcal.lib as lib

__all__ = [
    "__version__",
    "bifrost",
    "lib"
]
