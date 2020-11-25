# pyKVFinder information
_name = "pyKVFinder"
_version = "0.1"
_license = "GPL-3.0 License"

try:
    from .utils import *
    from .grid import *
    from .main import *
except SyntaxError:
    pass
