"""humanSLmapping package. `import hsl` installs the YEAST_SLs write guard (hsl.guard) and configures the vendored
upstream `ysl` package (third_party/ysl) to read data in place from YEAST_SLs (config/paths.yaml)."""
from hsl import paths as _paths  # noqa: F401  (side effect: guard + ysl path env)
