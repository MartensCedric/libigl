# Python binding: `winding_number_one_shot_2d`

A minimal [nanobind](https://github.com/wjakob/nanobind) module exposing
`igl::winding_number_one_shot_2d` (the One-Shot 2D generalized winding number)
to Python. It follows the `libigl/libigl-python-bindings` convention: each
`src/<name>.cpp` defines `bind_<name>(nb::module_&)`, and CMake auto-wires it
into `module.cpp`. Adding more bindings later is just dropping in more files.

Unlike the upstream bindings repo (which fetches a *pinned* libigl release),
this module builds against **the parent checkout** of this repo via
`add_subdirectory(..)`, so it picks up `winding_number_one_shot_2d` directly.

## Build

```bash
cd python
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

The compiled extension lands in `python/igl/`.

## Use

```python
import sys; sys.path.insert(0, "python")   # or run from the python/ dir
import igl
import numpy as np

# 4-Bézier approximation of a unit circle (CCW), 4 rows (P0..P3) per curve.
k = 0.5522847498307936
C = np.array([
    [ 1, 0], [ 1, k], [ k, 1], [ 0, 1],   # arc Q1: (1,0)->(0,1)
    [ 0, 1], [-k, 1], [-1, k], [-1, 0],   # arc Q2: (0,1)->(-1,0)
    [-1, 0], [-1,-k], [-k,-1], [ 0,-1],   # arc Q3: (-1,0)->(0,-1)
    [ 0,-1], [ k,-1], [ 1,-k], [ 1, 0],   # arc Q4: (0,-1)->(1,0)
], dtype=float)

Q = np.array([[0.0, 0.0], [5.0, 5.0]])    # inside -> ~1, outside -> ~0
print(igl.winding_number_one_shot_2d(C, Q))      # [1., 0.]
print(igl.winding_number_one_shot_2d(C, Q[0]))   # 1.0  (single-point overload)
```

## Caveat

The package dir is named `igl/` to match the upstream convention, so a
`python/` on `sys.path` will **shadow** a pip-installed `igl`. This is a
single-function experimental module; keep it off the global path if you also
use the full `igl` PyPI package.
