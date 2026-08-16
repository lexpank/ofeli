# Developer documentation

Ofeli (Open, Fast and Efficient Level set Implementation) demonstrates image
segmentation with discrete level-set curve evolution (implicit active
contours), based on the work of Y. Shi and W. C. Karl. The approach avoids
solving partial differential equations while retaining the topological-change
handling of level-set methods.

## Structure

The repository separates the Qt front ends from the image-processing core:

- `apps/gui/` contains the Qt Widgets application and `PixmapWidget`.
- `apps/cli/` contains the batch command-line front end.
- `src/core/` contains the numerical algorithms and `SegmentationSession`.
- `src/qt/` contains `QImage` conversion and atomic image output.
- `include/ofeli/` contains the public C++17 core headers.

Both front ends convert input into the shared `ofeli::Image` type and execute
the validated `ofeli::SegmentationSession`. The session owns the buffers
borrowed by the legacy `Filters` and active-contour implementations, applies
preprocessing in the GUI order, selects exactly one contour model, and returns
an owned `SegmentationResult`.

## Reusing the core

The core is installed as the static CMake target `Ofeli::Core`:

```cmake
find_package(Ofeli CONFIG REQUIRED)
add_executable(example main.cpp)
target_link_libraries(example PRIVATE Ofeli::Core)
```

`ofeli::Image` stores tightly packed grayscale or interleaved RGB bytes. The
public session API validates dimensions, masks, kernels, and numeric options
before constructing the legacy algorithms. The legacy classes remain
available through the installed headers for existing integrations, while new
front ends should use `SegmentationSession`.

## Historical references

1. Y. Shi and W. C. Karl, “A real-time algorithm for the approximation of
   level-set based curve evolution,” *IEEE Transactions on Image Processing*,
   17(5), 2008.
2. T. F. Chan and L. A. Vese, “Active contours without edges,” *IEEE
   Transactions on Image Processing*, 10(2), 2001.
3. V. Caselles, R. Kimmel, and G. Sapiro, “Geodesic active contours,”
   *International Journal of Computer Vision*, 22(1), 61–79, 1997.
4. P. Perona and J. Malik, “Scale-space and edge detection using anisotropic
   diffusion,” *IEEE Transactions on Pattern Analysis and Machine
   Intelligence*, 12(17), 1990.

Ofeli is distributed under the CeCILL v2 license; see `LICENSE`, `license.txt`,
and `licence.txt`.
