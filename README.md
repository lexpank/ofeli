# Ofeli

Ofeli (Open, Fast and Efficient Level set Implementation) is a Qt desktop application for image filtering and active-contour segmentation.

![Ofeli main window](screenshots/ofeli-main-window.png)

![Ofeli segmentation view in French](screenshots/ofeli-segmentation-fr.png)

## Requirements

- CMake 3.16 or newer
- Qt 6.5 or newer with Widgets, PrintSupport, and LinguistTools for the GUI
- A C++17 compiler

The project is configured for Linux, macOS, and Windows. CI currently covers
Linux core-only, full Qt 6.5, sanitizer, CLI, and install/consumer checks;
Windows and macOS remain local release checks.

## Build and test

~~~sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
~~~

`OFELI_BUILD_GUI` and `OFELI_BUILD_CLI` both default to `ON`. Turn either off to omit that front end.

For a Qt-free core build:

~~~sh
cmake -S . -B build-core \
  -DOFELI_BUILD_GUI=OFF -DOFELI_BUILD_CLI=OFF -DBUILD_TESTING=ON
cmake --build build-core --parallel
ctest --test-dir build-core --output-on-failure
~~~

If Qt is installed outside CMake's default search path, add `-DCMAKE_PREFIX_PATH=/path/to/Qt/6.5` to the GUI configure command.

## Desktop GUI

Build with `OFELI_BUILD_GUI=ON`, launch `Ofeli`, and use **Open** to load an
image. The **Settings** dialog selects the model, initialization, optional
preprocessing, colors, and display scale. **Start** runs the segmentation;
**Stop** cancels the current session. Use **Save** for the displayed result
or the preprocessing preview. The release checklist below covers the GUI
interactions that require a real desktop session.

## Command-line segmentation

Build the CLI without the desktop application:

~~~sh
cmake -S . -B build-cli \
  -DOFELI_BUILD_GUI=OFF -DOFELI_BUILD_CLI=ON -DBUILD_TESTING=ON
cmake --build build-cli --parallel
ctest --test-dir build-cli --output-on-failure
~~~

Run a segmentation with the v1 syntax:

~~~text
ofeli-cli segment INPUT --mask OUTPUT [--overlay OUTPUT]
    [--model chan-vese|geodesic] [--initial-mask INPUT] [--force]
~~~

Examples:

~~~sh
build-cli/ofeli-cli segment input.png --mask mask.png
build-cli/ofeli-cli segment input.png --mask mask.png \
  --overlay overlay.png --model geodesic --initial-mask initial.png --force
~~~

The CLI reads any format provided by the local Qt image plugins. Standard Qt
deployments commonly provide BMP, GIF, JPEG/JPG, PNG, PBM/PGM/PPM, TIFF, XBM,
and XPM; the exact set depends on the deployed plugins. Grayscale input is
normalized to 8-bit grayscale and other input to interleaved RGB888. Output
formats are selected from the filename suffix. `--mask` is required;
`--overlay`, `--model`, `--initial-mask`, and `--force` are optional. The CLI
refuses to overwrite outputs unless `--force` is present.

Exit codes are `0` for success/help/version, `2` for usage or output-path
errors, and `1` for image I/O or processing errors.

The French translation is compiled from `apps/gui/Ofeli_fr.ts` during the build; no generated `.qm` file is required in the source tree.

## Install and package

Install the built products, headers, and the `Ofeli::Core` CMake package:

~~~sh
cmake --install build --prefix /tmp/ofeli-install
/tmp/ofeli-install/bin/ofeli-cli --version
~~~

A consumer can use the installed Qt-free core with:

~~~cmake
find_package(Ofeli CONFIG REQUIRED)
add_executable(example main.cpp)
target_link_libraries(example PRIVATE Ofeli::Core)
~~~

Pass the install prefix through `CMAKE_PREFIX_PATH` when configuring that consumer.

## Qt Creator

Open the repository's `CMakeLists.txt` in Qt Creator and select a Qt 6.5+ CMake kit. Configure a fresh build directory, then build the `Ofeli` target. On an Apple Silicon Homebrew installation, set CMake's `CMAKE_PREFIX_PATH` to `/opt/homebrew/opt/qt` if Qt Creator does not find Qt automatically.

If Qt Creator reports `No rule to make target 'all'`, it is using an old qmake/Makefile build configuration. Remove that build configuration in Qt Creator and re-open `CMakeLists.txt`; do not open the old `.pro` project.

## Release smoke checklist

- Launch the app and restore/save window settings.
- Open grayscale and RGB images; run Chan-Vese and geodesic segmentation.
- Apply Gaussian, salt-and-pepper, and speckle noise.
- Switch English/French and restart; open the print dialog and preview/output.
- Check 100%, 150%, and 175% display scaling.

## License

Ofeli is distributed under the CeCILL v2 license. The complete English and French license texts are included in `license.txt` and `licence.txt`; source files retain the required CeCILL v2 notices.
