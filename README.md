# IdentifyFibers RC2 · 2D / 3D

Private test candidate for CellProfiler 4.2.8. The module supports **2D image** and **3D volume (test candidate)**. Pipelines from v1.0.1 load in 2D mode. The volumetric method uses Sato enhancement, centerline graphs and FIRE-inspired continuation.

[Manual em português](README.pt-BR.md) · [Five audits](AUDIT.md) · [2D manual](docs/manual.en.md)

## Testing the 2D update

Choose **2D image**, open the same pipeline and use the same acquisition and settings for both candidates. Time IdentifyFibers with the display window closed. Compare the fiber labels, centerlines and exported measurements. The included benchmark reproduces the comparison on microscopy data. [Timing and validation](PERFORMANCE.md)

## Run the supplied test

1. Copy `active_plugins/identifyfibers.py` and the whole `active_plugins/identifyfibers_engine/` folder into a separate plugins directory. Select that directory in CellProfiler Preferences and restart the application. Keep your usual plugins directory for switching back.
2. Open `demo/Fiber3D_demo.cppipe`. Add **only** `demo/fiber_phantom.ome.tif` to Images. Choose an output directory. Use Test Mode to inspect Z planes, then run Analyze Images to export CSVs.
3. The phantom has three fibers, including two that cross in the XY projection at different depths. The expected path lengths, sorted by length, are approximately **19.25, 20.50 and 23.69 µm**. The numerical reference is in `demo/expected_fibers.csv`.

![Test volume in XY and XZ](demo/phantom_views.png)

## Calibrate and inspect a confocal stack

Use a grayscale Z-stack for one channel and one time point. Enable **Process as 3D** in NamesAndTypes and enter relative X, Y and Z spacing. For 0.25 µm XY pixels and a 0.80 µm Z interval, use X = 1, Y = 1 and Z = 3.2. In IdentifyFibers choose **3D volume (test candidate)** and enter **0.25** for Y pixel size (µm). Zero marks missing calibration. Confirm all three exported `Fiber3D_VoxelSize…_um` values in the image CSV. Set the spacing for each acquisition separately.

| Setting | Interpretation |
|---|---|
| Minimum / maximum expected width (µm) | Choose the apparent widths you see in the original data. Decreasing the minimum can recover thin fibers and enhance noise. Increasing the maximum can favor thicker bundles and lengthens the neighborhood used to estimate branch direction. |
| Threshold correction factor | Start at 1. Reduce it when visible weak segments disappear, or increase it to reject background. Inspect several planes for unintended bridges. |
| Minimum traced length (µm) | A larger value removes shorter paths, including genuine short branches. The cutoff is applied after tracing. |
| Maximum crossing deflection | Smaller angles favor straight continuation. Larger angles admit tighter curves and can switch identity at a junction. |
| CellProfiler binary volume | Supply a registered binary image matching the input shape and spacing to define the foreground. The original grayscale input supplies intensity measurements. |
| Within objects | Extract fibers independently inside each parent. A fiber crossing two parents is measured in pieces. |

## How the volume becomes a set of objects

1. **Enhancement.** Valid intensities are scaled within the analysis region. The multiscale Sato tubular response uses a physical-coordinate Hessian and the geometric mean of the magnitudes of its two negative transverse eigenvalues. Otsu thresholds positive responses. When noise is estimable, a local contrast check on the original signal rejects trajectories with insufficient median support. External masks bypass the detection and contrast stages.
2. **Centerlines.** Foreground uses 26-connectivity. Already thin digital lines are preserved. Lee thinning produces a graph in thicker structures. The nearest continuation to a straight line is matched first at junctions, with directions measured in physical space. Every branch is consumed at most once. When a whole component thins away, a distance-guided geodesic can recover an elongated tube. The fallback requires at least a 2:1 ratio between the two largest principal standard deviations. Its objects carry a flag.
3. **Length.** Physical arc-length resampling and smoothing at the scale of one traversed voxel step reduce the bias of digital staircases. The endpoints of open paths are fixed. Both the resulting curve length and the unsmoothed voxel-step length are exported.
4. **Voxel ownership.** Accepted centerlines seed a geodesic partition within each foreground component. Physical step length and changes in input intensity guide the boundary. Each voxel in a seeded component has one owner. A short rejected spur attached to an accepted component can still contribute to its assigned volume. Unseeded components remain background.

## Measurements

All fiber features listed here carry the `Fiber3D_` prefix.

| Columns | Definition |
|---|---|
| `PathLength_um`, `VoxelPathLength_um`, `EndToEnd_um`, `Straightness` | Smoothed physical path length, raw digital path length, endpoint distance and endpoint distance divided by smoothed length. |
| `Volume_um3`, `ApparentWidth_um` | Assigned voxel count times voxel volume; mean of twice the physical distance from the traced axis to background, with the exterior treated as background. Junctions, unresolved bundles and optical blur affect apparent width. |
| `Azimuth_deg`, `Inclination_deg` | Principal trajectory axis. XY azimuth is modulo 180° (X right, Y increasing with rows). Inclination is signed relative to XY using that axial convention. A vertical fiber has undefined azimuth and 90° inclination. Both angles are undefined for a degenerate principal axis. |
| `IntensityMean`, `IntensityIntegrated` | Mean and sum of the input values loaded by CellProfiler across the region. The sum uses voxel values directly. |
| `TouchesImageEdge`, `FallbackTrace` | Flags for contact with a volume boundary and use of the alternate geodesic trace. |
| `Location_Center_X/Y/Z` (standard CellProfiler names) | Region centroid in voxel indices. |

Start with a crop containing an isolated fiber, a dense bundle and a projected crossing. Inspect the raw signal through Z. An unresolved bundle may have several possible biological interpretations. Fine axial details require adequate optical resolution and Z sampling. The FFT is a feature of the 2D workflow. The 3D controls reported here are synthetic; validation on biological stacks remains open.

**Methodological reference.** Stein AM, Vader DA, Jawerth LM, Weitz DA, Sander LM (2008). An algorithm for extracting the network geometry of three-dimensional collagen gels. *Journal of Microscopy* 232, 463–475. [DOI](https://doi.org/10.1111/j.1365-2818.2008.02141.x). This Python implementation combines FIRE-inspired continuation with Sato enhancement and graph-based extraction. [Attribution](NOTICE.txt).
