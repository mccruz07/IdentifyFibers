"""IdentifyFibers
==============

Identify and measure bright fibers in 2D grayscale images with Sato +
FIRE-inspired tracing. This independent Python implementation is based on
ideas from FIRE (Fiber Extraction; Stein et al., 2008). It is not the original
FIRE or CT-FIRE software and does not require MATLAB.

Getting started
---------------

1. Set the expected width range to the apparent thickness of your fibers in
   pixels. Width controls enhancement, not fiber length.
2. Start with Otsu and a correction factor of 1. Lower the factor slightly
   if faint segments disappear; increase it if background enters the objects.
3. Set Minimum fiber length to retain biologically relevant short fibers while
   rejecting small fragments. This filter does not reconnect broken fibers.
4. Compare centerlines and colored regions with the original image. Change
   one parameter at a time on representative images before processing a batch.

Reading the result
------------------

Centerlines describe the traced paths; regions contain the pixels assigned
for area and intensity measurements. Each final region pixel has one owner.
This does not guarantee that every traced path is a distinct biological fiber.
Inspect close parallel paths, junctions and bright unresolved bundles.

Geometric coverage is the percentage of the mask skeleton within 2 pixels of
an accepted path. It is a diagnostic, not accuracy. FiberUncoveredSkeleton
shows where coverage is missing.

Use Fiber_CenterlineLength for curved path length. Fiber_Length and the
FiberSummary length measurements retain the region major-axis definition.
Fiber_WidthMean is mask-based and may increase at junctions. FFT measures
orientation in the original image independently of the fiber objects.

Within objects analyzes each parent separately and clips fibers at its border.
An external CellProfiler binary image can define segmentation instead of Sato.
Advanced values remain active when hidden. Use ExportToSpreadsheet for
measurements and SaveImages for output images.

Reference and manual
--------------------

Stein AM, Vader DA, Jawerth LM, Weitz DA, Sander LM (2008).
An algorithm for extracting the network geometry of three-dimensional collagen
gels. Journal of Microscopy 232(3):463-475.
https://doi.org/10.1111/j.1365-2818.2008.02141.x

The GitHub-ready manuals (docs/manual.en.md and docs/manual.pt-BR.md) explain
parameter adjustment, processing stages, measurements and limitations.
"""

import numpy as np
from .geometry import label_overlay

from cellprofiler_core.image import Image
from cellprofiler_core.module import Module
from cellprofiler_core.object import Objects
from cellprofiler_core.constants.measurement import COLTYPE_FLOAT, COLTYPE_INTEGER
from cellprofiler_core.setting import Binary
from cellprofiler_core.setting.choice import Choice
from cellprofiler_core.setting.subscriber import ImageSubscriber, LabelSubscriber
from cellprofiler_core.setting.text import Float, ImageName, Integer, LabelName

from .core import FiberSettings, analyze_fibers, measure_neighbors, centerline_labels
from .fft_orientation import draw_distribution, fft_orientation, distribution_image, report_image


class IdentifyFibers(Module):
    # CellProfiler reads the class docstring, not the package docstring.
    __doc__ = __doc__
    module_name = "IdentifyFibers"
    category = "Object Processing"
    variable_revision_number = 14

    MODE_WHOLE = "Entire image"
    MODE_OBJECTS = "Within objects"

    FEATURES = (
        "Length", "CenterlineLength", "EndToEndLength", "Straightness",
        "Tortuosity", "Orientation", "WidthMean", "WidthMedian", "WidthMaximum",
        "CurvatureMean", "CurvatureMaximum", "IntensityMean", "IntensityMaximum",
        "IntensityIntegrated", "CenterY", "CenterX", "StartY", "StartX", "EndY", "EndX",
        "NearestCentroidDistance", "NeighborAlignment", "Isolated",
    )
    QC_FEATURES = ("CoveragePercent", "ReferenceSkeletonPixels", "UncoveredSkeletonPixels", "CoverageTolerancePixels", "BoundaryDirectionGuidance", "ResolvedSharedTails", "ResolvedSharedSegments")
    SUMMARY_FEATURES = (
        "TotalLength", "MeanLength", "MedianLength", "StdLength",
        "MeanWidth", "MedianWidth", "MaximumWidth", "MeanOrientation",
        "OrientationSD", "AlignmentCoefficient", "MeanStraightness",
        "MeanTortuosity", "MeanCurvature", "MeanIntensity",
        "IsolatedFiberCount", "MeanNearestCentroidDistance", "MeanNeighborAlignment",
    )

    def create_settings(self):
        self.input_image = ImageSubscriber("Select the input fiber image", "None")
        self.analysis_mode = Choice("Analyze fibers in", [self.MODE_WHOLE, self.MODE_OBJECTS])
        self.input_objects = LabelSubscriber("Select parent objects", "None")
        self.output_objects = LabelName("Name the output fiber objects", "Fibers")
        self.enhanced_image = ImageName("Name the enhanced image", "FiberEnhanced")
        self.skeleton_image = ImageName("Name the skeleton image", "FiberSkeleton")
        self.min_width = Float("Minimum expected fiber width (pixels)", 1.0, minval=0.1)
        self.max_width = Float("Maximum expected fiber width (pixels)", 8.0, minval=0.1)
        self.num_scales = Integer("Number of enhancement scales", 4, minval=1)
        self.threshold_method = Choice("Threshold method", ["Otsu", "Li", "Percentile", "Manual", "CellProfiler binary image"])
        self.threshold_value = Float("Manual threshold", 0.15, minval=0.0, maxval=1.0)
        self.threshold_percentile = Float("Threshold percentile", 75.0, minval=0.0, maxval=100.0)
        self.minimum_length = Float("Minimum fiber length (pixels)", 15.0, minval=0.0)
        self.minimum_area = Integer("Minimum segmented area (pixels²)", 10, minval=1)
        self.closing_radius = Integer("Closing radius (pixels)", 1, minval=0)
        self.crossing_max_angle = Float("Maximum crossing deflection (degrees)", 35.0, minval=0.0, maxval=90.0)
        self.tangent_points = Integer("Minimum points for branch direction", 5, minval=2)
        self.threshold_correction = Float("Threshold correction factor", 1.0, minval=0.001)
        self.threshold_lower = Float("Lower threshold bound", 0.0, minval=0.0, maxval=1.0)
        self.threshold_upper = Float("Upper threshold bound", 1.0, minval=0.0, maxval=1.0)
        self.binary_image = ImageSubscriber("Select binary fiber image from CellProfiler", "None")
        self.neighbor_radius = Float("Neighbor search radius (pixels)", 30.0, minval=0.01)
        self.use_fft = Binary("Calculate image-based FFT orientation?", False)
        self.fft_low = Float("FFT minimum frequency (cycles/pixel)", 0.02, minval=0.0, maxval=0.5)
        self.fft_high = Float("FFT maximum frequency (cycles/pixel)", 0.45, minval=0.0, maxval=0.5)
        self.fft_plot_name = ImageName("Name the FFT orientation distribution image", "FiberFFTDistribution")
        self.fft_report_name = ImageName("Name the FFT analysis report image", "FiberFFTReport")
        self.show_advanced = Binary("Show advanced settings?", False)
        self.uncovered_image = ImageName("Name the uncovered skeleton image", "FiberUncoveredSkeleton")
        self.boundary_direction_weight = Float("Direction guidance for object boundaries", 1.0, minval=0.0, maxval=5.0)
        descriptions = {'input_image': 'Choose a 2D grayscale image with bright fibers on a dark background. Intensity '
                'measurements use this input image, not the Sato enhancement.',
 'analysis_mode': 'Entire image analyzes all valid image pixels. Within objects analyzes each '
                  'parent separately; fibers crossing a parent boundary are clipped there.',
 'input_objects': 'Select regions already identified in the pipeline, such as tissue areas. Fiber '
                  'extraction and neighborhood measurements are restricted to each parent region.',
 'output_objects': 'Name the fiber objects used by downstream modules and exported measurements. '
                   'Each assigned region pixel belongs to one object.',
 'enhanced_image': 'Name the Sato response image used to inspect enhanced structures. It does not '
                   'replace the original input for intensity measurements.',
 'skeleton_image': 'Name the binary image of detected paths before length filtering. The colored '
                   'centerlines in the overview show accepted fibers only.',
 'min_width': 'Minimum apparent fiber thickness in pixels. Decreasing it includes finer structures '
              'and possibly noise; increasing it can miss thin fibers. This is an enhancement '
              'scale, not a strict measured-width cutoff.',
 'max_width': 'Maximum apparent fiber thickness in pixels. Increasing it favors thicker '
              'structures, including bundles; decreasing it can miss broad fibers. It does not '
              'control fiber length or impose a hard width cutoff.',
 'num_scales': 'Number of thickness scales sampled between the width limits. Increasing it refines '
               'scale sampling but takes longer; it does not guarantee more fibers. With 1, only '
               'the minimum scale is used.',
 'threshold_method': 'Otsu and Li estimate a threshold from Sato enhancement. Manual and '
                     'Percentile set it explicitly. CellProfiler binary image uses an external '
                     'mask; Sato is then shown for inspection only.',
 'threshold_value': 'Threshold from 0 to 1 in normalized enhancement. Decreasing it includes faint '
                    'signal and more background; increasing it can remove weak segments and '
                    'fragment fibers.',
 'threshold_percentile': 'Percentile of enhancement used as the threshold. Increasing it generally '
                         'retains fewer pixels; decreasing it includes weaker signal. This is not '
                         'the percentage of fibers to retain.',
 'minimum_length': 'Minimum traced path length in pixels. Increasing it removes short fibers and '
                   'fragments; decreasing it preserves short branches but may retain noise. It '
                   'does not join broken fibers.',
 'minimum_area': 'Minimum connected mask-component area in square pixels. Increasing it removes '
                 'small isolated components; decreasing it retains smaller structures. It does not '
                 'filter individual fibers within a connected network.',
 'closing_radius': 'Mask closing radius in pixels. Increasing it can fill small gaps but also '
                   'connect nearby fibers; decreasing it preserves separations. Zero disables '
                   'closing.',
 'crossing_max_angle': 'Maximum deviation from straight continuation, in degrees. Increasing it '
                       'permits sharper turns but can switch fibers at crossings; decreasing it '
                       'favors straighter continuity but can truncate curved fibers.',
 'tangent_points': 'Minimum number of points used to estimate branch direction; the window also '
                   'grows with local width. Increasing it smooths irregularities but can miss '
                   'tight curves; decreasing it makes direction more sensitive to noise.',
 'threshold_correction': 'Multiplies the Otsu or Li threshold. Below 1 retains more faint signal '
                         'and background; above 1 is more restrictive and may interrupt fibers. '
                         'Test small changes while inspecting the original image.',
 'threshold_lower': 'Lower bound of the automatic threshold, from 0 to 1. Increasing it matters '
                    'only when the calculated threshold falls below this bound and can exclude '
                    'faint signal. Decreasing it permits lower thresholds.',
 'threshold_upper': 'Upper bound of the automatic threshold, from 0 to 1. Decreasing it matters '
                    'only when the calculated threshold exceeds this bound and can include faint '
                    'signal and background. Increasing it permits stricter thresholds.',
 'binary_image': 'Select a binary mask with the same dimensions as the input, for example from '
                 'Threshold. Closing and minimum component area still apply. Sato does not '
                 'determine this mask.',
 'neighbor_radius': 'Centroid distance used for neighborhood measurements. Increasing it includes '
                    'more neighbors and fewer isolated fibers; decreasing it does the opposite. It '
                    'does not change segmentation or reconnect fibers.',
 'use_fft': 'Measure orientation from the Fourier power distribution of the original image. This '
            'analysis is independent of fiber objects and does not change fiber detection.',
 'fft_low': 'Minimum frequency in cycles/pixel. Increasing it excludes broader spatial variations; '
            'decreasing it includes larger-scale structures and possible background gradients. It '
            'affects FFT only.',
 'fft_high': 'Maximum frequency in cycles/pixel. Increasing it includes finer details and possibly '
             'noise; decreasing it excludes them. It affects FFT only.',
 'fft_plot_name': 'Name the whole-image orientation-distribution plot. Save it with SaveImages. '
                  'This name is retained in pipeline settings for compatibility.',
 'fft_report_name': 'Name the whole-image report containing the input, spectrum and orientation '
                    'plots. Save it with SaveImages; even in Within objects mode, this report '
                    'represents the whole image.',
 'show_advanced': 'Show additional controls. Hiding this section does not reset defaults: '
                  'configured values remain active.',
 'uncovered_image': 'Name the image marking reference-skeleton pixels farther than 2 pixels from '
                    'any accepted path. It highlights areas to inspect; it is not ground truth or '
                    'an automatic error map.',
 'boundary_direction_weight': 'Direction weight when dividing pixels between objects. Increasing '
                              'it favors growth across centerlines; decreasing it gives distance '
                              'and intensity more influence. Zero disables direction weighting. It '
                              'does not change paths or create exclusion gaps.'}
        for name, description in descriptions.items():
            getattr(self, name).doc = description

    def settings(self):
        return [
            self.input_image, self.analysis_mode, self.input_objects, self.output_objects,
            self.enhanced_image, self.skeleton_image, self.min_width, self.max_width,
            self.num_scales, self.threshold_method, self.threshold_value,
            self.threshold_percentile, self.minimum_length, self.minimum_area,
            self.closing_radius,
            self.crossing_max_angle, self.tangent_points, self.threshold_correction,
            self.threshold_lower, self.threshold_upper,
            self.binary_image, self.neighbor_radius,
            self.use_fft, self.fft_low, self.fft_high, self.fft_plot_name, self.fft_report_name, self.show_advanced, self.uncovered_image, self.boundary_direction_weight,
        ]

    def visible_settings(self):
        result = [self.input_image, self.analysis_mode]
        if self.analysis_mode.value == self.MODE_OBJECTS:
            result.append(self.input_objects)
        result += [self.output_objects, self.min_width, self.max_width, self.threshold_method]
        if self.threshold_method.value == "Manual":
            result.append(self.threshold_value)
        elif self.threshold_method.value == "Percentile":
            result.append(self.threshold_percentile)
        elif self.threshold_method.value == "CellProfiler binary image":
            result.append(self.binary_image)
        if self.threshold_method.value in ("Otsu", "Li"):
            result.append(self.threshold_correction)
        result += [self.minimum_length, self.use_fft]
        if self.use_fft.value:
            result.append(self.fft_report_name)
        result.append(self.show_advanced)
        if self.show_advanced.value:
            result += [self.enhanced_image, self.skeleton_image, self.uncovered_image, self.num_scales, self.boundary_direction_weight]
            if self.threshold_method.value in ("Otsu", "Li"):
                result += [self.threshold_lower, self.threshold_upper]
            if self.use_fft.value:
                result += [self.fft_low, self.fft_high]
            result += [self.minimum_area, self.closing_radius]
            result += [self.crossing_max_angle, self.tangent_points, self.neighbor_radius]
        return result

    def run(self, workspace):
        input_image = workspace.image_set.get_image(self.input_image.value, must_be_grayscale=True)
        labels = None
        if self.analysis_mode.value == self.MODE_OBJECTS:
            labels = workspace.object_set.get_objects(self.input_objects.value).segmented
        settings = FiberSettings(
            tracing_method="fire2d",
            min_width=self.min_width.value,
            max_width=self.max_width.value,
            num_scales=self.num_scales.value,
            threshold_method=self.threshold_method.value.lower() if self.threshold_method.value != "CellProfiler binary image" else "otsu",
            threshold_value=self.threshold_value.value,
            threshold_percentile=self.threshold_percentile.value,
            threshold_correction=self.threshold_correction.value if self.threshold_method.value in ("Otsu", "Li") else 1.0,
            threshold_lower=self.threshold_lower.value if self.threshold_method.value in ("Otsu", "Li") else 0.0,
            threshold_upper=self.threshold_upper.value if self.threshold_method.value in ("Otsu", "Li") else 1.0,
            minimum_length=self.minimum_length.value,
            boundary_direction_weight=self.boundary_direction_weight.value,
            minimum_area=self.minimum_area.value,
            closing_radius=self.closing_radius.value,
            crossing_max_angle=self.crossing_max_angle.value,
            tangent_points=self.tangent_points.value,
        )
        binary = None
        if self.threshold_method.value == "CellProfiler binary image":
            thresholded = workspace.image_set.get_image(self.binary_image.value, must_be_binary=True)
            binary = thresholded.pixel_data & thresholded.mask
        result = analyze_fibers(input_image.pixel_data, labels, input_image.mask, settings, segmentation_mask=binary)
        measure_neighbors(result.fibers, self.neighbor_radius.value)
        objects = Objects()
        objects.segmented = result.label_image
        objects.parent_image = input_image
        workspace.object_set.add_objects(objects, self.output_objects.value)
        workspace.image_set.add(self.enhanced_image.value, Image(result.enhanced, parent_image=input_image))
        workspace.image_set.add(self.skeleton_image.value, Image(result.skeleton, parent_image=input_image))
        workspace.image_set.add(self.uncovered_image.value,
                                Image(result.uncovered_skeleton, parent_image=input_image))
        self._add_measurements(workspace, result)
        if self.use_fft.value:
            input_data = input_image.pixel_data
            mask = input_image.mask
            kwargs = {"low_frequency": self.fft_low.value, "high_frequency": self.fft_high.value}
            angle, alignment, distribution = fft_orientation(input_data, mask, return_distribution=True, **kwargs)
            fft_plot = distribution_image(distribution, angle, alignment)
            workspace.image_set.add(self.fft_plot_name.value, Image(fft_plot))
            fft_report = report_image(input_data, mask, distribution, angle, alignment, self.fft_low.value, self.fft_high.value)
            workspace.image_set.add(self.fft_report_name.value, Image(fft_report))
            workspace.measurements.add_image_measurement("FiberFFT_Orientation", angle)
            workspace.measurements.add_image_measurement("FiberFFT_Alignment", alignment)
            if labels is not None:
                count = workspace.object_set.get_objects(self.input_objects.value).count
                values = [fft_orientation(input_data, mask & (labels == index), **kwargs)
                          for index in range(1, count + 1)]
                for i, feature in enumerate(("Orientation", "Alignment")):
                    workspace.measurements.add_measurement(self.input_objects.value, "FiberFFT_" + feature,
                                                          np.asarray([v[i] for v in values], dtype=float))
        if self.show_window:
            workspace.display_data.input_image = input_image.pixel_data
            workspace.display_data.enhanced = result.enhanced
            workspace.display_data.skeleton = result.skeleton
            workspace.display_data.fiber_labels = result.label_image
            workspace.display_data.fiber_overlay = label_overlay(
                result.label_image, image=input_image.pixel_data, alpha=0.45
            )
            paths = centerline_labels(result.fibers, input_image.pixel_data.shape)
            colors = label_overlay(paths)
            gray = np.clip(np.asarray(input_image.pixel_data, dtype=float), 0, 1)
            overlay = np.repeat(gray[..., None], 3, axis=2)
            overlay[paths > 0] = colors[paths > 0]
            workspace.display_data.centerline_overlay = overlay
            workspace.display_data.fiber_count = len(result.fibers)
            workspace.display_data.coverage_percent = result.coverage_percent
            workspace.display_data.boundary_direction_weight = self.boundary_direction_weight.value
            workspace.display_data.resolved_shared_tails = result.resolved_shared_tails
            workspace.display_data.resolved_shared_segments = result.resolved_shared_segments
            display_base = np.asarray(input_image.pixel_data, dtype=float)
            display_base = np.clip(display_base, 0, 1)
            qc = np.repeat(display_base[..., None], 3, axis=2)
            qc[result.label_image > 0] = [0.0, 1.0, 0.0]
            qc[result.rejected_short] = [1.0, 0.0, 1.0]
            workspace.display_data.qc_overlay = qc
            workspace.display_data.fft_values = (distribution, angle, alignment) if self.use_fft.value else None
            workspace.display_data.dimensions = input_image.dimensions
            workspace.display_data.summary_rows = self._make_summary_rows(result.fibers)
            workspace.display_data.setting_rows = self._make_setting_rows()

    @staticmethod
    def _linear_statistics(values):
        values = np.asarray(values, dtype=float)
        values = values[np.isfinite(values)]
        if values.size == 0:
            return ["—"] * 5
        return [
            f"{np.mean(values):.3f}",
            f"{np.median(values):.3f}",
            f"{np.std(values):.3f}",
            f"{np.min(values):.3f}",
            f"{np.max(values):.3f}",
        ]

    @staticmethod
    def _axial_orientation_statistics(values):
        """Return circular mean, circular SD, and axial alignment for 0–180° data."""

        values = np.asarray(values, dtype=float)
        values = values[np.isfinite(values)]
        if values.size == 0:
            return float("nan"), float("nan"), float("nan")
        doubled = np.deg2rad(2.0 * values)
        cosine = float(np.mean(np.cos(doubled)))
        sine = float(np.mean(np.sin(doubled)))
        resultant = float(np.hypot(cosine, sine))
        mean_angle = float((np.degrees(np.arctan2(sine, cosine)) / 2.0) % 180.0) if resultant > 1e-12 else float("nan")
        resultant_for_log = np.clip(resultant, np.finfo(float).tiny, 1.0)
        circular_sd = float(np.degrees(np.sqrt(-2.0 * np.log(resultant_for_log))) / 2.0)
        return mean_angle, circular_sd, resultant

    def _make_summary_rows(self, fibers):
        if not fibers:
            return [["Fiber count", "0", "objects"]]

        rows = [["Fiber count", str(len(fibers)), "—", "—", "—", "—", "objects"]]
        specifications = [
            ("Centerline length", "centerline_length", "px"),
            ("Object major axis", "length", "pixels"),
            ("End-to-end length", "end_to_end_length", "pixels"),
            ("Mean fiber width", "width_mean", "pixels"),
            ("Median fiber width", "width_median", "pixels"),
            ("Maximum fiber width", "width_maximum", "pixels"),
            ("Straightness", "straightness", "ratio"),
            ("Tortuosity", "tortuosity", "ratio"),
            ("Mean curvature", "curvature_mean", "pixel⁻¹"),
            ("Maximum curvature", "curvature_maximum", "pixel⁻¹"),
            ("Mean intensity", "intensity_mean", "a.u."),
            ("Maximum intensity", "intensity_maximum", "a.u."),
            ("Integrated intensity", "intensity_integrated", "a.u."),
        ]
        for label, attribute, unit in specifications:
            values = [getattr(fiber, attribute) for fiber in fibers]
            rows.append([label] + self._linear_statistics(values) + [unit])

        orientations = [fiber.orientation for fiber in fibers]
        angle_mean, angle_sd, alignment = self._axial_orientation_statistics(orientations)
        finite_orientations = np.asarray(orientations, dtype=float)
        finite_orientations = finite_orientations[np.isfinite(finite_orientations)]
        angle_median = np.median(finite_orientations) if finite_orientations.size else np.nan
        angle_min = np.min(finite_orientations) if finite_orientations.size else np.nan
        angle_max = np.max(finite_orientations) if finite_orientations.size else np.nan
        rows.append([
            "Fiber orientation (axial)",
            f"{angle_mean:.3f}", f"{angle_median:.3f}", f"{angle_sd:.3f}",
            f"{angle_min:.3f}", f"{angle_max:.3f}", "degrees",
        ])
        rows.append([
            "Axial alignment coefficient", f"{alignment:.3f}", "—", "—", "—", "—", "0–1",
        ])
        total_length = float(np.nansum([fiber.length for fiber in fibers]))
        rows.append(["Total fiber length", f"{total_length:.3f}", "—", "—", "—", "—", "pixels"])
        visible = {"Fiber count", "Centerline length", "Mean fiber width", "Fiber orientation (axial)", "Axial alignment coefficient"}
        return [[row[0], row[1], row[6]] for row in rows if row[0] in visible]

    def _make_setting_rows(self):
        return [
            ["Analysis region", self.analysis_mode.value],
            ["Expected width", f"{self.min_width.value:g}–{self.max_width.value:g} px"],
            ["Threshold", self.threshold_method.value],
            ["Minimum length", f"{self.minimum_length.value:g} px"],
            ["Method", "Sato + FIRE-inspired"],
            ["Maximum deflection", f"{self.crossing_max_angle.value:g}°"],
        ]

    def _summary_measurements(self, fibers):
        def values(attribute):
            result = np.asarray([getattr(fiber, attribute) for fiber in fibers], dtype=float)
            return result[np.isfinite(result)]

        def mean(attribute):
            data = values(attribute)
            return float(np.mean(data)) if data.size else float("nan")

        lengths = values("length")
        widths = values("width_mean")
        orientations = values("orientation")
        angle_mean, angle_sd, alignment = self._axial_orientation_statistics(orientations)
        return {
            "TotalLength": float(np.sum(lengths)) if lengths.size else 0.0,
            "MeanLength": float(np.mean(lengths)) if lengths.size else float("nan"),
            "MedianLength": float(np.median(lengths)) if lengths.size else float("nan"),
            "StdLength": float(np.std(lengths)) if lengths.size else float("nan"),
            "MeanWidth": float(np.mean(widths)) if widths.size else float("nan"),
            "MedianWidth": float(np.median(widths)) if widths.size else float("nan"),
            "MaximumWidth": float(np.max(values("width_maximum"))) if fibers else float("nan"),
            "MeanOrientation": angle_mean,
            "OrientationSD": angle_sd,
            "AlignmentCoefficient": alignment,
            "MeanStraightness": mean("straightness"),
            "MeanTortuosity": mean("tortuosity"),
            "MeanCurvature": mean("curvature_mean"),
            "MeanIntensity": mean("intensity_mean"),
            "IsolatedFiberCount": float(sum(fiber.isolated for fiber in fibers)),
            "MeanNearestCentroidDistance": mean("nearest_distance"),
            "MeanNeighborAlignment": mean("neighbor_alignment"),
        }

    def display(self, workspace, figure):
        """Two-row overview with native text and FFT axes (no raster scaling)."""
        data = workspace.display_data
        figure.set_subplots((3, 2), dimensions=data.dimensions)
        original_axes = figure.subplot_imshow_grayscale(
            0, 0, data.input_image, title="Input fiber image")
        figure.subplot_imshow_grayscale(
            1, 0, data.enhanced, title="Sato enhancement", sharexy=original_axes)
        figure.subplot_imshow_color(
            2, 0, data.fiber_overlay, title=f"Fiber regions (n={data.fiber_count})",
            sharexy=original_axes)
        figure.subplot_imshow_color(
            0, 1, data.centerline_overlay, title="Centerlines", sharexy=original_axes)
        ax = figure.subplot(1, 1)
        ax.set_axis_off()
        ax.set_title("Fiber measurements", fontsize=11)
        short_labels = {"Centerline length": "Mean path length",
                        "Mean fiber width": "Mean width",
                        "Fiber orientation (axial)": "Axial orientation",
                        "Axial alignment coefficient": "Alignment"}
        table_rows = [[short_labels.get(label, label), value,
                       {"pixels": "px", "degrees": "deg"}.get(unit, unit)]
                      for label, value, unit in data.summary_rows]
        table = ax.table(cellText=table_rows,
                         colLabels=["Measurement", "Mean / value", "Unit"],
                         colWidths=[.50, .30, .20], cellLoc="left",
                         bbox=[0., .26, 1., .64])
        table.auto_set_font_size(False)
        table.set_fontsize(9)
        for (row, col), cell in table.get_celld().items():
            cell.set_edgecolor("0.65")
            cell.set_linewidth(.6)
            cell.PAD = .06
            if row == 0:
                cell.set_facecolor("0.92")
                cell.get_text().set_weight("bold")
            else:
                cell.set_facecolor("white" if row % 2 else "0.97")
                if col == 1:
                    cell.get_text().set_ha("right")
        coverage = getattr(data, "coverage_percent", float("nan"))
        ax.text(.01, .12, f"Geometric coverage (2 px): {coverage:.1f}%\n"
                f"Boundary direction guidance: {getattr(data, 'boundary_direction_weight', 1.):g}\n"
                f"Shared traces resolved (tails / segments): "
                f"{getattr(data, 'resolved_shared_tails', 0)} / {getattr(data, 'resolved_shared_segments', 0)}",
                transform=ax.transAxes, fontsize=8, color="0.35", va="top")
        ax = figure.subplot(2, 1)
        if data.fft_values is not None:
            draw_distribution(ax, *data.fft_values, fontsize=9)
        else:
            ax.set_axis_off()
            ax.set_title("Analysis settings", fontsize=11)
            for index, (key, value) in enumerate(data.setting_rows):
                ax.text(0, .9 - .12 * index, f"{key}: {value}", transform=ax.transAxes, fontsize=9)


    def _add_measurements(self, workspace, result):
        measurements = workspace.measurements
        name = self.output_objects.value
        attrs = {
            "Length": "length", "CenterlineLength": "centerline_length", "EndToEndLength": "end_to_end_length", "Straightness": "straightness",
            "Tortuosity": "tortuosity", "Orientation": "orientation",
            "WidthMean": "width_mean", "WidthMedian": "width_median",
            "WidthMaximum": "width_maximum", "CurvatureMean": "curvature_mean",
            "CurvatureMaximum": "curvature_maximum", "IntensityMean": "intensity_mean",
            "IntensityMaximum": "intensity_maximum", "IntensityIntegrated": "intensity_integrated",
            "CenterY": "center_y", "CenterX": "center_x", "StartY": "start_y",
            "StartX": "start_x", "EndY": "end_y", "EndX": "end_x",
            "NearestCentroidDistance": "nearest_distance", "NeighborAlignment": "neighbor_alignment", "Isolated": "isolated",
        }
        for feature, attr in attrs.items():
            measurements.add_measurement(name, f"Fiber_{feature}", np.asarray([getattr(f, attr) for f in result.fibers]))
        measurements.add_measurement(name, "Number_Object_Number", np.arange(1, len(result.fibers) + 1))
        measurements.add_measurement(name, "Location_Center_X", np.asarray([f.center_x for f in result.fibers]))
        measurements.add_measurement(name, "Location_Center_Y", np.asarray([f.center_y for f in result.fibers]))
        measurements.add_image_measurement(f"Count_{name}", len(result.fibers))
        qc_values = (result.coverage_percent, result.reference_skeleton_pixels,
                     int(result.uncovered_skeleton.sum()), 2., self.boundary_direction_weight.value, result.resolved_shared_tails, result.resolved_shared_segments)
        for feature, value in zip(self.QC_FEATURES, qc_values):
            measurements.add_image_measurement("FiberQC_" + feature, value)
        image_summary = self._summary_measurements(result.fibers)
        for feature, value in image_summary.items():
            measurements.add_image_measurement(f"FiberSummary_{feature}", value)
        if self.analysis_mode.value == self.MODE_OBJECTS:
            parent_name = self.input_objects.value
            parents = np.asarray([f.parent_object for f in result.fibers], dtype=int)
            measurements.add_measurement(name, f"Parent_{parent_name}", parents)
            parent_objects = workspace.object_set.get_objects(parent_name)
            child_counts = np.bincount(parents, minlength=parent_objects.count + 1)[1:]
            measurements.add_measurement(parent_name, f"Children_{name}_Count", child_counts)
            parent_summaries = []
            for parent_number in range(1, parent_objects.count + 1):
                parent_fibers = [fiber for fiber in result.fibers if fiber.parent_object == parent_number]
                parent_summaries.append(self._summary_measurements(parent_fibers))
            for feature in self.SUMMARY_FEATURES:
                measurements.add_measurement(
                    parent_name,
                    f"FiberSummary_{feature}",
                    np.asarray([summary[feature] for summary in parent_summaries], dtype=float),
                )

    def get_measurement_columns(self, pipeline):
        name = self.output_objects.value
        columns = [(name, f"Fiber_{f}", COLTYPE_FLOAT) for f in self.FEATURES]
        columns += [
            (name, "Number_Object_Number", COLTYPE_INTEGER),
            (name, "Location_Center_X", COLTYPE_FLOAT),
            (name, "Location_Center_Y", COLTYPE_FLOAT),
            ("Image", f"Count_{name}", COLTYPE_INTEGER),
        ]
        columns += [("Image", "FiberQC_" + feature, COLTYPE_FLOAT) for feature in self.QC_FEATURES]
        columns += [("Image", f"FiberSummary_{feature}", COLTYPE_FLOAT) for feature in self.SUMMARY_FEATURES]
        if self.use_fft.value:
            columns += [("Image", "FiberFFT_" + feature, COLTYPE_FLOAT) for feature in ("Orientation", "Alignment")]
        if self.analysis_mode.value == self.MODE_OBJECTS:
            columns += [
                (name, f"Parent_{self.input_objects.value}", COLTYPE_INTEGER),
                (self.input_objects.value, f"Children_{name}_Count", COLTYPE_INTEGER),
            ]
            if self.use_fft.value:
                columns += [(self.input_objects.value, "FiberFFT_" + feature, COLTYPE_FLOAT) for feature in ("Orientation", "Alignment")]
            columns += [
                (self.input_objects.value, f"FiberSummary_{feature}", COLTYPE_FLOAT)
                for feature in self.SUMMARY_FEATURES
            ]
        return columns

    def upgrade_settings(self, setting_values, variable_revision_number, module_name):
        if variable_revision_number == 1:
            setting_values = list(setting_values) + ["1.0", "0.0", "1.0"]
            variable_revision_number = 2
        if variable_revision_number == 2:
            setting_values = list(setting_values) + ["None", "30.0"]
            variable_revision_number = 3
        if variable_revision_number == 3:
            setting_values = list(setting_values) + ["No", "0.02", "0.45"]
            variable_revision_number = 4
        if variable_revision_number == 4:
            setting_values = list(setting_values) + ["FiberFFTDistribution"]
            variable_revision_number = 5
        if variable_revision_number == 5:
            setting_values = list(setting_values) + ["FiberFFTReport"]
            variable_revision_number = 6
        if variable_revision_number == 6:
            setting_values = list(setting_values) + ["No"]
            variable_revision_number = 7
        if variable_revision_number == 7:
            setting_values = list(setting_values[:-1])
            variable_revision_number = 8
        if variable_revision_number == 8:
            setting_values = list(setting_values) + ["Yes"]
            variable_revision_number = 9
        if variable_revision_number == 9:
            setting_values = list(setting_values) + ["Skeleton (current)"]
            variable_revision_number = 10
        if variable_revision_number == 10:
            variable_revision_number = 11
        if variable_revision_number == 11:
            setting_values = list(setting_values) + ["FiberUncoveredSkeleton"]
            variable_revision_number = 12
        if variable_revision_number == 12:
            setting_values = list(setting_values) + ["1.0"]
            variable_revision_number = 13
        if variable_revision_number == 13:
            # Retire skeleton-only controls and the algorithm selector. All
            # remaining values retain their positions relative to each other.
            setting_values = [value for i, value in enumerate(setting_values)
                              if i not in (15, 16, 30)]
            variable_revision_number = 14
        return setting_values, variable_revision_number

    def get_categories(self, pipeline, object_name):
        if object_name == self.output_objects.value:
            result = ["Fiber", "Location", "Number"]
            if self.analysis_mode.value == self.MODE_OBJECTS:
                result.append("Parent")
            return result
        if self.analysis_mode.value == self.MODE_OBJECTS and object_name == self.input_objects.value:
            return ["Children", "FiberSummary"] + (["FiberFFT"] if self.use_fft.value else [])
        if object_name == "Image":
            return ["Count", "FiberSummary", "FiberQC"] + (["FiberFFT"] if self.use_fft.value else [])
        return []

    def get_measurements(self, pipeline, object_name, category):
        if object_name == "Image" and category == "FiberQC":
            return list(self.QC_FEATURES)
        if self.use_fft.value and category == "FiberFFT" and (object_name == "Image" or (self.analysis_mode.value == self.MODE_OBJECTS and object_name == self.input_objects.value)):
            return ["Orientation", "Alignment"]
        if object_name == self.output_objects.value and category == "Fiber":
            return list(self.FEATURES)
        if object_name == self.output_objects.value and category == "Number":
            return ["Object_Number"]
        if object_name == self.output_objects.value and category == "Location":
            return ["Center_X", "Center_Y"]
        if object_name == self.output_objects.value and category == "Parent":
            return [self.input_objects.value]
        if object_name == self.input_objects.value and category == "Children":
            return [f"{self.output_objects.value}_Count"]
        if category == "FiberSummary" and (
            object_name == "Image"
            or (self.analysis_mode.value == self.MODE_OBJECTS and object_name == self.input_objects.value)
        ):
            return list(self.SUMMARY_FEATURES)
        if object_name == "Image" and category == "Count":
            return [self.output_objects.value]
        return []

