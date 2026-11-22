"""CellProfiler adapter for the opt-in volume workflow."""

import numpy as np
from cellprofiler_core.constants.measurement import COLTYPE_FLOAT, COLTYPE_INTEGER
from cellprofiler_core.image import Image
from cellprofiler_core.object import Objects
from cellprofiler_core.setting import ValidationError
from cellprofiler_core.setting.choice import Choice
from cellprofiler_core.setting.text import Float
from .cellprofiler_module import IdentifyFibers as PlanarModule
from .volume import FEATURES, Settings3D, analyze_volume


class IdentifyFibers(PlanarModule):
    """Identify bright fibers in a planar image or a calibrated grayscale Z-stack.

    The 2D workflow is described in the English and Portuguese manuals. For
    the private 3D test candidate, select 3D volume here and Process as 3D
    in NamesAndTypes. Enter relative voxel sizes there and the physical size
    of one Y pixel here. Review individual Z planes, especially at junctions.

    The 3D workflow uses multiscale Sato enhancement, a centerline graph and
    FIRE-inspired directional continuation. A connected bundle may contain
    several unresolved fibers. Width is an apparent mask-based diameter.

    FIRE reference: Stein AM et al. (2008). An algorithm for extracting the
    network geometry of three-dimensional collagen gels. Journal of Microscopy
    232, 463-475. https://doi.org/10.1111/j.1365-2818.2008.02141.x
    """

    __doc__ = __doc__ + "\n\n2D workflow\n-----------\n\n" + PlanarModule.__doc__
    variable_revision_number = 15
    VOLUME = "3D volume (test candidate)"
    VOLUME_DEFAULTS = ("2D image", "0.0", "0.5", "2.0", "3.0", "Sato + Otsu")
    SPACING = ("VoxelSizeZ_um", "VoxelSizeY_um", "VoxelSizeX_um")

    def create_settings(self):
        super().create_settings()
        self.workflow = Choice(
            "Image dimensionality",
            ["2D image", self.VOLUME],
            doc="Select the geometry of the acquisition. The volume mode follows each fiber through depth.",
        )
        self.y_um = Float(
            "Y pixel size (µm)",
            0.0,
            minval=0.0,
            doc="Enter the physical height of one image pixel from the acquisition. Relative Z:Y:X spacing comes from NamesAndTypes. Zero means unconfigured.",
        )
        self.vol_min = Float(
            "Minimum expected width (µm)",
            0.5,
            minval=0.001,
            doc="Lower values emphasize fine fibers and may also enhance noise. Match the scale to the apparent thickness in the acquired volume.",
        )
        self.vol_max = Float(
            "Maximum expected width (µm)",
            2.0,
            minval=0.001,
            doc="Higher values include broader bundles and increase the distance used to estimate direction at junctions. Closely spaced strands may merge.",
        )
        self.vol_length = Float(
            "Minimum traced length (µm)",
            3.0,
            minval=0.0,
            doc="Increase to remove short fragments. Decrease to keep short biological branches. This cutoff applies after tracing.",
        )
        self.vol_seg = Choice(
            "3D segmentation source",
            ["Sato + Otsu", "CellProfiler binary volume"],
            doc="Use multiscale tubeness to identify bright structures, or supply a registered binary volume. Object intensities always come from the grayscale input.",
        )
        self.input_image.doc = "Select a bright-fiber grayscale image or Z-stack. Object intensity measurements use this original signal."

    def volumetric(self):
        return True

    def settings(self):
        return super().settings() + [
            self.workflow,
            self.y_um,
            self.vol_min,
            self.vol_max,
            self.vol_length,
            self.vol_seg,
        ]

    def upgrade_settings(self, setting_values, variable_revision_number, module_name):
        if variable_revision_number <= 14:
            values, revision = super().upgrade_settings(
                setting_values, variable_revision_number, module_name
            )
            if revision == 14:
                return list(values) + list(self.VOLUME_DEFAULTS), 15
        return setting_values, variable_revision_number

    def visible_settings(self):
        if self.workflow.value != self.VOLUME:
            return [self.workflow] + super().visible_settings()
        fields = [self.workflow, self.input_image, self.analysis_mode]
        if self.analysis_mode.value == self.MODE_OBJECTS:
            fields.append(self.input_objects)
        fields += [self.output_objects, self.y_um, self.vol_seg]
        if self.vol_seg.value == "Sato + Otsu":
            fields.extend([self.vol_min, self.vol_max, self.threshold_correction])
        else:
            fields.extend([self.binary_image, self.vol_max])
        fields.extend([self.vol_length, self.show_advanced])
        if self.show_advanced.value:
            fields.extend([self.crossing_max_angle, self.skeleton_image])
            if self.vol_seg.value == "Sato + Otsu":
                fields.append(self.enhanced_image)
        return fields

    def validate_module(self, pipeline):
        super().validate_module(pipeline)
        if self.workflow.value == self.VOLUME:
            if self.y_um.value <= 0:
                raise ValidationError(
                    "Enter the acquired Y pixel size in micrometres.", self.y_um
                )
            if (
                self.vol_seg.value == "Sato + Otsu"
                and self.vol_min.value > self.vol_max.value
            ):
                raise ValidationError(
                    "Maximum width must be at least minimum width.", self.vol_max
                )

    def _volume_schema(self):
        obj = self.output_objects.value
        schema = [(obj, "Fiber3D", key, COLTYPE_FLOAT) for key in FEATURES]
        schema += [(obj, "Location", "Center_" + axis, COLTYPE_FLOAT) for axis in "XYZ"]
        schema += [
            (obj, "Number", "Object_Number", COLTYPE_INTEGER),
            ("Image", "Count", obj, COLTYPE_INTEGER),
        ]
        schema += [("Image", "Fiber3D", key, COLTYPE_FLOAT) for key in self.SPACING]
        if self.analysis_mode.value == self.MODE_OBJECTS:
            parent = self.input_objects.value
            schema += [
                (obj, "Parent", parent, COLTYPE_INTEGER),
                (parent, "Children", obj + "_Count", COLTYPE_INTEGER),
            ]
        return schema

    def get_measurement_columns(self, pipeline):
        if self.workflow.value != self.VOLUME:
            return super().get_measurement_columns(pipeline)
        return [(obj, cat + "_" + f, typ) for obj, cat, f, typ in self._volume_schema()]

    def get_categories(self, pipeline, object_name):
        if self.workflow.value != self.VOLUME:
            return super().get_categories(pipeline, object_name)
        return sorted(
            {cat for obj, cat, _, _ in self._volume_schema() if obj == object_name}
        )

    def get_measurements(self, pipeline, object_name, category):
        if self.workflow.value != self.VOLUME:
            return super().get_measurements(pipeline, object_name, category)
        return [
            f
            for obj, cat, f, _ in self._volume_schema()
            if obj == object_name and cat == category
        ]

    def run(self, workspace):
        source = workspace.image_set.get_image(
            self.input_image.value, must_be_grayscale=True
        )
        if self.workflow.value != self.VOLUME:
            if source.dimensions != 2:
                raise ValueError("Select 3D volume in IdentifyFibers for this stack.")
            return super().run(workspace)
        if source.dimensions != 3 or source.pixel_data.ndim != 3:
            raise ValueError(
                "Select Process as 3D in NamesAndTypes and load a grayscale Z-stack."
            )
        if not np.isfinite(self.y_um.value) or self.y_um.value <= 0:
            raise ValueError("Set the physical size of one Y pixel in micrometres.")
        spacing = np.asarray(source.spacing) * self.y_um.value
        binary = None
        if self.vol_seg.value != "Sato + Otsu":
            bm = workspace.image_set.get_image(
                self.binary_image.value, must_be_binary=True
            )
            if (
                bm.dimensions != 3
                or bm.pixel_data.shape != source.pixel_data.shape
                or not np.allclose(bm.spacing, source.spacing)
            ):
                raise ValueError(
                    "The binary volume must have the same shape and relative spacing as the input."
                )
            binary = bm.pixel_data & bm.mask
        parents = None
        if self.analysis_mode.value == self.MODE_OBJECTS:
            parents = workspace.object_set.get_objects(self.input_objects.value)
            if parents.parent_image is not None and (
                parents.parent_image.dimensions != 3
                or not np.allclose(parents.parent_image.spacing, source.spacing)
            ):
                raise ValueError(
                    "Parent objects and input volume must use the same voxel spacing."
                )
            if parents.overlapping():
                raise ValueError(
                    "3D parent objects must have exclusive voxel ownership."
                )
        cfg = Settings3D(
            self.vol_min.value if binary is None else self.vol_max.value,
            self.vol_max.value,
            self.vol_length.value,
            self.threshold_correction.value,
            self.crossing_max_angle.value,
        )
        result = analyze_volume(
            source.pixel_data,
            spacing,
            cfg,
            source.mask,
            parents.segmented if parents is not None else None,
            binary,
        )
        obj = Objects()
        obj.segmented = result.labels
        obj.parent_image = source
        name = self.output_objects.value
        workspace.object_set.add_objects(obj, name)
        if binary is None:
            workspace.image_set.add(
                self.enhanced_image.value,
                Image(result.enhanced, dimensions=3, parent_image=source),
            )
        workspace.image_set.add(
            self.skeleton_image.value,
            Image(result.centers, dimensions=3, parent_image=source),
        )
        m = workspace.measurements
        fibers = result.fibers
        m.add_image_measurement("Count_" + name, len(fibers))
        for key, v in zip(self.SPACING, spacing):
            m.add_image_measurement("Fiber3D_" + key, v)
        for key in FEATURES:
            m.add_measurement(
                name, "Fiber3D_" + key, np.asarray([f[key] for f in fibers], float)
            )
        for i, axis in enumerate("ZYX"):
            m.add_measurement(
                name,
                "Location_Center_" + axis,
                np.asarray([f["centroid"][i] for f in fibers], float),
            )
        m.add_measurement(
            name, "Number_Object_Number", np.arange(1, len(fibers) + 1, dtype=int)
        )
        if parents is not None:
            ids = np.array([f["parent"] for f in fibers], dtype=int)
            m.add_measurement(name, "Parent_" + self.input_objects.value, ids)
            m.add_measurement(
                self.input_objects.value,
                "Children_" + name + "_Count",
                np.bincount(ids, minlength=parents.count + 1)[1:],
            )
        if self.show_window:
            d = workspace.display_data
            d.source = source.pixel_data
            d.labels = result.labels
            d.detection = result.enhanced if binary is None else binary
            d.centers = np.where(result.centers, result.labels, 0)
            d.count = len(fibers)

    def display(self, workspace, figure):
        if self.workflow.value != self.VOLUME:
            return super().display(workspace, figure)
        d = workspace.display_data
        figure.set_subplots((2, 2), dimensions=3)
        ax = figure.subplot_imshow_grayscale(0, 0, d.source, title="Input fiber volume")
        figure.subplot_imshow_grayscale(
            1, 0, d.detection, title=self.vol_seg.value, sharexy=ax
        )
        figure.subplot_imshow_labels(
            0, 1, d.labels, title=f"Fiber objects (n={d.count})", sharex=ax, sharey=ax
        )
        figure.subplot_imshow_labels(
            1, 1, d.centers, title="Centerlines", sharex=ax, sharey=ax
        )
