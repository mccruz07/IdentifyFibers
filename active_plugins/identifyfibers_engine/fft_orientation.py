"""Experimental Fourier-based axial image orientation, not FiberFit."""
import numpy as np


def fft_orientation(image, mask=None, low_frequency=0.02, high_frequency=0.45, return_distribution=False):
    """Return dominant axial angle and power-weighted axial resultant.

    Frequencies are cycles/pixel. Image x points right and y points down.
    Regions are cropped to their bounding rectangle; mask edges can bias FFT.
    """
    image = np.asarray(image, dtype=float)
    if image.ndim != 2 or not (0 <= low_frequency < high_frequency <= 0.5):
        raise ValueError("Expected a 2D image and 0 <= low < high <= 0.5 cycles/pixel")
    valid = np.isfinite(image) if mask is None else (np.asarray(mask, bool) & np.isfinite(image))
    if valid.shape != image.shape:
        raise ValueError("FFT mask must match image dimensions")
    if np.count_nonzero(valid) < 16:
        return (np.nan, np.nan, None) if return_distribution else (np.nan, np.nan)
    ys, xs = np.where(valid)
    box = np.s_[ys.min():ys.max()+1, xs.min():xs.max()+1]
    data, region = image[box], valid[box]
    if min(data.shape) < 4:
        return (np.nan, np.nan, None) if return_distribution else (np.nan, np.nan)
    centered = np.where(region, data - np.mean(data[region]), 0.0)
    window = np.outer(np.hanning(data.shape[0]), np.hanning(data.shape[1]))
    power = np.abs(np.fft.fft2(centered * window)) ** 2
    fy = np.fft.fftfreq(data.shape[0])[:, None]
    fx = np.fft.fftfreq(data.shape[1])[None, :]
    radius = np.hypot(fy, fx)
    use = (radius >= low_frequency) & (radius <= high_frequency) & (radius > 0)
    weights = np.where(use, power, 0.0)
    total = float(weights.sum())
    if not np.isfinite(total) or total <= np.finfo(float).eps:
        return (np.nan, np.nan, None) if return_distribution else (np.nan, np.nan)
    angle = np.arctan2(fy, fx) + np.pi / 2
    c = float(np.sum(weights * np.cos(2 * angle)) / total)
    s = float(np.sum(weights * np.sin(2 * angle)) / total)
    alignment = float(np.clip(np.hypot(c, s), 0, 1))
    orientation = float((np.degrees(np.arctan2(s, c)) / 2) % 180) if alignment > 1e-12 else np.nan
    if return_distribution:
        degrees = np.broadcast_to(np.mod(np.degrees(angle), 180), power.shape)
        histogram, edges = np.histogram(degrees[use], bins=36, range=(0, 180), weights=weights[use])
        return orientation, alignment, (edges, histogram / total)
    return orientation, alignment


def draw_distribution(ax, distribution, orientation, alignment, fontsize=10):
    """Draw into native axes: labels stay readable when the window resizes."""
    if distribution is None:
        ax.text(.5, .5, "No valid FFT signal", ha="center", va="center", transform=ax.transAxes)
    else:
        edges, weights = distribution
        ax.bar(edges[:-1], weights, width=np.diff(edges), align="edge", color="#428bca")
        if np.isfinite(orientation):
            ax.axvline(orientation, color="crimson", linestyle="--", linewidth=1.2)
            ax.text(.98, .96, f"Orientation: {orientation:.1f}°", ha="right", va="top",
                    transform=ax.transAxes, fontsize=fontsize,
                    bbox=dict(facecolor="white", edgecolor="none", alpha=.8))
    label = f"{alignment:.3f}" if np.isfinite(alignment) else "undefined"
    ax.set(xlim=(0, 180), xlabel="Orientation (degrees)", ylabel="FFT power fraction",
           title=f"FFT alignment: {label}")
    ax.set_xticks(np.arange(0, 181, 30))
    ax.tick_params(labelsize=fontsize)
    ax.xaxis.label.set_size(fontsize)
    ax.yaxis.label.set_size(fontsize)
    ax.title.set_size(fontsize + 1)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


def distribution_image(distribution, orientation, alignment):
    """Render the binned axial spectrum as a CellProfiler-compatible RGB image."""
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    fig = Figure(figsize=(8, 5), dpi=150)
    canvas = FigureCanvasAgg(fig)
    draw_distribution(fig.add_subplot(111), distribution, orientation, alignment, fontsize=14)
    fig.tight_layout()
    canvas.draw()
    return np.asarray(canvas.buffer_rgba())[..., :3].astype(np.float32) / 255.0


def report_image(image, mask, distribution, orientation, alignment, low, high):
    """Create an illustrative full-image FFT report, not a FiberFit fit."""
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    data = np.asarray(image, dtype=float)
    valid = np.asarray(mask, bool) & np.isfinite(data)
    fig = Figure(figsize=(12, 12), dpi=150)
    canvas = FigureCanvasAgg(fig)
    grid = fig.add_gridspec(3, 2, height_ratios=[1, 1, .24])
    ax = fig.add_subplot(grid[0, 0]); ax.imshow(data, cmap="gray"); ax.set_title("Input image"); ax.axis("off")
    ax = fig.add_subplot(grid[0, 1])
    if valid.sum() >= 16:
        ys, xs = np.where(valid)
        sl = np.s_[ys.min():ys.max()+1, xs.min():xs.max()+1]
        crop, m = data[sl], valid[sl]
        centered = np.where(m, crop - np.mean(crop[m]), 0)
        win = np.outer(np.hanning(crop.shape[0]), np.hanning(crop.shape[1]))
        power = np.fft.fftshift(np.abs(np.fft.fft2(centered * win)) ** 2)
        fy = np.fft.fftshift(np.fft.fftfreq(crop.shape[0]))[:, None]
        fx = np.fft.fftshift(np.fft.fftfreq(crop.shape[1]))[None, :]
        freq = np.hypot(fy, fx)
        power[(freq < low) | (freq > high)] = 0
        ax.imshow(np.log1p(power), cmap="gray")
    ax.set_title("Band-pass FFT power (log)"); ax.axis("off")
    polar = fig.add_subplot(grid[1, 0], projection="polar")
    hist = fig.add_subplot(grid[1, 1])
    if distribution is not None:
        edges, weights = distribution
        centers = (edges[:-1] + edges[1:]) / 2
        theta = np.deg2rad(centers)
        polar.plot(np.r_[theta, theta + np.pi, theta[0] + 2*np.pi], np.r_[weights, weights, weights[0]], color="black")
        if np.isfinite(orientation):
            for direction in (orientation, orientation + 180):
                polar.axvline(np.deg2rad(direction), color="crimson", linestyle="--")
            hist.axvline(orientation, color="crimson", linestyle="--", label=f"Mean {orientation:.1f}°")
            hist.legend()
        hist.bar(edges[:-1], weights, width=np.diff(edges), align="edge", color="#428bca")
    else:
        hist.text(.5, .5, "No valid FFT signal", ha="center", transform=hist.transAxes)
    polar.set_theta_zero_location("E"); polar.set_theta_direction(-1); polar.set_title("Axial angular distribution")
    hist.set(xlim=(0, 180), xlabel="Orientation (degrees)", ylabel="FFT power fraction", title="Angular distribution")
    ax = fig.add_subplot(grid[2, :]); ax.axis("off")
    ax.text(.5, .7, f"Axial mean orientation: {orientation:.2f}°     Axial alignment: {alignment:.3f}", ha="center", fontsize=12)
    ax.text(.5, .25, f"Frequency band: {low:g}–{high:g} cycles/pixel", ha="center", fontsize=10)
    hist.tick_params(labelsize=12)
    hist.xaxis.label.set_size(13); hist.yaxis.label.set_size(13)
    hist.title.set_size(15); polar.title.set_size(15)
    fig.tight_layout(); canvas.draw()
    return np.asarray(canvas.buffer_rgba())[..., :3].astype(np.float32) / 255.0
