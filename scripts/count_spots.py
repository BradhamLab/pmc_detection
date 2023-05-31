import json
import os
import sys
from collections import namedtuple
import logging
import itertools
from turtle import dot

import numpy as np
import xarray as xr
import skimage
from skimage import exposure, io, morphology, measure
from skimage.filters import threshold_otsu

import bigfish.detection as bf_detection
import bigfish.stack as bf_stack
import bigfish.plot as bf_plot

BoundingBox = namedtuple("BoundingBox", ["ymin", "ymax", "xmin", "xmax"])
try:
    sys.path.append(os.path.dirname(__file__))
except NameError:
    pass
import utils


def normalize_zstack(z_stack, bits):
    """Normalize z-slices in a 3D image using contrast stretching

    Parameters
    ----------
    z_stack : numpy.ndarray
        3 dimensional confocal FISH image.
    bits : int
        Bit depth of image.

    Returns
    -------
    numpy.ndarray
        Z-corrected image with each slice minimum and maximum matched
    """
    out = np.array(
        [
            exposure.rescale_intensity(
                x, in_range=(0, 2 ** bits - 1), out_range=(z_stack.min(), z_stack.max())
            )
            for x in z_stack
        ]
    )
    return skimage.img_as_uint(exposure.rescale_intensity(out))


def read_bit_img(img_file, bits=12):
    """Read an image and return as a 16-bit image."""
    img = exposure.rescale_intensity(
        io.imread(img_file),
        in_range=(0, 2 ** (bits) - 1)
        #         out_range=(0, )
    )
    return skimage.img_as_uint(img)


def select_signal(image, p_in_focus=0.75, margin_width=10):
    """
    Generate bounding box of FISH image to select on areas where signal is present.

    Parameters
    ----------
    image : np.ndarray
        3D FISH image
    p_in_focus : float, optional
        Percent of in-focus slices to retain for 2D projection, by default 0.75.
    margin_width : int, optional
        Number of pixels to pad selection by. Default is 10.

    Returns
    -------
    namedtuple
        minimum and maximum coordinate values of the bounding box in the xy plane
    """
    image = image.astype(np.uint16)
    focus = bf_stack.compute_focus(image)
    selected = bf_stack.in_focus_selection(image, focus, p_in_focus)
    projected_2d = bf_stack.maximum_projection(selected)
    foreground = np.where(projected_2d > threshold_otsu(projected_2d))
    limits = BoundingBox(
        ymin=max(foreground[0].min() - margin_width, 0),
        ymax=min(foreground[0].max() + margin_width, image.shape[1]),
        xmin=max(foreground[1].min() - margin_width, 0),
        xmax=min(foreground[1].max() + margin_width, image.shape[2]),
    )
    return limits


def crop_to_selection(img, bbox):
    """
    Crop image to selection defined by bounding box.

    Crops a 3D image to specified x and y coordinates.
    Parameters
    ----------
    img : np.ndarray
        3Dimensional image to crop
    bbox : namedtuple
        Tuple defining minimum and maximum coordinates for x and y planes.

    Returns
    -------
    np.ndarray
        3D image cropped to the specified selection.
    """
    return img[:, bbox.ymin : bbox.ymax, bbox.xmin : bbox.xmax]


def count_spots_in_labels(spots, labels):
    """
    Count the number of RNA molecules in specified labels.

    Parameters
    ----------
    spots : np.ndarray
        Coordinates in original image where RNA molecules were detected.
    labels : np.ndarray
        Integer array of same shape as `img` denoting regions to interest to quantify.
        Each separate region should be uniquely labeled.

    Returns
    -------
    dict
        dictionary containing the number of molecules contained in each labeled region.
    """
    assert spots.shape[1] == len(labels.shape)
    n_labels = np.unique(labels) - 1  # subtract one for backgroudn
    counts = {i: 0 for i in range(1, n_labels + 1)}
    for each in spots:
        if len(each) == 3:
            cell_label = labels[each[0], each[1], each[2]]
        else:
            cell_label = labels[each[0], each[1]]
        if cell_label != 0:
            counts[cell_label] += 1
    return counts


def preprocess_image(
    img, smooth_method="gaussian", sigma=7, whitehat=True, selem=None, stretch=99.99
):
    scaled = exposure.rescale_intensity(
        img, in_range=tuple(np.percentile(img, [0, stretch]))
    )
    if smooth_method == "log":
        smooth_func = bf_stack.log_filter
        to_smooth = bf_stack.cast_img_float64(scaled)
    elif smooth_method == "gaussian":
        smooth_func = bf_stack.remove_background_gaussian
        to_smooth = bf_stack.cast_img_uint16(scaled)
    else:
        raise ValueError(f"Unsupported background filter: {smooth_method}")
    if whitehat:
        f = lambda x, s: morphology.white_tophat(smooth_func(x, s), selem)
    else:
        f = lambda x, s: smooth_func(x, s)
    smoothed = np.stack([f(img_slice, sigma) for img_slice in to_smooth])
    return bf_stack.cast_img_float64(np.stack(smoothed))


def count_spots(
    smoothed_signal,
    cell_labels,
    voxel_size_nm,
    dot_radius_nm,
    smooth_method="gaussian",
    decompose_alpha=0.5,
    decompose_beta=1,
    decompose_gamma=5,
    verbose=False,
):

    if verbose:
        spot_radius_px = bf_detection.get_object_radius_pixel(
            voxel_size_nm=voxel_size_nm, object_radius_nm=dot_radius_nm, ndim=3
        )
        logging.info("spot radius (z axis): %0.3f pixels", spot_radius_px[0])
        logging.info("spot radius (yx plan): %0.3f pixels", spot_radius_px[-1])
    spots, threshold = bf_detection.detect_spots(
        smoothed_signal,
        return_threshold=True,
        voxel_size=voxel_size_nm,
        spot_radius=dot_radius_nm,
    )
    if verbose:
        logging.info("%d spots detected...", spots.shape[0])
        logging.info("plotting threshold optimization for spot detection...")
        bf_plot.plot_elbow(
            smoothed_signal,
            voxel_size=voxel_size_nm,
            spot_radius=dot_radius_nm,
        )
    decompose_cast = {
        "gaussian": bf_stack.cast_img_uint16,
        "log": bf_stack.cast_img_float64,
    }
    try:
        (
            spots_post_decomposition,
            dense_regions,
            reference_spot,
        ) = bf_detection.decompose_dense(
            decompose_cast[smooth_method](smoothed_signal),
            spots,
            voxel_size=voxel_size_nm,
            spot_radius=dot_radius_nm,
            alpha=decompose_alpha,  # alpha impacts the number of spots per candidate region
            beta=decompose_beta,  # beta impacts the number of candidate regions to decompose
            gamma=decompose_gamma,  # gamma the filtering step to denoise the image
        )
        logging.info(
            "detected spots before decomposition: %d\n"
            "detected spots after decomposition: %d",
            spots.shape[0],
            spots_post_decomposition.shape[0],
        )
        if verbose:
            print(
                f"detected spots before decomposition: {spots.shape[0]}\n"
                f"detected spots after decomposition: {spots_post_decomposition.shape[0]}\n"
                f"shape of reference spot for decomposition: {reference_spot.shape}"
            )
            bf_plot.plot_reference_spot(reference_spot, rescale=True)
    except RuntimeError:
        logging.warning("decomposition failed, using originally identified spots")
        spots_post_decomposition = spots
    n_labels = len(np.unique(cell_labels)) - 1
    counts = {i: 0 for i in range(1, n_labels + 1)}
    expression_3d = np.zeros_like(smoothed_signal)
    # get slices to account for cropping
    for each in spots_post_decomposition:
        spot_coord = tuple(each)
        cell_label = cell_labels[spot_coord]
        if cell_label != 0:
            counts[cell_label] += 1
    for region in measure.regionprops(cell_labels):
        expression_3d[region.slice][region.image] = counts[region.label]
    return counts, expression_3d


def average_intensity(smoothed_signal, cell_labels):
    n_labels = len(np.unique(cell_labels)) - 1
    intensities = {i: 0 for i in range(1, n_labels + 1)}
    z_normed_smooth = (smoothed_signal - smoothed_signal.mean()) / smoothed_signal.std()
    expression_3d = np.zeros_like(z_normed_smooth)
    for region in measure.regionprops(cell_labels, z_normed_smooth):
        intensities[region.label] = region.mean_intensity
        expression_3d[region.slice][region.image] = region.mean_intensity
    return intensities, expression_3d


def quantify_expression(
    fish_img,
    cell_labels,
    measures=["spots", "intensity"],
    voxel_size_nm=None,
    dot_radius_nm=None,
    whitehat=True,
    whitehat_selem=None,
    smooth_method="gaussian",
    smooth_sigma=1,
    decompose_alpha=0.5,
    decompose_beta=1,
    decompose_gamma=5,
    bits=12,
    crop_image=True,
    verbose=False,
):
    """
    Count the number of molecules in an smFISH image

    Parameters
    ----------
    fish_img : np.ndarray
        Image in which to perform molecule counting
    cell_labels : np.ndarray
        Integer array of same shape as `img` denoting regions to interest to quantify.
        Each separate region should be uniquely labeled.
    measures : list
        Measures to use to quantify expression. Possible values are "spots",
        "intensity", and ["spots", "intensity"]. Default is to measure both
    voxel_size_nm : tuple(float, int), None
        Physical dimensions of each voxel in ZYX order. Required if running spot
        counting.
    dot_radius_nm : tuple(float, int), None
        Physical size of expected dots. Required if running spot
        counting.
    whitehat : bool, optional
        Whether to perform white tophat filtering prior to image de-noising, by default True
    whitehat_selem : [int, np.ndarray], optional
        Structuring element to use for white tophat filtering.
    smooth_method : str, optional
        Method to use for image de-noising. Possible values are "log" and "gaussian" for
        Laplacian of Gaussians and Gaussian background subtraction, respectively. By default "log".
    smooth_sigma : [int, np.ndarray], optional
        Sigma value to use for smoothing function, by default 1
    decompose_alpha : float, optional
        Intensity percentile used to compute the reference spot, between 0 and 1.
        By default 0.7. For more information, see:
        https://big-fish.readthedocs.io/en/stable/detection/dense.html
    decompose_beta : int, optional
        Multiplicative factor for the intensity threshold of a dense region,
        by default 1. For more information, see:
        https://big-fish.readthedocs.io/en/stable/detection/dense.html
    decompose_gamma : int, optional
        Multiplicative factor use to compute a gaussian scale, by default 5.
        For more information, see:
        https://big-fish.readthedocs.io/en/stable/detection/dense.html
    bits : int, optional
        Bit depth of original image. Used for scaling image while maintaining
        ob
    crop_image : bool, optional
        Whether to crop signal. Default is True.
    verbose : bool, optional
        Whether to verbosely print results and progress.

    Returns
    -------
    (np.ndarray, dict)
        np.ndarray: positions of all identified mRNA molecules.
        dict: dictionary containing the number of molecules contained in each labeled region.
    """
    if (voxel_size_nm is None or dot_radius_nm is None) and "spots" in measures:
        raise ValueError(
            "Require `voxel_size_nm` and `dot_radius_nm` when performing spot counting."
        )
    if crop_image:
        if verbose:
            logging.info("Cropping image to signal")
        limits = select_signal(fish_img)
        if verbose:
            logging.info(
                "Cropped image to %d x %d",
                {limits.ymax - limits.ymin},
                {limits.xmax - limits.xmin},
            )
    else:
        # create BoudndingBox that selects whole image
        limits = BoundingBox(
            ymin=0, ymax=fish_img.shape[1], xmin=0, xmax=fish_img.shape[2]
        )
    if verbose:
        logging.info("Preprocessing image.")
    cropped_img = skimage.img_as_float64(
        exposure.rescale_intensity(
            crop_to_selection(fish_img, limits),
            in_range=(0, 2 ** bits - 1),
            out_range=(0, 1),
        )
    )
    smoothed = preprocess_image(
        cropped_img, smooth_method, smooth_sigma, whitehat, whitehat_selem, 99.99
    )

    cropped_labels = crop_to_selection(cell_labels, limits)
    quant = dict()
    if "spots" in measures:
        counts, counts_3d = count_spots(
            smoothed,
            cropped_labels,
            voxel_size_nm=voxel_size_nm,
            dot_radius_nm=dot_radius_nm,
            smooth_method=smooth_method,
            decompose_alpha=decompose_alpha,
            decompose_beta=decompose_beta,
            decompose_gamma=decompose_gamma,
            verbose=verbose,
        )
        quant["spots"] = counts
        if crop_image:  # match original shape if cropped
            counts_3d = utils.pad_to_shape(counts_3d, fish_img.shape)
    if "intensity" in measures:
        intensities, intense_3d = average_intensity(smoothed, cropped_labels)
        quant["intensity"] = intensities
        if crop_image:  # match original shape if cropped
            intense_3d = utils.pad_to_shape(intense_3d, fish_img.shape)
    if len(measures) == 2:
        expression_3d = np.stack([counts_3d, intense_3d])
    return (
        quant,
        expression_3d,
    )


def get_quant_measure(method):
    """Get method of quantification"""
    if method == "spots":
        return ["spots"]
    elif method == "intensity":
        return ["intensity"]
    elif method == "both":
        return ["spots", "intensity"]
    else:
        raise ValueError(f"Unrecognized method {method}.")


if __name__ == "__main__":
    import h5py
    import pandas as pd

    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        logging.basicConfig(filename=snakemake.log[0], level=logging.INFO)
        raw_img, dimensions = utils.read_image_file(
            snakemake.input["image"], as_nm=True
        )
        if dimensions is None and snakemake.input["image"].endswith(".h5"):
            with open(snakemake.params["dimensions"], "r") as handle:
                data = json.load(handle)
                dimensions = np.array([data[c] for c in "zyx"]) * (10 ** 3)
        labels = np.array(h5py.File(snakemake.input["labels"], "r")["image"])
        logging.info("%d labels detected.", len(np.unique(labels) - 1))
        start = int(snakemake.params["z_start"])
        stop = int(snakemake.params["z_end"])
        if stop < 0:
            stop = raw_img.shape[1]
        gene_params = snakemake.params["gene_params"]

        def has_probe_info(name, gene_params):
            if name not in gene_params.keys():
                logging.warning(
                    "No entry for %s found in gene parameters. Not quantifying signal",
                    name,
                )
                return False
            return True

        channels = {
            x: i
            for i, x in enumerate(snakemake.params["channels"].split(";"))
            if has_probe_info(x, gene_params)
        }
        if len(channels) < 0:
            raise ValueError(
                f"No quantification parameters provided for channels: {snakemake.params['channels'].replace(';', ', ')}"
            )
        genes = list(channels.keys())
        fish_exprs = {}
        summarized_images = [None] * len(channels)
        embryo = snakemake.wildcards["embryo"]
        measures = get_quant_measure(snakemake.params["quant_method"])
        for i, gene in enumerate(genes):
            logging.info(f"Quantifying {gene} signal...")
            fish_data = raw_img[channels[gene], start:stop, :, :]
            quant, image = quantify_expression(
                fish_data,
                labels,
                measures=measures,
                voxel_size_nm=dimensions.tolist(),
                dot_radius_nm=gene_params[gene]["radius"],
                whitehat=True,
                smooth_method="gaussian",
                smooth_sigma=7,
                verbose=True,
                bits=12,
                crop_image=snakemake.params["crop_image"],
            )
            for each in measures:
                fish_exprs[f"{gene}_{each}"] = quant[each]
            summarized_images[i] = image
        # write summarized expression images to netcdf using Xarray to keep
        # track of dims
        out_image = np.array(summarized_images)
        xr.DataArray(
            data=out_image,
            coords={"gene": genes, "measure": measures},
            dims=["gene", "measure", "Z", "Y", "X"],
        ).to_netcdf(snakemake.output["image"])
        exprs_df = pd.DataFrame.from_dict(fish_exprs)
        exprs_df.index.name = "label"
        physical_properties = (
            pd.DataFrame(
                measure.regionprops_table(
                    labels,
                    properties={
                        "label",
                        "centroid",
                        "area",
                        "equivalent_diameter",
                    },
                )
            )
            .rename(
                columns={
                    "centroid-0": "Z",
                    "centroid-1": "Y",
                    "centroid-2": "X",
                    "equivalent_diameter": "diameter",
                }
            )
            .set_index("label")
        )
        out = exprs_df.join(physical_properties)
        out["embryo"] = embryo
        out[
            [
                "embryo",
                "area",
                "diameter",
                "Z",
                "Y",
                "X",
            ]
            + [
                f"{gene}_{measure}"
                for (gene, measure) in itertools.product(genes, measures)
            ]
        ].to_csv(snakemake.output["csv"])
