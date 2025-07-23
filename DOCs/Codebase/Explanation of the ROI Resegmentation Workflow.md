# Explanation of the ROI Resegmentation Workflow

This document explains the logic and implementation of the `resegment_ROI_workflow`. This workflow is designed to perform cell segmentation on a specific Region of Interest (ROI) from your MERSCOPE data, which is significantly faster than processing the entire tissue sample.

## The Goal

The primary objective is to re-run segmentation (e.g., with a new custom model) on a small, user-defined rectangular region of the tissue without needing to process the entire, large mosaic image.

## The Challenge

The Vizgen Postprocessing Toolkit (`vpt`) is a powerful tool, but it does not have a simple, built-in command for ROI-based processing. Our debugging process revealed several challenges that required a custom workflow to solve:

1.  **No Direct ROI Input**: There is no argument in `vpt` to simply provide ROI coordinates and have it handle the rest.
2.  **Coordinate System Mismatch**: The coordinates you select in the MERSCOPE Visualizer are in **micron space**, but the segmentation tiles (`vpt`'s internal processing units) are defined in **pixel space**. Our initial attempts to transform the ROI to pixel space were incorrect. The final solution required inverting the logic to match the tool's internal behavior.
3.  **Inconsistent Tool Requirements**: The different sub-commands of `vpt` have very specific and different requirements:
    *   `vpt run-segmentation-on-tile` needs the **full** specification file and the **original index** of the tile to process.
    *   `vpt compile-tile-segmentation` needs a specification file that **only lists the tiles that were actually processed**. If it's given the full spec file, it will fail when it can't find the results for tiles outside the ROI.
4.  **Silent Failures & Unpredictable Naming**: We discovered that running `vpt` on multiple tiles in parallel could cause silent failures where the process would finish without error but fail to write the output file. The naming of the output files was also found to be unpredictable.

## The Solution

Our final workflow is a carefully orchestrated sequence of steps that solves all of these challenges. It is composed of two custom scripts that work together.

### The Files

#### 1. `roi_coords.json`

This is the simple configuration file where you define your target region. The coordinates must be in **micron space**, as copied from the MERSCOPE Visualizer.

```json
{
  "x_min": 10691.9,
  "y_min": 6190.08,
  "x_max": 11990.2,
  "y_max": 6940.36
}
```

#### 2. `filter_spec_and_get_indices_v2.py`

This Python script is the core of the custom logic. It performs several critical tasks:

1.  **Inverse Coordinate Transformation**: It reads the `micron_to_mosaic_pixel_transform.csv` file and calculates the **inverse** of the transformation matrix. It then uses this inverse matrix to convert the **pixel-space coordinates of each tile** into **micron-space**.
2.  **Filtering and Indexing**: It compares the transformed, micron-space tile boundaries to your original micron-space ROI. This is a much more robust comparison.
    *   It creates a **new, filtered specification file** (`segmentation_specification_roi.json`) that contains *only* the tiles that overlap with your ROI. This file is essential for the `compile` step.
    *   It prints the **original indices** of these overlapping tiles to the console. This list of indices is essential for the `run-on-tile` step.

#### 3. `run_roi_segmentation_custom.sh`

This is the main script that executes the entire workflow from start to finish. Here is a breakdown of its final, working logic:

*   **Step 1-2 (Prepare)**: It calls `vpt prepare-segmentation` to generate the full `segmentation_specification.json` for the entire sample.
*   **Step 3 (Filter & Get Indices)**: It runs our custom Python script (`filter_spec_and_get_indices_v2.py`). It saves the filtered specification file and captures the list of original tile indices into a shell variable (`$ROI_TILE_INDICES`).
*   **Step 4 (Sequential Run)**: It loops through the captured original tile indices **one by one**. We discovered that running the `vpt` commands in parallel was unreliable. This sequential loop guarantees that each tile is processed completely without interference.
*   **Step 5 (Rename Files)**: This is the robust replacement for the symlink workaround. It finds **whatever** `.parquet` files were created by the previous step, regardless of their unpredictable names, and **renames** them sequentially (`cell_0.parquet`, `cell_1.parquet`, etc.). This guarantees that the next step will find the files it needs with the correct names.
*   **Step 6 (Compile)**: It calls `vpt compile-tile-segmentation`, giving it the **filtered ROI specification file**. Because this file only lists the tiles that were processed, and the files have been renamed correctly, this step now succeeds.
*   **Step 7-8 (Analysis & VZG Update)**: It runs the final analysis steps on the compiled ROI data, producing the final `_roi_resegmented.vzg2` file.

### Diagnostic Tools

During our debugging, we also created `verify_segmentation_location.py`. This script is a valuable tool for plotting the output `cellpose_micron_space.parquet` file to visually confirm that the segmented cells are in the correct location.

## How to Use

1.  **Edit `roi_coords.json`**: Update the file with the micron coordinates of the region you want to segment.
2.  **Configure the Script**: Open `run_roi_segmentation_custom.sh` and ensure the `FOLDER`, `REGION`, `DATA_PATH`, and `CONFIG_NAME` variables are set correctly for your experiment.
3.  **Run the Script**: Execute the script from your terminal:
    ```bash
    bash Segmentation/resegment_ROI_workflow/run_roi_segmentation_custom.sh