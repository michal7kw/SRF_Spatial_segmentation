# Cell Segmentation Visualization Tool

This tool provides comprehensive visualization capabilities for cell segmentation results from the VPT (Vizgen Post-processing Tool) pipeline.

## Overview

The segmentation pipeline has processed the ROI and detected **56 cells** across **12 tiles**, with the following distribution:
- **Tile 4**: 49 cells (primary detection area)
- **Tile 8**: 7 cells (secondary detection area)
- **Other tiles**: 0 cells each

## Files Created

1. **`visualize_segmentation_results.py`** - Full-featured visualization script
2. **`visualize_segmentation_simple.py`** - Simplified, memory-efficient version
3. **`test_visualization.py`** - Test script for full version
4. **`test_simple_visualization.py`** - Test script for simplified version
5. **`README_visualization.md`** - This documentation file

## Requirements

The script requires the following Python packages:
```bash
pip install numpy pandas matplotlib seaborn scikit-image geopandas shapely
```

These should already be available in your `vpt` conda environment.

## Recommended Usage

**For most users, start with the simplified version:**

### 1. Quick Summary Report
```bash
python visualize_segmentation_simple.py --summary-only
```

### 2. Visualize Tiles with Cells
```bash
# View tile 4 (contains 49 cells) - automatically saves to ./segmentation_visualizations/
python visualize_segmentation_simple.py --tile-index 4

# View tile 8 (contains 7 cells) - automatically saves to ./segmentation_visualizations/
python visualize_segmentation_simple.py --tile-index 8

# View all tiles with cells - saves all to ./segmentation_visualizations/
python visualize_segmentation_simple.py --all-tiles
```

### 3. Custom Save Location
```bash
# Save visualizations to a specific directory
python visualize_segmentation_simple.py --tile-index 4 --save-dir ./my_output
```

### 4. Interactive Mode
```bash
# Interactive selection of tiles
python visualize_segmentation_simple.py
```

## Advanced Usage (Full Version)

**Note: The full version may have memory issues with large images. Use the simplified version if you encounter problems.**

### 1. ROI Overview Visualization
```bash
python visualize_segmentation_results.py
```
This shows the entire ROI with all detected cell boundaries overlaid on DAPI and PolyT channels.

### 2. Custom Channel Visualization
```bash
# Visualize with different channels
python visualize_segmentation_results.py --channels DAPI Pcp4 --z-layer 6
```

## Command Line Options

### Simplified Version (`visualize_segmentation_simple.py`)

| Option | Description | Default |
|--------|-------------|---------|
| `--output-dir` | Directory with segmentation results | `1task_cellpose2_p30-E165_R1_roi_analysis_parallel` |
| `--data-dir` | Directory with original images | `DATA/p30-E165/R1` |
| `--tile-index` | Specific tile to visualize | None (interactive mode) |
| `--channels` | Image channels to display | `DAPI PolyT` |
| `--z-layer` | Z-layer to visualize | `6` |
| `--save-dir` | Directory to save images | None (display only) |
| `--summary-only` | Only show summary report | False |
| `--all-tiles` | Visualize all tiles with cells | False |

### Full Version (`visualize_segmentation_results.py`)

| Option | Description | Default |
|--------|-------------|---------|
| `--output-dir` | Directory with segmentation results | `1task_cellpose2_p30-E165_R1_roi_analysis_parallel` |
| `--data-dir` | Directory with original images | `DATA/p30-E165/R1` |
| `--roi-coords` | ROI coordinates JSON file | `roi_coords.json` |
| `--tile-index` | Specific tile to visualize | None (shows ROI overview) |
| `--channels` | Image channels to display | `DAPI PolyT` |
| `--z-layer` | Z-layer to visualize | `6` |
| `--save-dir` | Directory to save images | None (display only) |
| `--summary-only` | Only show summary report | False |

## Visualization Features

### ROI Overview
- Shows the entire region of interest
- Overlays all detected cell boundaries
- Displays multiple channels (DAPI, PolyT) and RGB overlay
- Provides summary statistics

### Tile-Specific Views
- Focuses on individual tiles for detailed inspection
- Shows original images with segmentation overlays
- Useful for quality assessment of segmentation results
- Highlights areas with high/low cell density

### Interactive Mode
When running without `--tile-index`, the script will:
1. Show ROI overview first
2. Ask if you want to see individual tiles
3. Allow selection of specific tiles or all tiles

## Segmentation Results Summary

Based on your pipeline run:

```
ROI Coordinates:
  X: 10691.9 - 11990.2 microns
  Y: 6190.1 - 6940.4 microns

Tile Information:
  Total tiles: 12
  Tiles with cells: 2
  
Cell Detection:
  Total cells: 56
  Average per tile: 4.7 ± 13.5
  Range: 0-49 cells per tile
```

## Quality Assessment Tips

1. **Check High-Density Areas**: Focus on tile 4 which has 49 cells
   ```bash
   python visualize_segmentation_results.py --tile-index 4
   ```

2. **Verify Boundary Accuracy**: Look for:
   - Proper cell boundary detection
   - No over-segmentation (cells split incorrectly)
   - No under-segmentation (multiple cells merged)

3. **Assess Empty Tiles**: Check why 10 tiles have 0 cells:
   ```bash
   python visualize_segmentation_results.py --tile-index 0
   ```

4. **Compare Channels**: Ensure segmentation uses both DAPI and PolyT appropriately

## Troubleshooting

### Common Issues

1. **Memory Issues / Script Killed**:
   - Use the simplified version: `python visualize_segmentation_simple.py`
   - Focus on individual tiles instead of full ROI
   - The simplified version only loads tiles with detected cells

2. **"'Series' object has no attribute 'geometry'" Error**:
   - This is fixed in the updated versions
   - Use the simplified version for better error handling

3. **Import Errors**:
   - Ensure all required packages are installed in the `vpt` environment
   - Run: `pip install numpy pandas matplotlib scikit-image geopandas shapely`

4. **File Not Found**:
   - Verify that segmentation results exist in the expected directory
   - Check that the output directory name matches your pipeline run

### Getting Help

Run the test scripts to verify everything is working:
```bash
# Test simplified version (recommended)
python test_simple_visualization.py

# Test full version
python test_visualization.py
```

### Performance Tips

1. **Start with the simplified version** - it's more memory efficient
2. **Focus on tiles with cells** - tiles 4 and 8 in your case
3. **Use `--summary-only` first** to check results without loading images
4. **Save visualizations** instead of displaying them to reduce memory usage

## Output Files

When using `--save-dir`, the following files are created:
- `roi_overview_segmentation.png` - Full ROI visualization
- `tile_X_segmentation.png` - Individual tile visualizations

## Technical Details

- **Coordinate System**: Uses mosaic pixel coordinates
- **Image Enhancement**: Applies CLAHE normalization and contrast stretching
- **Segmentation Format**: Reads parquet files with polygon geometries
- **Visualization**: Uses matplotlib with customizable color schemes

## Next Steps

After visual inspection, you may want to:
1. Adjust segmentation parameters if results are unsatisfactory
2. Re-run the pipeline with modified settings
3. Proceed with downstream analysis using the detected cells
4. Export results for further analysis in other tools

---

*Generated for the SRF Spatial Segmentation Pipeline*