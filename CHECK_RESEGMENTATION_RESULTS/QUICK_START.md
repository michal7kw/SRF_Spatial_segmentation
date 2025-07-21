# Quick Start Guide - Cell Segmentation Visualization

## The Issue You Encountered

When you ran either:
```bash
python visualize_segmentation_simple.py --tile-index 4
python visualize_segmentation_results.py
```

The scripts processed successfully but you didn't see where the visualizations were saved. This is because:
1. You're running on a cluster without a display
2. The original versions tried to show plots with `plt.show()` but had no display
3. No `--save-dir` was specified, so nothing was saved

## Fixed Version

**Both scripts now automatically save all visualizations** to avoid this issue.

## Where Your Files Are Saved

### Default Location
When you run the script without `--save-dir`, files are automatically saved to:
```
./segmentation_visualizations/
```

### Your Results
After running:
```bash
python visualize_segmentation_simple.py --tile-index 4
```

You should find:
```
./segmentation_visualizations/tile_4_segmentation.png
```

## Quick Commands

### Simplified Version (Recommended)
```bash
# View tile 4 (49 cells) - saves to ./segmentation_visualizations/tile_4_segmentation.png
python visualize_segmentation_simple.py --tile-index 4

# View tile 8 (7 cells) - saves to ./segmentation_visualizations/tile_8_segmentation.png
python visualize_segmentation_simple.py --tile-index 8

# View both tiles - saves both files
python visualize_segmentation_simple.py --all-tiles
```

### Full Version (Now Fixed)
```bash
# ROI overview - saves to ./segmentation_visualizations/roi_overview_segmentation.png
python visualize_segmentation_results.py

# Specific tile - saves to ./segmentation_visualizations/tile_4_segmentation.png
python visualize_segmentation_results.py --tile-index 4
```

### Custom Save Location (Both Versions)
```bash
# Save to custom location
python visualize_segmentation_simple.py --tile-index 4 --save-dir ./my_results
python visualize_segmentation_results.py --tile-index 4 --save-dir ./my_results
```

## Check Your Results

```bash
# List generated files
ls -la segmentation_visualizations/

# Copy files to view locally (if needed)
scp user@cluster:path/to/segmentation_visualizations/*.png ./local_folder/
```

## What You'll See in the Visualizations

Each PNG file shows:
- **Left panel**: DAPI channel with red cell boundaries overlaid
- **Right panel**: PolyT channel with red cell boundaries overlaid
- **Title**: Shows tile number and cell count
- **High resolution**: 100 DPI for detailed inspection

The visualizations will help you assess:
- Cell boundary accuracy
- Segmentation quality
- Channel alignment
- Detection completeness

## Next Steps

1. **Check the generated PNG files** in `./segmentation_visualizations/`
2. **Review tile 4** (49 cells) for segmentation quality
3. **Review tile 8** (7 cells) for comparison
4. **Adjust segmentation parameters** if needed based on visual inspection