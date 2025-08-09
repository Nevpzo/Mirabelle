# Mirabelle
Interactive tool to manually define a polygonal mask and store it in netCDF format.

<p align="center"> <img src="/img/overview.png" align="center"> </p>

When working with climate data, it is often essential to create masks that isolate specific regions for targeted analysis. While pre-made masks sometimes exist, adapting them to match a given dataset’s grid and format can be time-consuming. This tool provides a fast, interactive way to draw custom polygonal masks directly on a map and export them in netCDF format, ready for use in climate data processing and analysis. While this version is centered on Greenland, it can be easily adapted to any region by providing a reference netCDF file with a lon/lat grid and changing the map extent.

The tool supports zooming, panning, undoing points, clearing the view, and displaying the computed area of each selection.

# Requirements
The program requires:
- `numpy`
- `matplotlib`
- `tkinter`
- `cartopy`
- `netCDF4`
- `shapely`
- `pyproj`

Install missing dependencies with:
```
pip install numpy matplotlib cartopy netCDF4 shapely pyproj
```

# Usage

The interface will display a stereographic projection of Greenland. You can then define polygons with the mouse:

**Mouse controls:**
- **Left Click**: Place a vertex of the polygon. Vertices are connected automatically.
- **Right Click**: Close the polygon (connect last vertex to the first).
- **Middle Click + Drag**: Pan the map.
- **Scroll Wheel**: Zoom in and out.
- **Ctrl + Right Click**: Remove a polygon.

**Buttons:**
- **Undo Last Point**: Removes the most recent vertex from the active polygon.
- **Reset View**: Restores the original zoom and position.
- **Clear View**: Removes all polygons and resets the map.
- **Save As**: Exports the current mask to a `.nc` file.
- **Clear on save** *(checkbox)*: Clears all polygons after saving.
- **Display area** *(checkbox)*: Shows polygon area in m² on the map.

**Keyboard shortcuts:**
- `Ctrl+S` → Save
- `Ctrl+Z` → Undo last point
- `Ctrl+R` → Reset view
- `Ctrl+Q` → Quit
- `Ctrl+Backspace` → Clear view

# Notes
- Make sure the path to the reference netCDF file is correct before running.  
- Multiple polygons can be drawn before saving. They will be combined into a single mask.
