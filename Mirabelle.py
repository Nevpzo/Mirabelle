"""
Last update: 17/08/2025
Version: Mirabelle 1.4

Description:
    - Interactive mapping tool to manually define polygonal masks
    and export them in netCDF format, ready for use with climate
    and weather models.

Update Note:
    - Removed polygon counter at the bottom of the screen
    - Added reference netCDF loading

Contact:
    - nnevpzo@gmail.com
    - GitHub: https://github.com/Nevpzo/Mirabelle/
    - Zenodo: https://doi.org/10.5281/zenodo.16789952
"""

# =================================================
# Imports
import tkinter as tk
from tkinter import messagebox
from tkinter import filedialog

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.patches import Polygon
from matplotlib.patches import Rectangle

import cartopy.crs as ccrs
import cartopy.feature as cfeature

from netCDF4 import Dataset
from shapely import vectorized
import shapely.geometry as shp
from shapely.ops import transform, unary_union
import pyproj

from datetime import datetime

# =================================================
# Initialize polygon, panning and map state

polygons = []
plotgons = []
pointsX = []
pointsY = []
polygon_color = "#ff595e"
is_poly_closed = False
panning = False
pan_start = None
map_extent = [-100, 100, -100, 100]

# =================================================
# Functions

def quit_me(e=None):
    global pointsX, polygons
    """Kills the Tkinter window"""
    if pointsX or polygons:
        answer = messagebox.askyesno("Quit", "You have unsaved changes. Quit anyways ?")
        if answer:
            root.quit()
            root.destroy()
    else:
        root.quit()
        root.destroy()


def polygon_to_netCDF(polygons, outfile):
    """
    Save user polygons as a netCDF mask file.
    Polygons are in lon/lat (PlateCarree), grid is transformed to projection coords
    before containment test.
    """

    global lon_grid, lat_grid, projection

    if not polygons:
        raise ValueError("No polygons to save.")

    if lon_grid is None or lat_grid is None:
        raise RuntimeError("lon_grid / lat_grid not loaded (globals).")

    if lon_grid.shape != lat_grid.shape:
        raise RuntimeError(
            f"lon_grid shape {lon_grid.shape} != lat_grid shape {lat_grid.shape}"
        )

    # Transform each polygon from lon/lat to projection coordinates
    projected_polys = []
    for poly in polygons:
        xs, ys = projection.transform_points(
            ccrs.PlateCarree(),
            np.array(poly.exterior.coords)[:, 0],
            np.array(poly.exterior.coords)[:, 1]
        )[:, :2].T
        projected_polys.append(type(poly)(zip(xs, ys)))

    # Merge transformed polygons into a single geometry
    combined_poly_proj = unary_union(projected_polys)

    # Transform the grid to projection coordinates
    proj_points = projection.transform_points(
        ccrs.PlateCarree(), lon_grid, lat_grid
    )
    x_proj, y_proj = proj_points[..., 0], proj_points[..., 1]

    # Containment test in projection space
    mask_bool = vectorized.contains(combined_poly_proj, x_proj, y_proj)

    if mask_bool.shape != lon_grid.shape:
        raise RuntimeError(
            f"Mask shape {mask_bool.shape} != lon/lat grid shape {lon_grid.shape}"
        )

    ny, nx = lon_grid.shape

    with Dataset(outfile, "w", format="NETCDF4") as root_grp:
        root_grp.description = "Mask of user-selected region"
        root_grp.source = "Generated from Mirabelle polygon selection tool"
        root_grp.history = f"Created on {datetime.now().isoformat()}"

        root_grp.createDimension("y", ny)
        root_grp.createDimension("x", nx)

        lon_var = root_grp.createVariable("lon", "f4", ("y", "x"))
        lat_var = root_grp.createVariable("lat", "f4", ("y", "x"))
        mask_var = root_grp.createVariable("mask", "i1", ("y", "x"))

        lon_var[:, :] = lon_grid
        lat_var[:, :] = lat_grid
        mask_var[:, :] = mask_bool.astype("i1")

        lon_var.units = "degrees_east"
        lon_var.standard_name = "longitude"
        lon_var.long_name = "Longitude"

        lat_var.units = "degrees_north"
        lat_var.standard_name = "latitude"
        lat_var.long_name = "Latitude"

        mask_var.units = "1"
        mask_var.standard_name = "land_mask"
        mask_var.long_name = "User-defined mask"

    print(f"Mask stored in {outfile.split('/')[-1]}")


def on_click(event):
    """Handles left (add vertex), right (close polygon), and middle (start panning) clicks"""
    global is_poly_closed, pointsX, pointsY, panning, pan_start, polygon_color, polygons, plotgons, areaTicked
    if event.inaxes is None:
        return

    if (event.button == 1) and (not(is_poly_closed)): # Left click: add vertex
            if len(pointsX) > 0:
                # Draw line from last point to current one
                ax.plot([pointsX[-1], event.xdata], [pointsY[-1], event.ydata], color=polygon_color)

            pointsX.append(event.xdata)
            pointsY.append(event.ydata)
            ax.plot(event.xdata, event.ydata, color=polygon_color, marker="+")
            canvas.draw()

    elif event.button == 3: # Right click
        if event.key == "control": # Ctrl click to remove polygon
            clean_map()
            for i in reversed(range(len(plotgons))):
                poly = plotgons[i]
                if poly.contains_point([event.x, event.y]):
                    if i < len(polygons):
                        polygons.pop(i)
                    plotgons.pop(i)
                    poly.remove()
                    canvas.draw()
                    break

        elif not (is_poly_closed) and (len(pointsX) > 2): # Normal click to close polygon
            # Fill the shape on the map
            maskedZone = list(zip(pointsX, pointsY))
            plotgon = Polygon(
                maskedZone, closed=True, facecolor=polygon_color, edgecolor="black", alpha=0.7
            )
        
            # Convert points to lon/lat
            polygon = shp.Polygon([
                ccrs.PlateCarree().transform_point(x, y, projection)
                for x, y in zip(pointsX, pointsY)
                ])

            colors = ["#ff595e", "#ffa165", "#ffca3a", "#83b731", "#1982C4", "#6A4C93"]
            polygon_color = colors[(colors.index(polygon_color) + 1) % 6]
            pointsX = []
            pointsY = []

            plotgons.append(plotgon)
            polygons.append(polygon)

            clean_map()

    elif event.button == 2: # Middle click: start panning
        panning = True
        pan_start = (event.x, event.y)


def on_button_release(event):
    """Ends panning mode when middle button released"""
    global panning
    if event.button == 2:
        panning = False


def on_motion(event):
    """Updates coordinates in window title; handles real-time panning"""
    global pan_start, panning
    if event.inaxes:
        lon, lat = ccrs.PlateCarree().transform_point(
            event.xdata, event.ydata, src_crs=projection
        )
        root.title(f"Mirabelle - Lat: {lat:.2f}, Lon: {lon:.2f}")

    if panning and event.inaxes:
        dx = event.x - pan_start[0]
        dy = event.y - pan_start[1]

        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()
        width = x1 - x0
        height = y1 - y0

        ax.set_xlim(
            x0 - dx * width / canvas.get_width_height()[0],
            x1 - dx * width / canvas.get_width_height()[0],
        )
        ax.set_ylim(
            y0 - dy * height / canvas.get_width_height()[1],
            y1 - dy * height / canvas.get_width_height()[1],
        )

        pan_start = (event.x, event.y)
        canvas.draw()


def on_scroll(event):
    """Zoom in/out around cursor"""
    if event.inaxes is None:
        return

    base_scale = 1.2
    scale_factor = base_scale if event.button == "down" else 1 / base_scale

    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    xdata = event.xdata
    ydata = event.ydata

    new_xlim = [
        xdata - (xdata - xlim[0]) * scale_factor,
        xdata + (xlim[1] - xdata) * scale_factor,
    ]
    new_ylim = [
        ydata - (ydata - ylim[0]) * scale_factor,
        ydata + (ylim[1] - ydata) * scale_factor,
    ]

    ax.set_xlim(new_xlim)
    ax.set_ylim(new_ylim)
    canvas.draw()


def undo_last_point(e=None):
    """Removes last added vertex from the active (unfinished) polygon only."""
    global pointsX, pointsY, polygon_color, plotgons, is_poly_closed

    if len(pointsX) == 0:
        return

    # Remove the last point
    pointsX.pop()
    pointsY.pop()

    # Redraw map
    clean_map()


def clean_map():
    global pointsX, pointsY, plotgons, polygons
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    ax.clear()
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS, linestyle=":")
    ax.gridlines(draw_labels=True)

    # Redraw completed (closed) polygons
    for i, poly in enumerate(plotgons):
        ax.add_patch(poly)

    # Redraw active (unfinished) polygon
    for i in range(len(pointsX)):
        ax.plot(pointsX[i], pointsY[i], color=polygon_color, marker="+")
        if i > 0:
            ax.plot([pointsX[i - 1], pointsX[i]], [pointsY[i - 1], pointsY[i]], color=polygon_color)
    canvas.draw()


def reset_view(e=None):
    """Restores initial map extent"""
    global map_extent
    ax.set_extent(map_extent, crs=ccrs.PlateCarree())
    canvas.draw()


def clear_view(e=None):
    """Resets map and clears all user input"""
    global plotgons, polygons, pointsX, pointsY, is_poly_closed, polygon_color, map_extent

    if e is None:
        answer = messagebox.askyesno(title='Confirmation',
                        message='Are you sure that you want to clear?')
    else:
        answer = True

    if answer:
        plotgons.clear()
        polygons.clear()
        pointsX.clear()
        pointsY.clear()
        is_poly_closed = False
        polygon_color = "#ff595e"

        ax.clear()
        ax.set_extent(map_extent, crs=ccrs.PlateCarree())
        ax.coastlines()
        ax.add_feature(cfeature.BORDERS, linestyle=":")
        ax.gridlines(draw_labels=True)
        canvas.draw()


def display_area():
    global plotgons
    if areaTicked.get() == 0:
        for txt in ax.texts:
            txt.remove()
        canvas.draw()
        return

    for i, poly in enumerate(plotgons):
        transformer = pyproj.Transformer.from_crs("EPSG:4326", pyproj.CRS("EPSG:3413"), always_xy=True).transform
        polygon_projected = transform(transformer, polygons[i])
        area = polygon_projected.area

        verts = poly.get_xy()
        x_center = verts[:, 0].mean()
        y_center = verts[:, 1].mean()

        ax.text(
            x_center, y_center,
            f"{area:.2e} m²",
            color='white', fontsize=7, weight='bold',
            ha='center', va='center',
            bbox=dict(facecolor='black', alpha=0.6, boxstyle='round,pad=0.2')
        )
    canvas.draw()


def save_as(e=None):
    """Save current polygon to user-specified netCDF file"""
    global is_poly_closed, pointsX, pointsY, polygons
    if len(polygons) == 0:
        messagebox.showwarning("Save Error", "No closed polygon to save.")
        return

    now = datetime.now()
    file_path = filedialog.asksaveasfilename(
        defaultextension=".nc",
        filetypes=[("NetCDF files", "*.nc")],
        title="Save mask as...",
        initialfile = f"mask-{now.year}{now.month}{now.day}_{now.hour}{now.minute}.nc"
    )
    if file_path:
        polygon_to_netCDF(polygons, file_path)
    else:
        return
    
    if clearTicked.get() == 1: # Clear if "clear on save" checkbox is ticked
        polygons.clear()
        plotgons.clear()
        clear_view(e="")
    else:                   # Reduce polygon opacity otherwise
        for plotgon in plotgons:
            plotgon.set_alpha(0.3)
        polygons.clear()
        canvas.draw()
        

def choose_projection(lon_min, lon_max, lat_min, lat_max):
    if lat_min >= 45:
        lat_center = (lat_min + lat_max) / 2
        lon_center = (lon_min + lon_max) / 2
        return ccrs.Stereographic(central_latitude=lat_center, central_longitude=lon_center)
    elif lat_max <= -45:
        lat_center = (lat_min + lat_max) / 2
        lon_center = (lon_min + lon_max) / 2
        return ccrs.SouthPolarStereo(central_longitude=lon_center)
    else:
        return ccrs.PlateCarree()


def get_var(names, nc_data):
    variables = nc_data.variables
    for name in names:
        if name in variables:
            return nc_data.variables[name][:]
    raise KeyError(f"No variable {names} found.")


def load_netcdf():
    global ax, fig, canvas, map_extent, lon_grid, lat_grid, projection

    clear_view(e="")

    file_path = filedialog.askopenfilename(
        title="Select base NetCDF file",
        filetypes=[("NetCDF files", "*.nc *.nc4")]
    )
    if not file_path:
        return

    try:
        nc_data = Dataset(file_path, "r")

        lon = get_var(["lon", "longitude", "long", "x"], nc_data)
        lat = get_var(["lat", "latitude", "y"], nc_data)

        lon_grid = np.array(lon)
        lat_grid = np.array(lat)

        if lon_grid.ndim == 1 and lat_grid.ndim == 1:
            lon_grid, lat_grid = np.meshgrid(lon, lat)

        lon_min, lon_max = np.min(lon), np.max(lon)
        lat_min, lat_max = np.min(lat), np.max(lat)

        projection = choose_projection(lon_min, lon_max, lat_min, lat_max)
        map_extent = [lon_min, lon_max, lat_min, lat_max]

        fig.clear()
        ax = fig.add_subplot(projection=projection)
        ax.set_extent(map_extent, crs=ccrs.PlateCarree())
        ax.coastlines()
        ax.add_feature(cfeature.BORDERS, linestyle=":")
        ax.gridlines(draw_labels=True)

        canvas.draw()

    except Exception as e:
        messagebox.showerror("Error", f"Unable to load file: {e}")


def show_help():
    help_win = tk.Toplevel(root)
    help_win.title("Help")
    help_win.geometry("500x400")
    help_win.resizable(False, False)

    button_frame = tk.Frame(help_win)
    button_frame.pack(side=tk.TOP, fill=tk.X, pady=5)

    text_area = tk.Text(help_win, wrap=tk.WORD, font=("Arial", 14))
    text_area.pack(expand=True, fill=tk.BOTH, padx=10, pady=10)

    def update_text(content):
        text_area.config(state=tk.NORMAL)
        text_area.delete(1.0, tk.END)
        text_area.insert(tk.END, content)
        text_area.config(state=tk.DISABLED)

    def show_mouse():
        update_text(
            "Mouse Controls:\n"
            "- Left Click: Place a point on the map. If two or more points are placed,\n"
            "  lines will automatically connect them.\n"
            "- Right Click: Close the polygon by connecting the last point to the first.\n"
            "- Middle Click: Click and drag to move the map around.\n"
            "- Scroll Wheel: Zoom in and out on the map.\n"
            "- Ctrl + Right Click: Removes the clicked polygon."
        )

    def show_buttons():
        update_text(
            "Buttons:\n"
            "- Load netCDF: Load the reference netCDF file to mask a region on.\n"
            "- Undo: Remove the last point and its connecting line.\n"
            "- Reset View: Reset the map to its original zoom level and position.\n"
            "- Clear View: Remove all points, lines, and polygons, and reset the view.\n"
            "- Save As: Saves the selection to a netCDF file.\n"
            "- Clear on save: If checked, saving the selection clears the map. If uncheck: reduces the opacity of drawn polygons.\n"
            "- Display area: If checked, displays the area of currently drawn polygons."
        )

    def show_shortcuts():
        update_text(
            "Shortcuts:\n"
            "- Ctrl+S: Save\n"
            "- Ctrl+Z: Undo\n"
            "- Ctrl+R: Reset view\n"
            "- Ctrl+Q: Close app\n"
            "- Ctrl+BackSpace: Clear view"
        )

    # Buttons to switch help categories
    tk.Button(button_frame, text="Mouse Controls", command=show_mouse).pack(side=tk.LEFT, padx=5)
    tk.Button(button_frame, text="Buttons", command=show_buttons).pack(side=tk.LEFT, padx=5)
    tk.Button(button_frame, text="Shortcuts", command=show_shortcuts).pack(side=tk.LEFT, padx=5)

    # Show mouse help by default
    show_mouse()

# =================================================
# Initialize Tkinter Window and Canvas

root = tk.Tk()
root.title("Mirabelle")
root.geometry("600x800")
root.protocol("WM_DELETE_WINDOW", quit_me)

# Projection and map setup
projection = ccrs.PlateCarree()
fig, ax = plt.subplots(facecolor="w", subplot_kw=dict(projection=projection))

ax.set_extent([-100, 100, -100, 100], crs=ccrs.PlateCarree())
ax.coastlines()
ax.add_feature(cfeature.BORDERS, linestyle=":")
ax.gridlines(draw_labels=True)

# =================================================
# Features

# Button controls
button_frame = tk.Frame(root)
button_frame.pack(side=tk.TOP, fill=tk.X, padx=5, pady=5)

clear_frame = tk.Frame(root)
clear_frame.pack(side=tk.TOP, fill=tk.X, padx=5, pady=5)

clear_btn = tk.Button(button_frame, text="Clear", command=clear_view, height = 1, width = 12, font=('Helvetica', 11))
clear_btn.pack(side=tk.LEFT, padx=2)

reset_btn = tk.Button(button_frame, text="Reset View", command=reset_view, height = 1, width = 12, font=('Helvetica', 11))
reset_btn.pack(side=tk.LEFT, padx=2)

undo_btn = tk.Button(button_frame, text="Undo Last Point", command=undo_last_point, height = 1, width = 12, font=('Helvetica', 11))
undo_btn.pack(side=tk.LEFT, padx=2)

help_btn = tk.Button(button_frame, text="Help", command=show_help, height=1, width=12, font=('Helvetica', 11))
help_btn.pack(side=tk.RIGHT, padx=2)

save_btn = tk.Button(button_frame, text="Save As", command=save_as, height=1, width=12, font=('Helvetica', 11))
save_btn.pack(side=tk.RIGHT, padx=2)

clearTicked = tk.IntVar()
clearChk = tk.Checkbutton(clear_frame, text='Clear on save', variable=clearTicked, font=('Helvetica', 12))
clearChk.pack(side=tk.LEFT)

areaTicked = tk.IntVar()
areaChk = tk.Checkbutton(clear_frame, text='Display area', variable=areaTicked, font=('Helvetica', 12), command=lambda: display_area())
areaChk.pack(side=tk.RIGHT)

baseCDFBtn = tk.Button(root, text='Load NetCDF', command=load_netcdf, height=1, width=12, font=('Helvetica', 11))
baseCDFBtn.pack(side=tk.BOTTOM)

# Mouse event and shortcuts bindings
canvas = FigureCanvasTkAgg(fig, master=root)
canvas.draw()
canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
canvas.mpl_connect("button_press_event", on_click)
canvas.mpl_connect("button_release_event", on_button_release)
canvas.mpl_connect("scroll_event", on_scroll)
canvas.mpl_connect("motion_notify_event", on_motion)
root.bind('<Control-s>', save_as)
root.bind('<Control-z>', undo_last_point)
root.bind('<Control-BackSpace>', clear_view)
root.bind('<Control-r>', reset_view)
root.bind('<Control-q>', quit_me)

# =================================================
# Mainloop


root.mainloop()
