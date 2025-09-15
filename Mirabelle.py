"""
Last update: 09/09/2025
Version: Mirabelle 1.6

Description:
    - Interactive mapping tool to draw polygonal masks interactively
    and export them as NetCDF for geophysical models or GeoJSON for GIS use.

Update Note:
    - Added option to export polygons as geojson file
    - Changed layout to wrap buttons in a menubar

Contact:
    - nnevpzo@gmail.com
    - GitHub: https://github.com/Nevpzo/Mirabelle/
    - Zenodo: https://doi.org/10.5281/zenodo.16789951
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

import cartopy.crs as ccrs
import cartopy.feature as cfeature

from netCDF4 import Dataset
from shapely import vectorized
import shapely.geometry as shp
from shapely.ops import transform, unary_union
import pyproj

from datetime import datetime
import json

# =================================================
# Utility functions

def polygon_to_netCDF(polygons, outfile, lon_grid, lat_grid, projection):
    """
    Save user polygons as a netCDF mask file.
    Polygons are in lon/lat (PlateCarree), grid is transformed to projection coords
    before containment test.
    """
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

        root_grp.createDimension("ny", ny)
        root_grp.createDimension("nx", nx)

        lon_var = root_grp.createVariable("lon", "f4", ("ny", "nx"))
        lat_var = root_grp.createVariable("lat", "f4", ("ny", "nx"))
        mask_var = root_grp.createVariable("mask", "i1", ("ny", "nx"))

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


def polygon_to_geojson(polygons, outfile):
    feature_collection = {
        "type": "FeatureCollection",
        "features": []
    }
    for poly in polygons:
        feature = {
            "type": "Feature",
            "properties": {},
            "geometry": shp.mapping(poly)
        }
        feature_collection["features"].append(feature)

    with open(outfile, "w") as f:
        json.dump(feature_collection, f, indent=2)


def choose_projection(lon_min, lon_max, lat_min, lat_max):
    """
    Choose a projection based on reference netCDF extent
    """
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

# =================================================
# App

class MirabelleApp:
    # -------------------
    # Initialisation
    # -------------------
    def __init__(self, root):
        # Store Tk root
        self.root = root
        self.root.title("Mirabelle")
        self.root.geometry("600x800")
        self.root.protocol("WM_DELETE_WINDOW", self.quit_me)
        self.root.option_add("*Menu*Font", "Helvetica 12")

        # Initial state
        self.polygons = []
        self.plotgons = []
        self.pointsX = []
        self.pointsY = []
        self.polygon_color = "#ff595e"
        self.is_poly_closed = False
        self.panning = False
        self.pan_start = None
        self.map_extent = [-100, 100, -100, 100]
        self.projection = ccrs.PlateCarree()
        self.lon_grid = None
        self.lat_grid = None

        # Setup buttons and menus
        self._build_ui()

        # Setup matplotlib figure
        self.fig, self.ax = plt.subplots(
            facecolor="w",
            subplot_kw=dict(projection=self.projection)
        )
        self.ax.set_extent(self.map_extent, crs=ccrs.PlateCarree())
        self.ax.coastlines()
        self.ax.add_feature(cfeature.BORDERS, linestyle=":")
        self.ax.gridlines(draw_labels=True)

        # Embed in Tk
        self.canvas = FigureCanvasTkAgg(self.fig, master=root)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Bind mouse events
        self.canvas.mpl_connect("button_press_event", self.on_click)
        self.canvas.mpl_connect("button_release_event", self.on_button_release)
        self.canvas.mpl_connect("scroll_event", self.on_scroll)
        self.canvas.mpl_connect("motion_notify_event", self.on_motion)


    def _build_ui(self):
        """Create menubar and checkboxes."""
        # -------------------
        # Menubar
        # -------------------
        menubar = tk.Menu(self.root)

        # File menu
        file_menu = tk.Menu(menubar, tearoff=0)
        file_menu.add_command(label="Load NetCDF", command=self.load_netcdf)

        # Save As submenu
        save_menu = tk.Menu(file_menu, tearoff=0)
        save_menu.add_command(label="NetCDF", command=self.save_as_netcdf)
        save_menu.add_command(label="GeoJSON", command=self.save_as_geojson)
        file_menu.add_cascade(label="Save As", menu=save_menu)

        file_menu.add_separator()
        file_menu.add_command(label="Help", command=self.show_help)
        file_menu.add_separator()
        file_menu.add_command(label="Quit", command=self.quit_me)
        menubar.add_cascade(label="File", menu=file_menu)

        # Edit menu
        menubar.add_command(label="Clear View", command=self.clear_view)
        menubar.add_command(label="Reset View", command=self.reset_view)
        menubar.add_command(label="Undo", command=self.undo_last_point)

        # Attach menubar to root
        self.root.config(menu=menubar)


        # -------------------
        # Checkboxes
        # -------------------
        clear_frame = tk.Frame(self.root)
        clear_frame.pack(side=tk.TOP, fill=tk.X, padx=5, pady=5)

        self.clearTicked = tk.IntVar()
        clearChk = tk.Checkbutton(clear_frame, text='Clear on save',
                                variable=self.clearTicked, font=('Helvetica', 12))
        clearChk.pack(side=tk.LEFT)

        self.areaTicked = tk.IntVar()
        areaChk = tk.Checkbutton(clear_frame, text='Display area',
                                variable=self.areaTicked, font=('Helvetica', 12),
                                command=self.display_area)
        areaChk.pack(side=tk.RIGHT)


        # -------------------
        # Shortcuts
        # -------------------
        self.root.bind('<Control-s>', self.save_as_netcdf)
        self.root.bind('<Control-Shift-S>', self.save_as_geojson)
        self.root.bind('<Control-z>', self.undo_last_point)
        self.root.bind('<Control-BackSpace>', self.clear_view)
        self.root.bind('<Control-r>', self.reset_view)
        self.root.bind('<Control-q>', self.quit_me)

    def quit_me(self, e=None):
        """Kills the Tkinter window"""
        if self.pointsX or self.polygons:
            answer = messagebox.askyesno("Quit", "You have unsaved changes. Quit anyways ?")
            if answer:
                self.root.quit()
                self.root.destroy()
        else:
            self.root.quit()
            self.root.destroy()


    # -------------------
    # Data handling
    # -------------------
    def load_netcdf(self):
        self.clear_view(e="")

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

            self.lon_grid = np.array(lon)
            self.lat_grid = np.array(lat)

            if self.lon_grid.ndim == 1 and self.lat_grid.ndim == 1:
                self.lon_grid, self.lat_grid = np.meshgrid(lon, lat)

            lon_min, lon_max = np.min(lon), np.max(lon)
            lat_min, lat_max = np.min(lat), np.max(lat)

            self.projection = choose_projection(lon_min, lon_max, lat_min, lat_max)
            self.map_extent = [lon_min, lon_max, lat_min, lat_max]

            self.fig.clear()
            self.ax = self.fig.add_subplot(projection=self.projection)
            self.ax.set_extent(self.map_extent, crs=ccrs.PlateCarree())
            self.ax.coastlines()
            self.ax.add_feature(cfeature.BORDERS, linestyle=":")
            self.ax.gridlines(draw_labels=True)

            self.canvas.draw()

        except Exception as e:
            messagebox.showerror("Error", f"Unable to load file: {e}")


    def save_as_netcdf(self, e=None):
        """Save current polygon to user-specified netCDF file"""
        if len(self.polygons) == 0:
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
            polygon_to_netCDF(self.polygons, file_path, self.lon_grid, self.lat_grid, self.projection)
        else:
            return

        if self.clearTicked.get() == 1: # Clear if "clear on save" checkbox is ticked
            self.polygons.clear()
            self.plotgons.clear()
            self.clear_view(e="")
        else:                           # Reduce polygon opacity otherwise
            for plotgon in self.plotgons:
                plotgon.set_alpha(0.3)
            self.polygons.clear()
            self.canvas.draw()

    def save_as_geojson(self, e=None):
        """Save current polygon to user-specified netCDF file"""
        if len(self.polygons) == 0:
            messagebox.showwarning("Save Error", "No closed polygon to save.")
            return

        now = datetime.now()
        file_path = filedialog.asksaveasfilename(
            defaultextension=".geojson",
            filetypes=[("Geojson files", "*.geojson")],
            title="Save mask as...",
            initialfile = f"mask-{now.year}{now.month}{now.day}_{now.hour}{now.minute}.geojson"
        )
        if file_path:
            polygon_to_geojson(self.polygons, file_path)
        else:
            return

        if self.clearTicked.get() == 1: # Clear if "clear on save" checkbox is ticked
            self.polygons.clear()
            self.plotgons.clear()
            self.clear_view(e="")
        else:                           # Reduce polygon opacity otherwise
            for plotgon in self.plotgons:
                plotgon.set_alpha(0.3)
            self.polygons.clear()
            self.canvas.draw()

    # -------------------
    # User interaction
    # -------------------
    def on_click(self, event):
        """Handles left (add vertex), right (close polygon), and middle (start panning) clicks"""
        if event.inaxes is None:
            return

        if (event.button == 1) and (not(self.is_poly_closed)): # Left click: add vertex
                if len(self.pointsX) > 0:
                    # Draw line from last point to current one
                    self.ax.plot([self.pointsX[-1], event.xdata], [self.pointsY[-1], event.ydata], color=self.polygon_color)

                self.pointsX.append(event.xdata)
                self.pointsY.append(event.ydata)
                self.ax.plot(event.xdata, event.ydata, color=self.polygon_color, marker="+")
                self.canvas.draw_idle()

        elif event.button == 3: # Right click
            if event.key == "control": # Ctrl click to remove polygon
                self.clean_map()
                for i in reversed(range(len(self.plotgons))):
                    poly = self.plotgons[i]
                    if poly.contains_point([event.x, event.y]):
                        if i < len(self.polygons):
                            self.polygons.pop(i)
                        self.plotgons.pop(i)
                        poly.remove()
                        self.canvas.draw()
                        break

            elif not (self.is_poly_closed) and (len(self.pointsX) > 2): # Normal click to close polygon
                # Fill the shape on the map
                maskedZone = list(zip(self.pointsX, self.pointsY))
                plotgon = Polygon(
                    maskedZone, closed=True, facecolor=self.polygon_color, edgecolor="black", alpha=0.7
                )

                # Convert points to lon/lat
                polygon = shp.Polygon([
                    ccrs.PlateCarree().transform_point(x, y, self.projection)
                    for x, y in zip(self.pointsX, self.pointsY)
                    ])

                colors = ["#ff595e", "#ffa165", "#ffca3a", "#83b731", "#1982C4", "#6A4C93"]
                self.polygon_color = colors[(colors.index(self.polygon_color) + 1) % 6]
                self.pointsX = []
                self.pointsY = []

                self.plotgons.append(plotgon)
                self.polygons.append(polygon)

                self.clean_map()

        elif event.button == 2: # Middle click: start panning
            self.panning = True
            self.pan_start = (event.x, event.y)


    def on_button_release(self, event):
        """Ends panning mode when middle button released"""
        if event.button == 2:
            self.panning = False


    def on_motion(self, event):
        """Updates coordinates in window title; handles real-time panning"""
        if event.inaxes:
            lon, lat = ccrs.PlateCarree().transform_point(
                event.xdata, event.ydata, src_crs=self.projection
            )
            root.title(f"Mirabelle - Lat: {lat:.2f}, Lon: {lon:.2f}")

        if self.panning and event.inaxes:
            dx = event.x - self.pan_start[0]
            dy = event.y - self.pan_start[1]

            x0, x1 = self.ax.get_xlim()
            y0, y1 = self.ax.get_ylim()
            width = x1 - x0
            height = y1 - y0

            self.ax.set_xlim(
                x0 - dx * width / self.canvas.get_width_height()[0],
                x1 - dx * width / self.canvas.get_width_height()[0],
            )
            self.ax.set_ylim(
                y0 - dy * height / self.canvas.get_width_height()[1],
                y1 - dy * height / self.canvas.get_width_height()[1],
            )

            self.pan_start = (event.x, event.y)
            self.canvas.draw()


    def on_scroll(self, event):
        """Zoom in/out around cursor"""
        if event.inaxes is None:
            return

        base_scale = 1.2
        scale_factor = base_scale if event.button == "down" else 1 / base_scale

        xlim = self.ax.get_xlim()
        ylim = self.ax.get_ylim()

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

        self.ax.set_xlim(new_xlim)
        self.ax.set_ylim(new_ylim)
        self.canvas.draw()


    # -------------------
    # Polygon manipulation
    # -------------------
    def clean_map(self):
        xlim = self.ax.get_xlim()
        ylim = self.ax.get_ylim()

        self.ax.clear()
        self.ax.set_xlim(xlim)
        self.ax.set_ylim(ylim)
        self.ax.coastlines()
        self.ax.add_feature(cfeature.BORDERS, linestyle=":")
        self.ax.gridlines(draw_labels=True)

        # Redraw completed (closed) polygons
        for i, poly in enumerate(self.plotgons):
            self.ax.add_patch(poly)

        # Redraw active (unfinished) polygon
        for i in range(len(self.pointsX)):
            self.ax.plot(self.pointsX[i], self.pointsY[i], color=self.polygon_color, marker="+")
            if i > 0:
                self.ax.plot([self.pointsX[i - 1], self.pointsX[i]], [self.pointsY[i - 1], self.pointsY[i]], color=self.polygon_color)
        self.canvas.draw()


    def undo_last_point(self, e=None):
        """Removes last added vertex from the active (unfinished) polygon only."""
        if len(self.pointsX) == 0:
            return

        # Remove the last point
        self.pointsX.pop()
        self.pointsY.pop()

        # Redraw map
        self.clean_map()


    def clear_view(self, e=None):
        """Resets map and clears all user input"""
        if e is None:
            answer = messagebox.askyesno(title='Confirmation',
                            message='Are you sure that you want to clear?')
        else:
            answer = True

        if answer:
            self.plotgons.clear()
            self.polygons.clear()
            self.pointsX.clear()
            self.pointsY.clear()
            self.is_poly_closed = False
            self.polygon_color = "#ff595e"

            self.ax.clear()
            self.ax.set_extent(self.map_extent, crs=ccrs.PlateCarree())
            self.ax.coastlines()
            self.ax.add_feature(cfeature.BORDERS, linestyle=":")
            self.ax.gridlines(draw_labels=True)
            self.canvas.draw()


    def reset_view(self, e=None):
        self.ax.set_extent(self.map_extent, crs=ccrs.PlateCarree())
        self.canvas.draw()


    def display_area(self):
        if self.areaTicked.get() == 0:
            for txt in self.ax.texts:
                txt.remove()
            self.canvas.draw()
            return

        for i, poly in enumerate(self.plotgons):
            transformer = pyproj.Transformer.from_crs("EPSG:4326", "EPSG:6933", always_xy=True).transform
            polygon_projected = transform(transformer, self.polygons[i])
            area = polygon_projected.area

            verts = poly.get_xy()
            x_center = verts[:, 0].mean()
            y_center = verts[:, 1].mean()

            self.ax.text(
                x_center, y_center,
                f"{area:.2e} m²",
                color='white', fontsize=7, weight='bold',
                ha='center', va='center',
                bbox=dict(facecolor='black', alpha=0.6, boxstyle='round,pad=0.2')
            )
        self.canvas.draw()


    # -------------------
    # Helper
    # -------------------
    def show_help(self):
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
                "- Left Click: Place a point on the map. If two or more points are placed, lines will automatically connect them.\n"
                "- Right Click: Close the polygon by connecting the last point to the first.\n"
                "- Middle Click: Click and drag to move the map around.\n"
                "- Scroll Wheel: Zoom in and out of the map.\n"
                "- Ctrl + Right Click: Remove the clicked polygon."
            )

        def show_menu():
            update_text(
                "Menu:\n"
                "- File > Load NetCDF: Load the reference netCDF file to mask a region on.\n"
                "- File > Save As > NetCDF: Save the selection to a NetCDF file.\n"
                "- File > Save As > GeoJSON: Export the selection as a GeoJSON file.\n"
                "- File > Help: Show this help window.\n"
                "- File > Quit: Exit the application.\n\n"
                "- Clear View: Remove all points, lines, and polygons, and reset the view.\n"
                "- Reset View: Reset the map to its original zoom level and position.\n"
                "- Undo: Remove the last point and its connecting line.\n"
                "- Clear on save (checkbox): If checked, saving clears the map. If unchecked, reduces the opacity of drawn polygons.\n"
                "- Display area (checkbox): If checked, displays the area of currently drawn polygons."
            )

        def show_shortcuts():
            update_text(
                "Keyboard Shortcuts:\n"
                "- Ctrl+S: Save as NetCDF\n"
                "- Ctrl+Shift+S: Save as GeoJSON\n"
                "- Ctrl+Z: Undo last point\n"
                "- Ctrl+R: Reset view\n"
                "- Ctrl+Backspace: Clear view\n"
                "- Ctrl+Q: Quit application"
            )


        # Buttons to switch help categories
        tk.Button(button_frame, text="Mouse Controls", command=show_mouse).pack(side=tk.LEFT, padx=5)
        tk.Button(button_frame, text="Menu", command=show_menu).pack(side=tk.LEFT, padx=5)
        tk.Button(button_frame, text="Shortcuts", command=show_shortcuts).pack(side=tk.LEFT, padx=5)

        # Show mouse help by default
        show_mouse()

# =================================================
# Mainloop

root = tk.Tk()
app = MirabelleApp(root)
root.mainloop()
