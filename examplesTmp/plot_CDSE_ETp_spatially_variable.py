#%%
import requests
import datetime
import xarray as xr
import geopandas as gpd
import rioxarray
from pathlib import Path
import os
import matplotlib.pyplot as plt


# ----------------------------
# Paths
# ----------------------------
rootPath = Path(os.getcwd())
data_dir = rootPath / "CDSE_data"
data_dir.mkdir(exist_ok=True)

aoi_shp = rootPath / "CDSE_data/agramon.geojson"  # set your AOI path

nc_files = list(data_dir.glob("*.nc"))

# ----------------------------
# Open NetCDF with xarray/rioxarray
# ----------------------------
ds = xr.open_dataset(nc_files[0])
ds = ds.rio.write_crs("EPSG:4326")  # assign CRS
ds.rio.resolution()


# ----------------------------
# Extract potential evaporation (ETp)
# ----------------------------
ssm = ds["ssm"]
ssm.time

# ----------------------------
# Optional: clip to AOI
# ----------------------------
if aoi_shp.exists():
    shape = gpd.read_file(aoi_shp).to_crs("EPSG:4326")
    shape.geometry = shape.geometry.buffer(0.1)  # 0.01 deg ~ 1 km
    ssm_clip = ssm.isel(time=0).rio.clip(shape.geometry, shape.crs)
else:
    ssm_clip = ssm.isel(time=0)


print(ssm.rio.crs)   # raster CRS
print(shape.crs)           # AOI CRS


# ----------------------------
# Plot
# ----------------------------

fig, ax = plt.subplots(figsize=(8,6))
# etp_clip.rio.plot(ax=ax, cmap="viridis")
ssm.isel(time=0).plot.imshow(ax=ax, cmap="viridis")
ax.set_title("ERA5 Potential Evaporation (ETp, mm/day)")
plt.show()


fig, ax = plt.subplots(figsize=(8,6))
# etp_clip.rio.plot(ax=ax, cmap="viridis")
ssm_clip.plot.imshow(ax=ax, cmap="viridis")
ax.set_title("ERA5 Potential Evaporation (ETp, mm/day)")
plt.show()

# # ----------------------------
# # Save daily ETp NetCDF for CATHY atmbc
# # ----------------------------
# output_file = rootPath / "ETp_ERA5_2022_Q1.nc"
# ssm_clip.to_netcdf(output_file)
# print(f"Daily ETp saved to {output_file}")
# etp_clip.rio.to_raster("etp_clip.tif")


