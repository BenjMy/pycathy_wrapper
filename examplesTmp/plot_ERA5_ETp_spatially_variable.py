"""
ERA5 ETp download, process, and plot (download only if not exists)
"""

import os
from pathlib import Path
import zipfile

import cdsapi
import xarray as xr
import rioxarray as rxr
import matplotlib.pyplot as plt
import geopandas as gpd

# ----------------------------
# Paths
# ----------------------------
rootPath = Path(os.getcwd())
download_path = rootPath / "era5_download.zip"
data_dir = rootPath / "era5_data"
data_dir.mkdir(exist_ok=True)

# Optional: AOI shapefile
aoi_shp = rootPath / "era5_data/agramon.geojson"  # set your AOI path

# ----------------------------
# Check if download exists
# ----------------------------
if not download_path.exists():
    # ----------------------------
    # CDS API download request
    # ----------------------------
    dataset = "reanalysis-era5-single-levels"
    request = {
        'product_type': ['reanalysis'],
        'variable': [
                    # '2m_temperature',
                     # 'total_precipitation',
                     # 'evaporation',
                     'potential_evaporation'
                     ],
        'year': ['2022'],
        'month': [f"{i:02d}" for i in range(1, 4)],  # Jan-Mar
        'day': [f"{i:02d}" for i in range(1, 32)],   # all days
        'time': [f"{i:02d}:00" for i in range(24)],  # hourly
        'data_format': 'netcdf',
        'download_format': 'zip',
        'area': [43.82, -9.37, 35.95, 3.39]  # N, W, S, E
    }

    client = cdsapi.Client()
    print("Downloading ERA5 data...")
    client.retrieve(dataset, request).download(str(download_path))
    print("Download complete.")
else:
    print(f"Download already exists: {download_path}")

# ----------------------------
# Unzip (if folder empty)
# ----------------------------
nc_files = list(data_dir.glob("*.nc"))
if len(nc_files) == 0:
    print("Extracting zip...")
    with zipfile.ZipFile(download_path, 'r') as zip_ref:
        zip_ref.extractall(data_dir)
    print("Extraction complete.")
    nc_files = list(data_dir.glob("*.nc"))
else:
    print("Data already extracted.")

if len(nc_files) == 0:
    raise FileNotFoundError("No .nc files found in extracted folder.")

# ----------------------------
# Open NetCDF with xarray/rioxarray
# ----------------------------
ds = xr.open_dataset(nc_files[0])
ds = ds.rio.write_crs("EPSG:4326")  # assign CRS
ds.rio.resolution()


# ----------------------------
# Extract potential evaporation (ETp)
# ----------------------------
etp = ds["pev"]
etp.valid_time

# Resample to daily sum
etp_daily = etp.resample(valid_time="1D").sum()

# ----------------------------
# Optional: clip to AOI
# ----------------------------
if aoi_shp.exists():
    shape = gpd.read_file(aoi_shp).to_crs("EPSG:4326")
    shape.geometry = shape.geometry.buffer(0.1)  # 0.01 deg ~ 1 km
    etp_clip = etp_daily.isel(valid_time=0).rio.clip(shape.geometry, shape.crs)
else:
    etp_clip = etp_daily.isel(valid_time=0)


print(etp_daily.rio.crs)   # raster CRS
print(shape.crs)           # AOI CRS


# ----------------------------
# Plot
# ----------------------------
fig, ax = plt.subplots(figsize=(8,6))
# etp_clip.rio.plot(ax=ax, cmap="viridis")
etp_clip.plot.imshow(ax=ax, cmap="viridis")
ax.set_title("ERA5 Potential Evaporation (ETp, mm/day)")
plt.show()

# ----------------------------
# Save daily ETp NetCDF for CATHY atmbc
# ----------------------------
output_file = rootPath / "ETp_ERA5_2022_Q1.nc"
etp_daily.to_netcdf(output_file)
print(f"Daily ETp saved to {output_file}")
etp_clip.rio.to_raster("etp_clip.tif")
