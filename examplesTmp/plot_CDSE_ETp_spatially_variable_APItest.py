#%%
import requests
import datetime
import xarray as xr
import geopandas as gpd
import rioxarray

# -----------------------------
# 1️⃣ CDSE Authentication
# -----------------------------
def get_cdse_access_token(username, password):
    token_url = "https://identity.dataspace.copernicus.eu/auth/realms/CDSE/protocol/openid-connect/token"
    data = {
        "client_id": "cdse-public",
        "username": username,
        "password": password,
        "grant_type": "password"
    }
    response = requests.post(token_url, data=data)
    response.raise_for_status()
    return response.json()["access_token"]

# -----------------------------
# 2️⃣ Download function
# -----------------------------
def download_clms_file(token, s3_path, out_file):
    base_url = "https://dataspace.copernicus.eu"
    download_url = f"{base_url}{s3_path}"
    headers = {"Authorization": f"Bearer {token}"}

    with requests.get(download_url, headers=headers, stream=True) as r:
        r.raise_for_status()
        with open(out_file, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
    print(f"Saved to {out_file}")

# -----------------------------
# 3️⃣ Generate S3 paths for multiple dates
# -----------------------------
def generate_s3_paths(start_date, end_date):
    """Generate daily S3 paths based on CLMS naming convention."""
    paths = []
    current = start_date
    while current <= end_date:
        y = current.year
        m = f"{current.month:02d}"
        d = f"{current.day:02d}"
        filename = f"c_gls_SSM1km_{y}{m}{d}0000_CEURO_S1CSAR_V1.2.1_nc"
        s3_path = f"/eodata/CLMS/bio-geophysical/surface_soil_moisture/ssm_europe_1km_daily_v1/{y}/{m}/{d}/{filename}"
        paths.append((current, s3_path))
        current += datetime.timedelta(days=1)
    return paths

# -----------------------------
# 4️⃣ Main workflow
# -----------------------------
if __name__ == "__main__":
    # Credentials
    username = "benjamin.mary@ica.csic.es"
    password = "mxUCu~k,5'47"

    # Access token
    token = get_cdse_access_token(username, password)

    # Date range
    start_date = datetime.date(2026, 2, 15)
    end_date = datetime.date(2026, 2, 20)  # example: 6 days

    # Optional AOI for clipping
    aoi_path = "era5_data/agramon.geojson"  # set your AOI path
    aoi = gpd.read_file(aoi_path).to_crs("EPSG:4326")

    # Loop over dates
    for date, s3_path in generate_s3_paths(start_date, end_date):
        print(date)
        out_file = f"ssm_{date}.nc"
        download_clms_file(token, s3_path, out_file)

        # Clip to AOI
        ds = xr.open_dataset(out_file)
        ds.rio.write_crs("EPSG:4326", inplace=True)
        ds_clip = ds.rio.clip(aoi.geometry, aoi.crs)
        ds_clip.to_netcdf(f"ssm_{date}_clip.nc")
        print(f"Processed {date}")
