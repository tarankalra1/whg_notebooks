"""
@author: justinrogers

These are a series of tools for the coastal analysis, can be added to as needed.
"""
import numpy as np
import pandas as pd
import xarray as xr
from scipy.interpolate import griddata
from scipy.optimize import minimize
from scipy.stats import genextreme


def nearest_point_cartesian(x_grid, y_grid, point):
    """
    Find the nearest point to a given point in a grid defined by 2D X and Y points.

    Args:
        x_grid (2D numpy array): X-coordinates of the grid points
        y_grid (2D numpy array): Y-coordinates of the grid points
        point (tuple): (x, y) coordinates of the target point

    Returns:
        tuple: (x, y) coordinates of the nearest point in the grid
    """
    dist = np.sqrt((x_grid - point[0]) ** 2 + (y_grid - point[1]) ** 2)
    idx = np.unravel_index(np.argmin(dist), dist.shape)
    return idx

def get_nearest_point(LON: np.array, LAT: np.array, lon_pt: float, lat_pt: float):
    """
    This function takes input lat/lon arrays and point lat_pt/lon_pt
    and returns index of nearest point index and distance.
    """
    # function assumes -180 < lon < +180 deg
    LON2, lon_pt2 = normalize_lon(LON, lon_pt)
    dists = np.sqrt((lon_pt2 - LON2) ** 2 + (lat_pt - LAT) ** 2)
    indx = np.argmin(dists)
    dist = dists[indx]
    return indx, dist

def get_dist(LON: np.array, LAT: np.array, lon_pt: float, lat_pt: float) -> float:
    """
    This function takes input lat/lon points and point lat_pt/lon_pt
    and returns distance.
    """
    # function assumes -180 < lon < +180 deg
    LON2, lon_pt2 = normalize_lon_scalar(LON, lon_pt)
    dist = np.sqrt((lon_pt2 - LON2) ** 2 + (lat_pt - LAT) ** 2)
    return dist


def get_dists(LON: np.array, LAT: np.array, lon_pt: float, lat_pt: float) -> float:
    """
    This function takes input lat/lon vectors and point lat_pt/lon_pt
    and returns distance.
    """
    # function assumes -180 < lon < +180 deg
    LON2, lon_pt2 = normalize_lon(LON, lon_pt)
    dist = np.sqrt((lon_pt2 - LON2) ** 2 + (lat_pt - LAT) ** 2)
    return dist


def get_nearest_point_dist(
    LON: np.array,
    LAT: np.array,
    lon_pt: float,
    lat_pt: float,
    min_obs_stations: int,
) -> float:
    """
    This function takes input lat/lon arrays and point lat_pt/lon_pt, and number of
    stations. It radius of all stations to the pt, and radius of obs_station indx.
    """
    # function assumes -180 < lon < +180 deg
    LON2, lon_pt2 = normalize_lon(LON, lon_pt)
    dists = np.sqrt((lon_pt2 - LON2) ** 2 + (lat_pt - LAT) ** 2)
    max_radius = np.sort(dists)[min_obs_stations - 1] + 1e-1
    return max_radius, dists


def normalize_lon(LON: np.array, lon_pt: float) -> np.array:
    LON2 = LON[:] + 0.0
    lon_pt2 = lon_pt + 0.0
    if lon_pt2 > 180:
        lon_pt2 = lon_pt2 - 360.0
    if np.max(LON2) > 180.0:
        LON2[LON2 > 180.0] = LON2[LON2 > 180.0] - 360.0
    # if lon point is close to -180/+180 transition point
    # flip to 0 < lon < 360
    if np.abs(lon_pt2) > 150.0:
        if lon_pt2 < 0:
            lon_pt2 = lon_pt2 + 360.0
        LON2[LON2 < 0] = LON2[LON2 < 0] + 360.0
    return LON2, lon_pt2


def normalize_lon_scalar(LON: np.array, lon_pt: float) -> np.array:
    LON2 = LON + 0.0
    lon_pt2 = lon_pt + 0.0
    if lon_pt2 > 180:
        lon_pt2 = lon_pt2 - 360.0
    if LON2 > 180.0:
        LON2 = LON2 - 360.0
    # if lon point is close to -180/+180 transition point
    # flip to 0 < lon < 360
    if np.abs(lon_pt2) > 150.0:
        if lon_pt2 < 0:
            lon_pt2 = lon_pt2 + 360.0
        if LON2 < 0:
            LON2 = LON2 + 360.0
    return LON2, lon_pt2


def gaussian(x: np.array, amplitude: float, mean: float, stddev: float) -> np.array:
    """
    This function takes array x, with defined properties mean, standard deviation,
    amplitude, and returns array y.
    """
    return amplitude * np.exp(-(((x - mean) / 4 / stddev) ** 2))


def ocean_basin(ocean: pd.DataFrame, lon: np.array, lat: np.array) -> np.array:
    """
    This function takes a pandas dataframe of ocean classification points, and
    interplates to the provided lat/lon points
    """
    indx = np.isfinite(ocean["id"].values[:])
    lon_ocean = ocean["lon"].values[indx]
    lat_ocean = ocean["lat"].values[indx]
    ocean_ID = ocean["id"].values[indx]
    ocean_NAME = ocean["Ocean"].values[indx]
    # get unique id and names
    ocean_id_list = []
    ocean_name_list = []
    for i, ii in enumerate(np.unique(ocean_ID)):
        indx = np.argwhere(ocean_ID == ii)
        ocean_id_list.append(ocean_ID[indx[0, 0]])
        ocean_name_list.append(ocean_NAME[indx[0, 0]])
    ocean_id_list = np.array(ocean_id_list)
    ocean_name_list = np.array(ocean_name_list)
    print(f"interplate ocean classification")
    ocean_id = griddata((lon_ocean[:], lat_ocean[:]), ocean_ID[:], (lon, lat), method="nearest")
    # assign names to interpolated
    ocean_name = []
    for i, ii in enumerate(ocean_id):
        ocean_name.append(ocean_name_list[ii == ocean_id_list][0])
    ocean_name = np.array(ocean_name)
    return ocean_id, ocean_name


def interp_tc_rate(TC: xr, lon: np.array, lat: np.array) -> np.array:
    tc_rate = TC.tc_rate.values[:]
    lon_tc = TC.lon.values[:]
    lat_tc = TC.lat.values[:]
    lon_tc[lon_tc > 180.0] = lon_tc[lon_tc > 180.0] - 360.0
    [LON_tc, LAT_tc] = np.meshgrid(lon_tc, lat_tc)
    LON_tc = np.reshape(LON_tc, (-1, 1))
    LAT_tc = np.reshape(LAT_tc, (-1, 1))
    tc_rate = np.reshape(tc_rate, (-1, 1))
    print(f"interplate TC probability grid to points")
    tc_interp = griddata(
        (LON_tc[:, 0], LAT_tc[:, 0]), tc_rate[:, 0], (lon, lat), method="nearest"
    )

    return tc_interp


def interp_mdt(mdt: xr, lon: np.array, lat: np.array) -> np.array:
    MDT = mdt.mdt[0, :, :].values
    lon_mdt = mdt.longitude.values[:]
    lat_mdt = mdt.latitude.values[:]
    lon_mdt[lon_mdt > 180.0] = lon_mdt[lon_mdt > 180.0] - 360.0
    [LON_mdt, LAT_mdt] = np.meshgrid(lon_mdt, lat_mdt)
    LON_mdt = np.reshape(LON_mdt, (-1, 1))
    LAT_mdt = np.reshape(LAT_mdt, (-1, 1))
    MDT = np.reshape(MDT, (-1, 1))
    # only keep water values
    indx = ~np.isnan(MDT)
    LON_mdt = LON_mdt[indx]
    LAT_mdt = LAT_mdt[indx]
    MDT = MDT[indx]
    print(f"interplate MDT grid to points")
    mdt_pt = griddata((LON_mdt, LAT_mdt), MDT, (lon, lat), method="nearest")

    return mdt_pt


def interp_vdatum(vdatum: pd, lon: np.array, lat: np.array) -> np.array:
    lon_vd = vdatum["lon"].values[:]
    lat_vd = vdatum["lat"].values[:]
    offset = vdatum["offset"].values[:]
    lon_vd[lon_vd > 180.0] = lon_vd[lon_vd > 180.0] - 360.0
    print(f"interplate vdatum to points")
    vdatum_pt = griddata((lon_vd, lat_vd), offset, (lon, lat), method="linear")
    return vdatum_pt


def gumbel_cor(x, param_obs, param_model):
    # standard Gumbel GUM(mu,sigma)
    # zcor = zobs = (x - mu)/beta
    zmodel = (x - param_model[0]) / param_model[1]
    x_corr = zmodel * param_obs[1] + param_obs[0]
    return x_corr


def gev_cor(x, param_obs, param_model):
    # x = data
    # param_obs = [mu, sigma, xi], desired parameters
    # param_model = [mu, sigma, xi], parameters of x
    # Let x be a GEV(mu,sigma,xi) then
    # y = log(1 + xi*((x-mu)/sigma)))/xi
    # is a "standardized" Gumbel, i.e. a Gumbel with location 0, scale 1, and shape 0.
    # You can go the other way as well. If y is a standardized Gumbel, then
    # x = mu + sigma*((exp(xi*y)-1)/xi) is a GEV(mu,sigma,xi)
    # if xi =0, this is the gumbel distribution
    # GEV is only defined where 1 + xi*((x-mu)/sigma > 0
    if param_model[2] == 0 or param_obs[2] == 0:
        x_corr = gumbel_cor(x, param_obs, param_model)
    else:
        zmodel = (
            np.log(1 + param_model[2] * (x - param_model[0]) / param_model[1]) / param_model[2]
        )
        x_corr = (
            param_obs[0] + param_obs[1] * (np.exp(param_obs[2] * zmodel) - 1) / param_obs[2]
        )
        # check GEV limits and use gumbel for that case
        x_corr = np.where(
            1 + param_model[2] * (x - param_model[0]) / param_model[1] <= 0,
            gumbel_cor(x, param_obs, param_model),
            x_corr,
        )
    return x_corr


def fit_gev_old(wl, gev_shape_lims):
    (shape, loc, scale) = genextreme.fit(wl)
    shape = -shape  # keep normal definition
    if (shape < gev_shape_lims[0]) or (shape > gev_shape_lims[1]):
        # logger.info("shape limits exceeded, recompute...")
        shape = np.maximum(shape, gev_shape_lims[0])
        shape = np.minimum(shape, gev_shape_lims[1])
        shape, loc, scale = genextreme.fit(wl, f0=-shape)
        shape = -shape  # keep normal definition
    return shape, loc, scale


def fit_gev_set_shape_old(wl, shape):
    (shape, loc, scale) = genextreme.fit(wl, f0=-shape)
    shape = -shape  # keep normal definition
    return shape, loc, scale


def fit_gev(wl, gev_shape_lims):
    args = genextreme._fitstart(wl)
    x0, func, restore, args = genextreme._reduce_func(args, {})
    res = minimize(func, x0, args=(wl,), method="BFGS")
    # bnds = ((gev_shape_lims[0], gev_shape_lims[1]), (None, None),(None, None))
    # res = minimize(func,x0,args=(wl,),method='BFGS',bounds=bnds)

    # res.x would give you the parameter estimates
    # res.hess_inv would give you the Hessian.
    # res.x = [shape, loc, scale]
    shape = -res.x[0]
    loc = res.x[1]
    scale = res.x[2]
    if (shape < gev_shape_lims[0]) or (shape > gev_shape_lims[1]):
        # logger.info("shape limits exceeded, recompute...")
        shape = np.maximum(shape, gev_shape_lims[0])
        shape = np.minimum(shape, gev_shape_lims[1])
        # shape, loc, scale = genextreme.fit(wl, f0=-shape)
        x0, func, restore, args = genextreme._reduce_func(args, {"fc": -shape})
        res = minimize(func, x0, args=(wl,), method="BFGS")
        loc = res.x[0]
        scale = res.x[1]
    return shape, loc, scale


def fit_gev_set_shape(wl, shape):
    # (shape,loc,scale) = genextreme.fit(wl, f0=-shape)
    # shape = -shape # keep normal definition
    args = genextreme._fitstart(wl)
    x0, func, restore, args = genextreme._reduce_func(args, {"fc": -shape})
    res = minimize(func, x0, args=(wl,), method="BFGS")
    loc = res.x[0]
    scale = res.x[1]
    return shape, loc, scale


def fit_gev_hessian(wl, gev_shape_lims):
    args = genextreme._fitstart(wl)
    x0, func, restore, args = genextreme._reduce_func(args, {})
    res = minimize(func, x0, args=(wl,), method="BFGS")

    # res.x would give you the parameter estimates
    # res.hess_inv would give you the Hessian.
    # res.x = [shape, loc, scale]
    # this is a full 3x3 matrix
    hess_inv = np.reshape(
        res.hess_inv,
        (
            np.product(
                res.hess_inv.shape,
            )
        ),
    )
    shape = -res.x[0]
    loc = res.x[1]
    scale = res.x[2]
    if (shape < gev_shape_lims[0]) or (shape > gev_shape_lims[1]):
        print("shape limits exceeded, recompute...")
        shape = np.maximum(shape, gev_shape_lims[0])
        shape = np.minimum(shape, gev_shape_lims[1])
        x0, func, restore, args = genextreme._reduce_func(args, {"fc": -shape})
        res = minimize(func, x0, args=(wl,), method="BFGS")
        # this is a 2x2 matrix, add nans to places 5 to 9
        hess_inv = hess_inv + np.nan
        hess_inv[0:4] = np.reshape(
            res.hess_inv,
            (
                np.product(
                    res.hess_inv.shape,
                )
            ),
        )
        loc = res.x[0]
        scale = res.x[1]
        # bnds = ((-gev_shape_lims[1], -gev_shape_lims[0]), (None, None),(None, None))
        # res = minimize(func,x0,args=(wl,),method='L-BFGS-B',bounds=bnds)
        # shape = -res.x[0]
    return shape, loc, scale, hess_inv
