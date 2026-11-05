"""Compare five station radiation cycles, following radiation_offline_spect.

Run: python analysis/teste_time.py --date 20250814
Requires numpy, pandas, xarray, h5netcdf, matplotlib, duckdb, and pyhdf.
Use --help to override the model results and reference dataset paths.

RRTMG SW uses valid time + radt/2 (30 minutes); model LW uses valid time.
ERA5 one-hour accumulations ending at HH+1 are converted from J/m2 to W/m2
and plotted at HH:30, matching CERES hourly midpoints. SONDA retains its
native UTC samples and is all-sky; surface model/reference fluxes are clear-sky.
IrradPy REST2v5 with AOD550=0 retains its instantaneous UTC sample times.
Missing references remain NaN with warnings. All-zero corrected SW outputs
are shown and flagged, not replaced; rerun RRTMG to regenerate those files.
"""

import argparse
from datetime import datetime, timedelta
from io import BytesIO
from pathlib import Path
import warnings
from zipfile import ZipFile, is_zipfile

import duckdb
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from pyhdf.SD import SD, SDC


PROJECT_DIR = Path(__file__).resolve().parents[1]
STATIONS = {
    "CPA": (-22.69, -45.006, "SP"),
    "CTS": (-25.508, -49.319, "SP"),
    "ORN": (-23.0, -49.844, "SP"),
    "PTR": (-9.069, -40.320, "PTR"),
    "CGR": (-20.447, -54.839, "CGR"),
}
EXPERIMENTS = {"rrtmg_orig": "RRTMG original", "rrtmg_corr": "RRTMG corrected"}
COLORS = {
    "RRTMG original": "#1f77b4",
    "RRTMG corrected": "#222222",
    "ERA5": "#b13a5b",
    "CERES": "#7551a8",
    "SONDA": "#e09f3e",
    "IrradPy AOD=0": "#009E73",
}
FIELDS = {
    "sw_surface": ("sw_dn_clear", "ssrdc", "adj_clr_sw_dn", 6),
    "sw_toa": ("sw_dn", "tisr", "toa_sw_insol", None),
    "lw_surface": ("lw_dn_clear", "strdc", "adj_clr_lw_dn", 6),
    "lw_toa": ("lw_up_clear", "ttrc", "adj_clr_lw_up", 1),
}


def nearest_index(coordinates, target, longitude=False):
    distance = np.asarray(coordinates, dtype=float) - target
    if longitude:
        distance = (distance + 180.0) % 360.0 - 180.0
    return int(np.argmin(np.abs(distance)))


def utc_hours(times, day):
    return (np.asarray(times, dtype="datetime64[ns]") - np.datetime64(day)) / np.timedelta64(1, "h")


def read_model(path, latitude, longitude, day):
    with xr.open_dataset(path) as dataset:
        if dataset.attrs.get("radiation_scheme") != "rrtmg":
            raise ValueError(f"{path}: expected radiation_scheme='rrtmg'")
        point = dataset.isel(
            latitude=nearest_index(dataset.latitude.values, latitude),
            longitude=nearest_index(dataset.longitude.values, longitude, longitude=True),
        )
        times = point.valid_time.values
        hours = utc_hours(times, day)
        if times.size != 24 or not np.array_equal(hours, np.arange(24)):
            raise ValueError(f"{path}: expected 24 hourly valid times on {day:%Y-%m-%d}")
        levels = point.t.notnull().sum("pressure_level")
        if (levels.values < 4).any():
            raise ValueError(f"{path}: missing or invalid atmospheric columns")
        fluxes = {}
        for field, (variable, _, _, _) in FIELDS.items():
            selected = point[variable].isel(half_level=levels if field.endswith("surface") else 0)
            values = np.asarray(selected.values, dtype=float)
            if values.shape != (24,) or not np.isfinite(values).all():
                raise ValueError(f"{path}: invalid {variable} fluxes at selected interfaces")
            fluxes[field] = (hours + (0.5 if field.startswith("sw") else 0.0), values)
        grid = (float(point.latitude.item()), float(point.longitude.item()))
    return times, grid, fluxes


def open_era5(path):
    if is_zipfile(path):
        with ZipFile(path) as archive:
            source = BytesIO(archive.read("data_stream-oper_stepType-accum.nc"))
        return xr.open_dataset(source, engine="h5netcdf")
    return xr.open_dataset(path)


def read_era5(dataset, latitude, longitude, times, day):
    longitude = (longitude + 180.0) % 360.0 - 180.0 if dataset.longitude.min().item() < 0 else longitude % 360.0
    point = dataset.sel(latitude=latitude, longitude=longitude)
    end_times = times + np.timedelta64(1, "h")
    missing = end_times[~np.isin(end_times, dataset.valid_time.values)]
    if missing.size:
        warnings.warn(
            f"ERA5 ({latitude}, {longitude}): missing accumulation end times {missing}; "
            "these hours remain NaN",
            stacklevel=2,
        )
    result = {}
    for field, (_, variable, _, _) in FIELDS.items():
        values = np.asarray(point[variable].reindex(valid_time=end_times).values, dtype=float) / 3600.0
        if field == "lw_toa":
            values = -values
        result[field] = (utc_hours(times, day) + 0.5, values)
    return result


def read_ceres(path):
    dataset = SD(str(path), SDC.READ)
    try:
        latitudes = dataset.select("latitude")[:]
        longitudes = dataset.select("longitude")[:]
        hours = np.asarray(dataset.select("gmt_hr_index")[:], dtype=float) - 0.5
        if not np.array_equal(hours, np.arange(24) + 0.5):
            raise ValueError(f"{path}: expected 24 hourly GMT intervals")
        levels = dataset.select("level")[:]
        result = {}
        for station, (latitude, longitude, _) in STATIONS.items():
            ilat = nearest_index(latitudes, latitude)
            ilon = nearest_index(longitudes, longitude, longitude=True)
            result[station] = {}
            for field, (_, _, name, level) in FIELDS.items():
                variable = dataset.select(name)
                expected_dims = (
                    ("gmt_hr_index", "latitude", "longitude") if level is None
                    else ("level", "gmt_hr_index", "latitude", "longitude")
                )
                if tuple(variable.dimensions()) != expected_dims:
                    raise ValueError(f"{path}: unexpected dimensions for {name}")
                values = np.asarray(
                    variable[:, ilat, ilon] if level is None
                    else variable[int(np.flatnonzero(levels == level).item()), :, ilat, ilon],
                    dtype=float,
                )
                attributes = variable.attributes()
                lower, upper = attributes["valid_range"]
                invalid = (
                    ~np.isfinite(values) | (values == attributes["_FillValue"])
                    | (values < lower) | (values > upper)
                )
                if invalid.any():
                    warnings.warn(f"CERES {station} {name}: {invalid.sum()} invalid hours", stacklevel=2)
                    values[invalid] = np.nan
                result[station][field] = (hours, values)
        return result
    finally:
        dataset.end()


def read_sonda(path, day):
    with duckdb.connect() as connection:
        observations = connection.execute(
            """
            SELECT acronym, CAST(timestamp AS TIMESTAMP) AS timestamp, glo_avg
            FROM read_parquet(?)
            WHERE acronym IN (SELECT unnest(?))
                AND CAST(timestamp AS TIMESTAMP) >= ?
                AND CAST(timestamp AS TIMESTAMP) < ?
            ORDER BY acronym, timestamp
            """,
            [str(path), list(STATIONS), day, day + timedelta(days=1)],
        ).fetchdf()
    result = {}
    for station in STATIONS:
        selected = observations.loc[observations.acronym.eq(station)]
        if selected.empty:
            raise ValueError(f"SONDA: no observations for {station} on {day:%Y-%m-%d}")
        values = selected.glo_avg.to_numpy(dtype=float)
        invalid = ~np.isfinite(values) | (values < 0)
        if invalid.any():
            warnings.warn(f"SONDA {station}: {invalid.sum()} invalid GHI samples", stacklevel=2)
            values[invalid] = np.nan
        result[station] = (utc_hours(selected.timestamp.to_numpy(), day), values)
    return result


def read_irradpy(directory, day):
    result = {}
    for station, (latitude, longitude, _) in STATIONS.items():
        path = directory / f"rest2_{station}_{day:%Y%m%d}.csv"
        frame = pd.read_csv(path)
        times = pd.to_datetime(frame["time_utc"], utc=True, errors="raise")
        hours = utc_hours(times.dt.tz_localize(None).to_numpy(), day)
        if (
            not frame["station"].eq(station).all()
            or not np.allclose(frame["latitude"], latitude)
            or not np.allclose(frame["longitude"], longitude)
            or not frame["aod550"].eq(0.0).all()
            or not frame["model"].eq("IrradPy REST2v5").all()
            or len(hours) != 24
            or not np.array_equal(hours, np.arange(24) + 0.5)
        ):
            raise ValueError(f"{path}: unexpected station, AOD, model, or timestamps")
        values = frame["ghi_W_m2"].to_numpy(dtype=float)
        if not np.isfinite(values).all() or (values < 0).any():
            raise ValueError(f"{path}: invalid IrradPy GHI")
        result[station] = (hours, values)
    return result


def plot_cycles(panels, radiation, output_dir, date, zoom=False):
    fields = ("sw_surface", "sw_toa") if radiation == "sw" else ("lw_surface", "lw_toa")
    titles = {
        "sw_surface": "Surface SW down (GHI)",
        "sw_toa": "TOA incoming SW",
        "lw_surface": "Surface LW down",
        "lw_toa": "TOA LW up",
    }
    fig, axes = plt.subplots(2, 5, figsize=(23, 9), sharex=True, squeeze=False, layout="constrained")
    for column, (station, series) in enumerate(panels.items()):
        for row, field in enumerate(fields):
            axis = axes[row, column]
            for source, (hours, values) in series[field].items():
                label = f"{source} (all-sky)" if source == "SONDA" else source
                axis.plot(
                    hours, values, color=COLORS[source], label=label,
                    marker=None if source == "SONDA" else "o", markersize=3,
                    linewidth=2.2 if source == "RRTMG corrected" else 1.0 if source == "SONDA" else 1.7,
                    linestyle="--" if source == "RRTMG corrected" else "-",
                    zorder=10 if source == "RRTMG corrected" else 2,
                )
            if field.startswith("sw") and np.all(series[field]["RRTMG corrected"][1] == 0):
                axis.text(
                    0.02, 0.94, "Corrected run: all-zero SW\nRerun required",
                    transform=axis.transAxes, va="top", color="#b22222", fontsize=9,
                )
            axis.set_xlim((12, 18) if zoom else (0, 24))
            axis.set_xticks(np.arange(12, 19) if zoom else np.arange(0, 25, 3))
            if zoom:
                visible = [
                    value for hours, values in series[field].values()
                    for hour, value in zip(hours, values)
                    if 12 <= hour <= 18 and np.isfinite(value)
                ]
                if not visible:
                    raise ValueError(f"{station}: no finite {field} samples between 12 and 18 UTC")
                padding = max((max(visible) - min(visible)) * 0.05, 1.0)
                axis.set_ylim(min(visible) - padding, max(visible) + padding)
            axis.set_title(f"{station}: {titles[field]}")
            axis.set_ylabel("Flux (W m$^{-2}$)")
            axis.grid(alpha=0.3)
            if row == 1:
                axis.set_xlabel("Hour (UTC)")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=len(labels))
    fig.suptitle(
        f"Diurnal radiation cycles: {date}\n"
        "Surface models, ERA5 and CERES: clear-sky; SONDA: all-sky. "
        "SW/ERA5/CERES: hourly midpoints; model LW: valid time."
    )
    suffix = "_12_18UTC" if zoom else ""
    path = output_dir / f"{radiation}_diurnal_cycle_{date}_CPA_CTS_ORN_PTR_CGR{suffix}.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)
    print(f"Saved {path}")
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--date", default="20250814", help="Day in YYYYMMDD format")
    parser.add_argument("--results-dir", type=Path)
    parser.add_argument("--era5-path", type=Path)
    parser.add_argument("--sonda-path", type=Path, default=PROJECT_DIR.parent / "dados_sonda/Solarimetrica.parquet")
    parser.add_argument("--ceres-path", type=Path)
    parser.add_argument("--irradpy-dir", type=Path, help="Directory containing zero-AOD REST2 station CSV files")
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent / "figures")
    args = parser.parse_args()
    day = datetime.strptime(args.date, "%Y%m%d")
    results = args.results_dir or PROJECT_DIR / "rrtmg_offline" / f"results_{args.date}"
    era5_path = args.era5_path or (
        PROJECT_DIR.parent / "radiation_offline_spect/rrtmg_offline"
        / f"era5_single_level_radiation_{args.date}.nc"
    )
    ceres_path = args.ceres_path or (
        PROJECT_DIR.parent / "global_data/CERES/CER_SYN1deg"
        / f"CER_SYN1deg-1Hour_Terra-Aqua-NOAA20_Edition4B_415412.{args.date}"
    )
    ceres = read_ceres(ceres_path)
    sonda = read_sonda(args.sonda_path, day)
    irradpy_dir = args.irradpy_dir or (
        PROJECT_DIR.parent / "radiation_offline_spect/analysis" / f"rest2_results_{args.date}_aod0"
    )
    irradpy = read_irradpy(irradpy_dir, day)
    panels = {}
    records = []
    with open_era5(era5_path) as era5:
        for station, (latitude, longitude, region) in STATIONS.items():
            series = {field: {} for field in FIELDS}
            reference_times = reference_grid = None
            for experiment, label in EXPERIMENTS.items():
                path = results / region / experiment / "output.nc"
                times, grid, model = read_model(path, latitude, longitude, day)
                if reference_times is not None:
                    if not np.array_equal(times, reference_times) or grid != reference_grid:
                        raise ValueError(f"{station}: experiments have different times or selected grid points")
                reference_times, reference_grid = times, grid
                for field in FIELDS:
                    series[field][label] = model[field]
                if experiment == "rrtmg_corr" and np.all(model["sw_toa"][1] == 0):
                    warnings.warn(
                        f"{station}: corrected RRTMG SW is zero all day in {path}; "
                        "rerun after removing the driver's coszr reset",
                        stacklevel=2,
                    )
            for source, reference in (
                ("ERA5", read_era5(era5, *reference_grid, reference_times, day)),
                ("CERES", ceres[station]),
            ):
                for field in FIELDS:
                    series[field][source] = reference[field]
            series["sw_surface"]["SONDA"] = sonda[station]
            series["sw_surface"]["IrradPy AOD=0"] = irradpy[station]
            panels[station] = series
            print(f"{station}: station=({latitude}, {longitude}), model/ERA5 grid={reference_grid}")
            for field, datasets in series.items():
                for source, (hours, values) in datasets.items():
                    for hour, value in zip(hours, values):
                        records.append({
                            "station": station, "date": args.date, "hour_utc": hour,
                            "field": field, "source": source, "flux_W_m2": value,
                        })
    args.output_dir.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame.from_records(records)
    csv_path = args.output_dir / f"diurnal_cycles_{args.date}_CPA_CTS_ORN_PTR_CGR.csv"
    frame.to_csv(csv_path, index=False)
    print(f"Saved {csv_path}")
    plot_cycles(panels, "sw", args.output_dir, args.date)
    plot_cycles(panels, "sw", args.output_dir, args.date, zoom=True)
    plot_cycles(panels, "lw", args.output_dir, args.date)


if __name__ == "__main__":
    main()