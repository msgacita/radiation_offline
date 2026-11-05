#!/usr/bin/env python3
"""
plot_rad_offline_output.py
---------------------------
Plota as variaveis do output.nc gerado pelo rad_offline_driver
(RRTMG/ecRad forcado por ERA5), sem depender do ncview -- que nao
lida bem com o _FillValue em variaveis double (trava com
"Numeric conversion not representable" em alguns builds/versoes,
mesmo com o atributo declarado corretamente no arquivo).

O xarray, ao abrir o arquivo, ja converte automaticamente qualquer
valor igual a _FillValue em NaN (mask_and_scale=True e o padrao),
entao todo o codigo abaixo trata "sem dado" naturalmente via NaN,
sem nenhum tratamento manual do sentinela.

Uso basico:
    python3 plot_rad_offline_output.py output.nc

Uso com opcoes:
    python3 plot_rad_offline_output.py output.nc \
        --outdir ./plots --time 0 --lat -22.9 --lon -43.2

Gera, dentro de --outdir (default: ./plots_<nome_do_arquivo>):
    - perfis_verticais.png   : perfil medio (media espacial) de cada
                                variavel com nivel vertical (t, q, z,
                                clwc, ciwc, crwc, cswc, cc, o3,
                                rthratensw, rthratenlw)
    - fluxos_verticais.png   : perfil medio de sw_up/sw_dn/lw_up/lw_dn
                                por half_level
    - mapas_superficie.png   : mapas horizontais de sp e skt
    - mapas_niveis.png       : mapas horizontais de cada variavel de
                                perfil em um nivel proximo da
                                superficie (900 hPa por padrao)
    - perfil_coluna_<lat>_<lon>.png : perfil vertical completo (T, q,
                                nuvem, fluxos) em UMA coluna especifica,
                                se --lat/--lon forem passados -- util
                                para inspecionar a coluna que estava
                                dando problema

Requisitos: xarray, netCDF4 (ou h5netcdf), matplotlib, numpy.
    pip install --user xarray netCDF4 matplotlib numpy
"""

import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")  # roda sem display (cluster/ssh sem X11)
import matplotlib.pyplot as plt
import xarray as xr


PROFILE_VARS = ["t", "q", "z", "clwc", "ciwc", "crwc", "cswc", "cc", "o3",
                 "rthratensw", "rthratenlw"]
FLUX_VARS = ["sw_up", "sw_dn", "lw_up", "lw_dn"]
SURFACE_VARS = ["sp", "skt"]


def open_dataset(path):
    """Abre o arquivo com mask_and_scale ligado -- _FillValue vira NaN
    automaticamente, nao precisa de nenhum tratamento manual aqui."""
    ds = xr.open_dataset(path, mask_and_scale=True)
    return ds


def _present(ds, names):
    return [v for v in names if v in ds.variables]


def plot_vertical_profiles(ds, outdir, itime):
    varnames = _present(ds, PROFILE_VARS)
    if not varnames:
        return
    ncols = 3
    nrows = int(np.ceil(len(varnames) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3.2 * nrows),
                              sharey=True)
    axes = np.atleast_1d(axes).ravel()

    plev = ds["pressure_level"].values

    for ax, name in zip(axes, varnames):
        da = ds[name].isel(valid_time=itime)
        # media espacial ignorando NaN (fill/below-ground)
        profile = da.mean(dim=("latitude", "longitude"), skipna=True).values
        n_valid_frac = 1.0 - np.isnan(da.values).mean(axis=(-1, -2))

        ax.plot(profile, plev, marker="o", markersize=2, linewidth=1)
        ax.invert_yaxis()  # pressao decrescendo pra cima (TOA no topo)
        ax.set_title(f"{name}  [{da.attrs.get('units','')}]", fontsize=9)
        ax.set_xlabel("media espacial")
        ax.grid(alpha=0.3)

        # onde a fracao de colunas validas cai (topografia/missing),
        # marca em cinza pra deixar visivel onde ha menos dado
        low_cov = n_valid_frac < 0.5
        if low_cov.any():
            ax.axhspan(plev[low_cov].min(), plev[low_cov].max(),
                       color="grey", alpha=0.15,
                       label="cobertura < 50% das colunas")
            ax.legend(fontsize=6, loc="lower right")

    for ax in axes[len(varnames):]:
        ax.axis("off")

    axes[0].set_ylabel("pressao (hPa)")
    fig.suptitle(f"Perfis verticais (media espacial) -- valid_time index {itime}",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = os.path.join(outdir, "perfis_verticais.png")
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"  -> {out}")


def plot_flux_profiles(ds, outdir, itime):
    varnames = _present(ds, FLUX_VARS)
    if not varnames:
        return
    fig, ax = plt.subplots(figsize=(5.5, 5))
    hlev = np.arange(1, ds["half_level"].size + 1)

    for name in varnames:
        da = ds[name].isel(valid_time=itime)
        profile = da.mean(dim=("latitude", "longitude"), skipna=True).values
        ax.plot(profile, hlev, marker="o", markersize=2, linewidth=1, label=name)

    ax.invert_yaxis()
    ax.set_xlabel("W m-2 (media espacial)")
    ax.set_ylabel("half_level (1 = TOA)")
    ax.set_title(f"Fluxos radiativos -- valid_time index {itime}")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    out = os.path.join(outdir, "fluxos_verticais.png")
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"  -> {out}")


def _pcolormesh_map(ax, ds, da, title, cmap="viridis"):
    lon = ds["longitude"].values
    lat = ds["latitude"].values
    pcm = ax.pcolormesh(lon, lat, da.values, shading="auto", cmap=cmap)
    ax.set_title(title, fontsize=9)
    ax.set_xlabel("lon")
    ax.set_ylabel("lat")
    plt.colorbar(pcm, ax=ax, shrink=0.85)


def plot_surface_maps(ds, outdir, itime):
    varnames = _present(ds, SURFACE_VARS)
    if not varnames:
        return
    fig, axes = plt.subplots(1, len(varnames), figsize=(6 * len(varnames), 4.5))
    axes = np.atleast_1d(axes)
    for ax, name in zip(axes, varnames):
        da = ds[name].isel(valid_time=itime)
        _pcolormesh_map(ax, ds, da, f"{name} [{da.attrs.get('units','')}]")
    fig.suptitle(f"Variaveis de superficie -- valid_time index {itime}")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    out = os.path.join(outdir, "mapas_superficie.png")
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"  -> {out}")


def plot_level_maps(ds, outdir, itime, target_hpa=900.0):
    varnames = _present(ds, PROFILE_VARS)
    if not varnames or "pressure_level" not in ds.coords:
        return
    plev = ds["pressure_level"].values
    ilev = int(np.argmin(np.abs(plev - target_hpa)))
    actual_hpa = plev[ilev]

    ncols = 3
    nrows = int(np.ceil(len(varnames) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 3.8 * nrows))
    axes = np.atleast_1d(axes).ravel()

    for ax, name in zip(axes, varnames):
        da = ds[name].isel(valid_time=itime, pressure_level=ilev)
        _pcolormesh_map(ax, ds, da, f"{name} @ {actual_hpa:.0f} hPa")

    for ax in axes[len(varnames):]:
        ax.axis("off")

    fig.suptitle(f"Mapas no nivel mais proximo de {target_hpa:.0f} hPa "
                 f"(real: {actual_hpa:.0f} hPa) -- valid_time index {itime}",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = os.path.join(outdir, "mapas_niveis.png")
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"  -> {out}")


def plot_single_column(ds, outdir, itime, lat0, lon0):
    """Perfil vertical completo em UMA coluna especifica -- util pra
    inspecionar de perto a coluna/hora que estava dando o erro
    original no ncview."""
    col = ds.sel(latitude=lat0, longitude=lon0, method="nearest")
    real_lat = float(col["latitude"].values)
    real_lon = float(col["longitude"].values)
    col = col.isel(valid_time=itime)

    fig, axes = plt.subplots(1, 3, figsize=(13, 5))

    plev = ds["pressure_level"].values
    for name in _present(ds, ["t", "q", "clwc", "ciwc", "crwc", "cswc", "cc", "o3"]):
        pass  # (t/q plotados a parte por escala; nuvem/o3 no segundo eixo)

    ax = axes[0]
    if "t" in col:
        ax.plot(col["t"].values, plev, "-o", ms=2, color="tab:red", label="t (K)")
    ax.invert_yaxis()
    ax.set_xlabel("K")
    ax.set_ylabel("pressao (hPa)")
    ax.set_title("Temperatura")
    ax.grid(alpha=0.3)

    ax = axes[1]
    for name, color in [("clwc", "tab:blue"), ("ciwc", "tab:cyan"),
                         ("crwc", "tab:green"), ("cswc", "tab:purple")]:
        if name in col:
            ax.plot(col[name].values, plev, "-o", ms=2, color=color, label=name)
    ax.invert_yaxis()
    ax.set_xlabel("kg kg-1")
    ax.set_title("Condensado (nuvem/chuva/neve)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

    ax = axes[2]
    hlev = np.arange(1, ds["half_level"].size + 1)
    for name, color in [("lw_up", "tab:orange"), ("lw_dn", "tab:red"),
                         ("sw_up", "tab:blue"), ("sw_dn", "tab:cyan")]:
        if name in col:
            ax.plot(col[name].values, hlev, "-o", ms=2, color=color, label=name)
    ax.invert_yaxis()
    ax.set_xlabel("W m-2")
    ax.set_title("Fluxos radiativos")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

    fig.suptitle(f"Coluna lat={real_lat:.2f} lon={real_lon:.2f} "
                 f"-- valid_time index {itime}")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    out = os.path.join(
        outdir, f"perfil_coluna_{real_lat:.2f}_{real_lon:.2f}.png")
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"  -> {out}")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("file", help="caminho para output.nc")
    p.add_argument("--outdir", default=None,
                   help="diretorio de saida (default: ./plots_<nome_do_arquivo>)")
    p.add_argument("--time", type=int, default=0,
                   help="indice do valid_time a plotar (default: 0)")
    p.add_argument("--level-hpa", type=float, default=900.0,
                   help="pressao alvo (hPa) para os mapas por nivel (default: 900)")
    p.add_argument("--lat", type=float, default=None,
                   help="latitude da coluna para o perfil detalhado")
    p.add_argument("--lon", type=float, default=None,
                   help="longitude da coluna para o perfil detalhado")
    args = p.parse_args()

    if not os.path.isfile(args.file):
        sys.exit(f"Arquivo nao encontrado: {args.file}")

    outdir = args.outdir or f"./plots_{os.path.splitext(os.path.basename(args.file))[0]}"
    os.makedirs(outdir, exist_ok=True)

    print(f"Abrindo {args.file} ...")
    ds = open_dataset(args.file)

    ntime = ds.sizes.get("valid_time", 1)
    if args.time >= ntime:
        sys.exit(f"--time {args.time} fora do intervalo (valid_time tem {ntime} passos)")

    print(f"Gerando figuras em {outdir}/ (valid_time index {args.time} de {ntime}) ...")
    plot_vertical_profiles(ds, outdir, args.time)
    plot_flux_profiles(ds, outdir, args.time)
    plot_surface_maps(ds, outdir, args.time)
    plot_level_maps(ds, outdir, args.time, target_hpa=args.level_hpa)

    if args.lat is not None and args.lon is not None:
        plot_single_column(ds, outdir, args.time, args.lat, args.lon)

    ds.close()
    print("Pronto.")


if __name__ == "__main__":
    main()
