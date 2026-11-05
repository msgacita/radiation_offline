# Dados do ERA5 necessários — guia de download

Este documento lista exatamente o que baixar do Copernicus Climate
Data Store (CDS) para alimentar o `rrtmg_offline_driver`, e nada mais
— cada variável abaixo corresponde a um `id_*` lido em
`src/era5_bridge.F90`.

Pré-requisito: conta gratuita no CDS (https://cds.climate.copernicus.eu)
e o pacote `cdsapi` instalado (`pip install cdsapi`), com o arquivo
`~/.cdsapirc` configurado com sua chave.

---

## 1. Arquivo de níveis de pressão

Dataset: **`reanalysis-era5-pressure-levels`**

| Variável CDS (nome completo) | Nome curto lido pelo bridge | Uso | Obrigatória? |
|---|---|---|---|
| `temperature` | `t` | temperatura | **sim** |
| `specific_humidity` | `q` | vapor d'água | **sim** |
| `geopotential` | `z` | altura geopotencial (para p8w/t8w e dz8w) | **sim** |
| `specific_cloud_liquid_water_content` | `clwc` | nuvem líquida | não* |
| `specific_cloud_ice_water_content` | `ciwc` | gelo de nuvem | não* |
| `specific_rain_water_content` | `crwc` | chuva | não* |
| `specific_snow_water_content` | `cswc` | neve | não* |
| `fraction_of_cloud_cover` | `cc` | fração de nuvem | não* |
| `ozone_mass_mixing_ratio` | `o3` | ozônio real (ver nota) | não* |

\* Se omitidas, o `era5_bridge` usa um fallback (zero para
hidrometeoros/ozônio; fração de nuvem diagnosticada de forma
grosseira a partir do conteúdo de água) — mas a qualidade do resultado
cai bastante. **Recomendo baixar todas.**

Níveis de pressão recomendados (os 37 padrão do CDS — o código lida
com qualquer subconjunto, mas mais níveis = melhor resolução vertical
para o RRTMG):

```
1, 2, 3, 5, 7, 10, 20, 30, 50, 70, 100, 125, 150, 175, 200, 225, 250,
300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 775, 800, 825, 850,
875, 900, 925, 950, 975, 1000
```

### Requisição Python (cdsapi)

```python
import cdsapi

client = cdsapi.Client()
client.retrieve(
    "reanalysis-era5-pressure-levels",
    {
        "product_type": "reanalysis",
        "data_format": "netcdf",
        "variable": [
            "temperature",
            "specific_humidity",
            "geopotential",
            "specific_cloud_liquid_water_content",
            "specific_cloud_ice_water_content",
            "specific_rain_water_content",
            "specific_snow_water_content",
            "fraction_of_cloud_cover",
            "ozone_mass_mixing_ratio",
        ],
        "pressure_level": [
            "1","2","3","5","7","10","20","30","50","70","100","125",
            "150","175","200","225","250","300","350","400","450",
            "500","550","600","650","700","750","775","800","825",
            "850","875","900","925","950","975","1000",
        ],
        "year": "2024",
        "month": "07",
        "day": ["15"],
        "time": [f"{h:02d}:00" for h in range(24)],
        # área opcional: [norte, oeste, sul, leste] em graus.
        # Restrinja à sua região de interesse para acelerar o download
        # e reduzir o tamanho do arquivo — o driver lê uma coluna por
        # vez de qualquer forma.
        "area": [-15, -50, -30, -40],  # exemplo: sudeste do Brasil
        "grid": ["0.25", "0.25"],
    },
    "era5_pressure_level.nc",
)
```

---

## 2. Arquivo de nível único (superfície)

Dataset: **`reanalysis-era5-single-levels`**

| Variável CDS (nome completo) | Nome curto lido pelo bridge | Uso | Obrigatória? |
|---|---|---|---|
| `surface_pressure` | `sp` | pressão de superfície (fronteira inferior) | **sim** |
| `skin_temperature` | `skt` | temperatura de pele | **sim** |
| `forecast_albedo` | `fal` | albedo de superfície | não (fallback 0.20) |
| `sea_ice_cover` | `siconc` | fração de gelo marinho | não (fallback 0) |
| `snow_depth` | `sd` | profundidade de neve (equiv. em água) | não (fallback 0) |
| `land_sea_mask` | `lsm` | máscara terra/mar (invariante) | não (fallback = terra) |
| `geopotential` | `z` | altura da superfície (invariante) | não (fallback 0) |

### Requisição Python (cdsapi)

```python
import cdsapi

client = cdsapi.Client()
client.retrieve(
    "reanalysis-era5-single-levels",
    {
        "product_type": "reanalysis",
        "data_format": "netcdf",
        "variable": [
            "surface_pressure",
            "skin_temperature",
            "forecast_albedo",
            "sea_ice_cover",
            "snow_depth",
            "land_sea_mask",
            "geopotential",
        ],
        "year": "2024",
        "month": "07",
        "day": ["15"],
        "time": [f"{h:02d}:00" for h in range(24)],
        "area": [-15, -50, -30, -40],
        "grid": ["0.25", "0.25"],
    },
    "era5_single_level.nc",
)
```

**Nota sobre `land_sea_mask` e a `geopotential` de superfície:** são
campos invariantes no tempo. Pedir 24 cópias idênticas (uma por hora)
não é errado, só desperdiça espaço — se preferir, baixe-os separados,
num arquivo à parte com um único horário, e adapte o caminho no
`era5_bridge` (ou simplesmente aceite os fallbacks: terra em todo
ponto, altura de superfície = 0 — aceitável se seu domínio for baixo/
plano, ruim em terreno de montanha).

---

## 3. Formato do arquivo — atenção a uma armadilha real do CDS

O CDS mudou o formato de conversão GRIB→NetCDF em **novembro de 2024**.
O `era5_bridge.F90` foi escrito para aceitar **os dois formatos**
(antigo e novo), mas vale saber a diferença caso algo pareça estranho:

| | Antigo (`netcdf_legacy`) | Novo (padrão desde nov/2024) |
|---|---|---|
| Dimensão de tempo | `time` | `valid_time` |
| Dimensão de nível | `level` | `pressure_level` |
| Época de referência do tempo | `1900-01-01` | `1970-01-01` |

Se quiser forçar o formato antigo (mais simples, mais testado por
outras ferramentas como o WRF), adicione `"data_format": "netcdf_legacy"`
ao dicionário da requisição em vez de `"netcdf"`. O bridge funciona com
qualquer um dos dois — essa escolha é só uma questão de preferência.

---

## 4. Tamanho esperado do download

Para uma região regional (ex. ~15°×10°, grade 0.25°, 37 níveis, 9
variáveis, 24 horas): tipicamente **algumas centenas de MB por dia**
no arquivo de níveis de pressão. Para períodos longos (meses/anos),
considere baixar em blocos mensais e rodar o driver em lote, um mês de
cada vez.

---

## 5. Depois de baixar

```
./rrtmg_offline_driver era5_pressure_level.nc era5_single_level.nc \
    <ilon> <ilat> <itime_start> <itime_end> saida.csv
```

`ilon`/`ilat`/`itime_start`/`itime_end` são **índices** (1-based)
dentro do arquivo NetCDF, não coordenadas em graus nem datas — use
`ncdump -h era5_pressure_level.nc` para ver as dimensões e escolher os
índices certos, ou adapte o driver para aceitar lat/lon/data
diretamente se for usar em lote (fica como próximo passo natural).

---

## 6. Ordem dos níveis de pressão: crescente ou decrescente, os dois funcionam

O CDS pode entregar o arquivo de níveis de pressão com a coordenada
`pressure_level` em ordem **crescente** (1, 2, 3, ..., 1000 hPa — topo
primeiro) ou **decrescente** (1000, 975, ..., 1 hPa — superfície
primeiro), dependendo de como o pedido foi montado. As duas ordens são
detectadas automaticamente pelo `era5_bridge.F90` (comparando o primeiro
e o último valor de `pressure_level` na abertura do arquivo) e
normalizadas internamente para topo-primeiro antes de qualquer cálculo
-- você não precisa se preocupar com qual ordem seu arquivo específico
tem.

Quando a ordem decrescente é detectada, aparece este aviso (informativo,
não é erro):

```
era5_bridge: pressure_level was stored surface-to-TOA in this file --
reversing to TOA-to-surface order internally (see era5_bridge.F90).
```

**Isso corrige um bug real** encontrado ao rodar contra um arquivo real
do ERA5 com essa ordem: sem a detecção, o topo da atmosfera era
confundido com a superfície, produzindo uma extrapolação de pressão
inválida (maior que a própria superfície) logo nas duas primeiras
interfaces, disparando a rotina interna de reversão de níveis do
ecRad (`radiation_reverse`) de forma inconsistente com o resto do
código e, em seguida, uma falha de segmentação. Validado após a
correção: os mesmos dados, nas duas ordens, produzem resultado
numericamente idêntico (até a última casa decimal) tanto para o ecRad
quanto para o RRTMG.
