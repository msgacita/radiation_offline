# rad_offline — RRTMG e ecRad rodando standalone, escolhidos por namelist, forçados por ERA5

Driver Fortran standalone que roda física de radiação real — **RRTMG**
(`module_ra_rrtmg_sw.F`/`module_ra_rrtmg_lw.F`, extraídos de
`github.com/MPAS-Dev/MPAS-Model` sem modificação de física) **ou**
**ecRad** (extraído do tarball `ecrad-1.5.0_source` fornecido, com o
bug de compilação e o mesmo bug de `cloud%allocate` corrigidos — ver
`DADOS_ECRAD.md`) — escolhido em tempo de execução via namelist,
alimentado por dados horários do ERA5 (níveis de pressão) em vez de
uma previsão do MPAS.

## Status (atualizado após validação com dados reais)

- **ecRad: validado numericamente de ponta a ponta, inclusive em
  produção num cluster HPC real.** Com os dados reais do ecRad
  (repositório oficial, esquema ECCKD) e um NetCDF sintético de teste,
  o driver produziu taxas de aquecimento radiativo fisicamente
  plausíveis, tanto no ambiente de desenvolvimento quanto rodando de
  fato no cluster do usuário via `ftn` (Cray/HPE) — ver `DADOS_ECRAD.md`
  para os números e todo o processo de depuração. Isso incluiu
  encontrar e corrigir dois bugs reais: (1) `cloud%allocate` com
  `ntype` presente anula `q_liq`/`q_ice` — também corrigido
  retroativamente no acoplamento MONAN entregue antes
  (`radiation_driver_monan_mpas.F90`); (2) a **primeira chamada** ao
  ecRad dentro do processo retorna valores fisicamente absurdos em
  pelo menos um ambiente (Cray/`ftn`), corrigido por um *workaround*
  de "chamada de aquecimento" descartável antes do laço real — a
  causa raiz (provavelmente dentro do `libradiation.a` pré-compilado,
  não no código deste projeto) não foi identificada, ver
  `DADOS_ECRAD.md`.
- **RRTMG: validado numericamente de ponta a ponta.** A causa do que eu
  tinha descrito como "incompatibilidade de formato binário" era, na
  verdade, três bugs meus, todos no meu próprio código de driver (não no
  RRTMG extraído do MPAS-Model, que continua intocado): (1) faltava a flag
  de compilação `-fconvert=big-endian` no `Makefile` deste projeto — os
  arquivos `RRTMG_SW_DATA_DBL`/`RRTMG_LW_DATA_DBL` do WRF são gravados em
  *big-endian* por convenção, e sem essa flag a leitura não dá erro
  de I/O mas retorna valores sem sentido físico (exatamente o sintoma
  observado antes); (2) a variável `kme` estava definida como `nlev` em
  vez de `nlev+1` — convenção real do WRF em que os arrays de nível cheio
  sempre reservam um nível de memória extra, usado pelos arrays de
  interface (`p8w`/`t8w`); (3) os arrays de diagnóstico espectral
  opcionais (`swupflx`/`lwupflx` etc., nunca lidos por este driver) foram
  removidos da chamada, já que seu tamanho real depende de camadas
  sintéticas extras que o RRTMG adiciona internamente acima do topo do
  modelo, não do `kme+2` que a interface declarada sugere. Corrigidos os
  três, o RRTMG produziu taxas de resfriamento de onda longa fisicamente
  plausíveis (`-6×10⁻⁶` a `-3.8×10⁻⁵` K/s), na mesma ordem de grandeza e
  estrutura vertical do resultado já obtido com o ecRad — ver
  `DADOS_ECRAD.md` para os detalhes completos.
- **A leitura do ERA5 foi validada** contra um NetCDF sintético,
  incluindo o caso de borda de nível abaixo da superfície.

## Pré-requisitos

- `gfortran`, `libnetcdff-dev` (ou equivalente)
- Dois arquivos NetCDF do ERA5 (CDS) — ver `DADOS_ERA5.md`
- Para `radiation_scheme='rrtmg'`: `RRTMG_SW_DATA.DBL` e
  `RRTMG_LW_DATA.DBL` no diretório de execução
- Para `radiation_scheme='ecrad'`: a fonte do ecRad 1.5.0 extraída em
  `../ecrad_work/ecrad-1.5.0_source/` (um nível acima deste projeto —
  ver "Build" abaixo, o `make` compila o ecRad sozinho a partir daí) e
  os dados de óptica de gases/nuvens/aerossóis — ver `DADOS_ECRAD.md`
  (**recomendação: esquema ECCKD**, não RRTMG-IFS/RADRRTM)

## Build

```
make          # compila o ecRad (se ainda não estiver compilado -- ver
              # abaixo) e gera build/rad_offline_driver (RRTMG + ecRad)
              # e build/rrtmg_offline_driver (só RRTMG, por linha de
              # comando, sem depender do ecRad)
make test     # compila e roda um teste isolado da leitura do ERA5
              # (não precisa de nenhuma tabela de dados)
make clean
```

**O `make` compila o ecRad automaticamente**, chamando o `Makefile`
dele mesmo (`cd ecrad_work/ecrad-1.5.0_source && make PROFILE=gfortran`
por trás dos panos) — este projeto não duplica o sistema de build do
ecRad, só invoca. Antes de compilar, o alvo `ecrad-lib` aplica os três
arquivos de `ecrad_patches/` (correção do bug de compilação e do bug
de `cloud%allocate`, ver `DADOS_ECRAD.md`) por cima da fonte extraída
em `../ecrad_work/ecrad-1.5.0_source/radiation/`. **Isso é só a
compilação** — os dados de óptica (`data/`, ver `DADOS_ECRAD.md`) são
necessários em tempo de **execução**, não de build, e continuam
precisando ser copiados manualmente (não vêm no tarball-fonte que
vocês enviaram).

Se sua árvore do ecRad estiver em outro lugar, ou quiser usar um
`PROFILE` diferente de `gfortran`, ajuste `ECRAD_DIR`/`ECRAD_PROFILE`
no `Makefile` ou passe na linha de comando, ex.:
`make ECRAD_DIR=/caminho/para/ecrad-1.5.0_source ECRAD_PROFILE=intel`.

## Uso — driver unificado (recomendado)

```
cd build/
./rad_offline_driver <namelist.txt>
```

Exemplo de namelist (ver `namelist.offline.sample`):

```fortran
&offline_control
  radiation_scheme     = 'ecrad'        ! 'rrtmg' ou 'ecrad'
  path_pressure_level  = 'era5_pressure_level.nc'
  path_single_level    = 'era5_single_level.nc'
  ilon_start           = 1              ! índices 1-based dentro dos
  ilon_end             = 1              ! arquivos NetCDF de entrada
  ilat_start           = 1              ! (não graus/datas) -- uma
  ilat_end             = 1              ! coluna única é
                                         ! ilon_start=ilon_end e
                                         ! ilat_start=ilat_end; uma
                                         ! área é um intervalo maior
                                         ! que 1 em qualquer um dos dois
  itime_start          = 1
  itime_end            = 24
  path_output          = 'output.nc'
  ecrad_namelist_path  = 'config.nam'   ! só lido se radiation_scheme='ecrad'
/
```

Para compatibilidade com namelists antigos (só uma coluna), ainda são
aceitos `ilon`/`ilat` no lugar de `ilon_start=ilon_end`/`ilat_start=
ilat_end` — `ilon=3` equivale a `ilon_start=3, ilon_end=3`.

O driver processa `(ilon_end-ilon_start+1) × (ilat_end-ilat_start+1)`
colunas por horário, todas escritas no mesmo NetCDF de saída (uma
única rodada cobre a área inteira, sem precisar de laço externo). Como
cada coluna pode ter um número diferente de níveis acima da superfície
(`nlev_valid` varia com a topografia/pressão de superfície), a
dimensão `pressure_level`/`half_level` do arquivo de saída é sempre
fixada no total de níveis do **arquivo de entrada**, e cada
coluna/horário escreve só os níveis que tem de verdade, deixando o
resto com o valor de preenchimento padrão do NetCDF (`_FillValue`) —
nenhuma coluna ou horário é descartado por ter uma profundidade
diferente das demais.

**Saída: um único NetCDF (`path_output`, use extensão `.nc`)**, no mesmo
formato para os dois esquemas — dimensões `longitude`/`latitude`
(tamanho igual ao intervalo pedido — 1 para coluna única, maior para
área), `pressure_level` (níveis cheios, mesmos valores do arquivo de
entrada), `half_level` (interfaces,
`pressure_level`+1) e `valid_time` (um passo por horário processado).
Contém tanto as variáveis de entrada ecoadas (`t`, `q`, `z`, `clwc`,
`ciwc`, `crwc`, `cswc`, `cc`, `o3`, mais `sp`/`skt` do arquivo de
nível único) quanto a saída radiativa: `rthratensw`/`rthratenlw`
(K/s, por camada) e `sw_up`/`sw_dn`/`lw_up`/`lw_dn` (W/m², por
interface — índice 1 é o topo, o último índice é a superfície, e os
índices intermediários são o fluxo em cada fronteira de camada, então
topo/superfície/camadas saem todos do mesmo array). Um atributo global
`radiation_scheme` registra qual esquema gerou o arquivo. `level_k=1`
é sempre o nível mais alto do ERA5 nos dois esquemas, então dá pra
comparar RRTMG e ecRad diretamente (`xarray`/`ncdiff`/etc.) rodando o
mesmo namelist duas vezes só trocando `radiation_scheme` (e o
`path_output`, para não sobrescrever).

## Uso — driver antigo (só RRTMG, sem namelist, saída CSV)

Mantido para quem só quer RRTMG e prefere a interface por linha de
comando (não linka contra o ecRad, saída CSV simples, sem os perfis de
fluxo):

```
./rrtmg_offline_driver <pressure_level.nc> <single_level.nc> \
    <ilon> <ilat> <itime_start> <itime_end> <saida.csv>
```

## Estrutura de arquivos

| Arquivo | Origem |
|---|---|
| `src/module_ra_rrtmg_sw.F`, `module_ra_rrtmg_lw.F` | MPAS-Model oficial, intocados exceto o patch de ozônio (ver abaixo) |
| `src/module_ra_rrtmg_vinterp.F` | MPAS-Model oficial, intocado |
| `src/mpas_atmphys_rrtmg_swinit.F`, `_lwinit.F` | MPAS-Model oficial, intocados |
| `src/mpas_stubs_offline.F90` | stubs mínimos (MPI/log/unidades) para compilar RRTMG fora do MPAS |
| `src/radconst_offline.F90` | porte fiel de `radconst()` (declinação/constante solar) |
| `src/era5_bridge.F90` | leitura NetCDF do ERA5 e montagem da coluna |
| `src/era5_ecrad_bridge.F90` | ponte ERA5 → tipos nativos do ecRad (agora também extrai os perfis de fluxo `sw_up`/`sw_dn`/`lw_up`/`lw_dn`) |
| `src/output_netcdf.F90` | escreve o NetCDF de saída (entrada ecoada + tendências + perfis de fluxo), mesmo formato para os dois esquemas |
| `src/namelist_offline.F90` | leitura do `&offline_control` |
| `src/rad_offline_driver.F90` | **programa principal unificado** (RRTMG ou ecRad via namelist) |
| `src/rrtmg_offline_driver.F90` | driver antigo, só RRTMG, por linha de comando |
| `src/test_era5_bridge.F90` | teste isolado da leitura do ERA5 |
| `data/test_*.cdl` | NetCDF sintético de teste, coluna única (gerar com `ncgen`) |
| `data/test_area_*.cdl` | NetCDF sintético de teste, área 2×2 pontos (gerar com `ncgen`) |
| `ecrad_patches/radiation_interface.F90`, `module_ra_ecrad_sw.F90`, `module_ra_ecrad_lw.F90` | correções aplicadas automaticamente pelo `make` sobre a fonte do ecRad em `../ecrad_work/ecrad-1.5.0_source/radiation/` antes de compilar — ver `DADOS_ECRAD.md` |

## Modificações deliberadas (documentadas no código com `DELIBERATE MODIFICATION`)

1. **Ozônio no RRTMG**: quando `o3input≠2` (fora da climatologia CAM),
   o RRTMG do MPAS ignora silenciosamente qualquer perfil real passado
   em `o33d`. Como o ERA5 já tem ozônio real, um patch pequeno em
   `module_ra_rrtmg_sw.F`/`_lw.F` usa `o33d` diretamente nesse branch
   em vez do perfil padrão embutido.
2. **`radt=0.0`** no lugar de um valor arbitrário: `rrtmg_swrad`
   calcula o ângulo zenital centrado em `xtime+radt/2` (o truque de
   "média para frente" do WRF entre chamadas de radiação espaçadas).
   Como este driver já fornece a hora exata de cada passo do ERA5 e
   quer o ângulo instantâneo (não uma média), `radt` tem que ser zero
   — um valor não-zero deslocava silenciosamente o ângulo em `radt/2`
   minutos (bug encontrado e corrigido durante a validação).

Nenhuma outra linha de física foi alterada em nenhum dos dois esquemas.

## Correções do lado do driver (build e dimensionamento de arrays)

Estas não são modificações no RRTMG extraído — são bugs no código deste
projeto (`rad_offline_driver.F90`/`rrtmg_offline_driver.F90` e o
`Makefile`), encontrados só ao rodar contra os dados binários reais do
RRTMG (`RRTMG_SW_DATA_DBL`/`RRTMG_LW_DATA_DBL`):

1. **Faltava `-fconvert=big-endian`** no `Makefile`. Esses arquivos de
   dados são gravados em *big-endian* por convenção do WRF, e sem essa
   flag a leitura não dá erro de I/O mas retorna valores sem sentido
   físico — o sintoma que antes eu tinha descrito, incorretamente, como
   uma "incompatibilidade de formato binário" dos dados em si.
2. **`kme` definido errado** (`=nlev`, deveria ser `=nlev+1`). O RRTMG
   segue a convenção real do WRF em que os arrays de nível cheio
   (`t3d`, `qv3d` etc., declarados `kms:kme` na interface) sempre têm um
   nível de memória a mais do que o fisicamente necessário — é esse
   nível extra, compartilhado pela mesma variável `kme`, que os arrays
   de interface (`p8w`/`t8w`) usam para o valor de topo da atmosfera.
   Com `kme=nlev`, a escrita em `p8w(1,kme,1)` ficava um índice além do
   array realmente alocado.
3. **Arrays de diagnóstico espectral com dimensionamento dinâmico**
   (`swupflx`/`swdnflx`/`lwupflx`/`lwdnflx` e as variantes `*c` de céu
   limpo). O RRTMG estende a coluna com camadas sintéticas extras acima
   do topo do modelo (incrementos fixos de ~4~hPa até ~0~hPa) para
   preencher esses arrays de diagnóstico por nível espectral — o
   tamanho real necessário depende da pressão de topo de cada coluna,
   não do `kme+2` que a interface declarada sugere. Numa primeira
   versão eles tinham sido removidos por completo (já que não eram
   lidos), mas agora que o driver escreve o perfil de fluxo completo no
   NetCDF de saída (`sw_up`/`sw_dn`/`lw_up`/`lw_dn`), eles voltaram —
   com `kme` calculado *antes* de qualquer alocação, usando a mesma
   fórmula que o RRTMG usa internamente para estimar quantas camadas
   extras vai precisar, a partir da pressão de topo extrapolada da
   própria coluna do ERA5.

## Limitações conhecidas

- **Sem aerossóis** nos dois esquemas — o ERA5 não tem o equivalente
  aos aerossóis "water-/ice-friendly" do Thompson (RRTMG) nem uma
  fonte pronta para a razão de mistura que o ecRad espera.
- **Sem graupel** (`qg3d=0`) — o IFS/ERA5 não produz essa categoria.
- **Raio efetivo de nuvem fixo** (RRTMG: diagnosticado internamente,
  `has_reqc/reqi/reqs=0`; ecRad: 10/40 mícron fixos) — o ERA5 não
  fornece isso diretamente.
- **Emissividade de superfície fixa em 0.98** nos dois esquemas.
- Testado só com `reanalysis-era5-pressure-levels` (37 níveis fixos),
  não com o arquivo de nível de modelo (137 níveis, acesso via MARS).
- **ecRad**: sobreposição de nuvem sem ajuste por coluna/latitude
  (diferente do que foi feito no acoplamento MONAN); aerossóis não
  conectados (mesma limitação documentada em `radiation_driver_monan_mpas.F90`).
- Testado só contra dados sintéticos noturnos até agora (ambos os
  esquemas deram `rthratensw=0` corretamente, mas nunca se viu de fato
  um valor de onda curta diferente de zero nesta cadeia) — falta rodar
  contra um dia real do ERA5 para validar o ciclo diurno.
