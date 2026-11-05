# Dados do ecRad necessários

O ecRad precisa de um conjunto de arquivos de coeficientes ópticos
(gases, nuvens, aerossóis) carregados em tempo de execução a partir do
diretório apontado por `directory_name` no `config.nam` (padrão:
`"data"`, relativo ao diretório de onde o executável é chamado). Esses
arquivos **não vêm no tarball de código-fonte** que você enviou
(`ecrad-1_5_0_source_tar.gz` — a pasta `data/` está vazia nele).

> **Nota de escopo**: este documento é sobre os dados do **ecRad**
> (`RADRRTM`/`RADSRTM`, esquema RRTMG-IFS interno do ecRad, ou os
> arquivos NetCDF do ECCKD). Isso é um caso **diferente** dos dados do
> **RRTMG standalone** extraído do MPAS-Model
> (`RRTMG_SW_DATA_DBL`/`RRTMG_LW_DATA_DBL`, do repositório do WRF) — a
> incompatibilidade que eu tinha reportado ali **foi resolvida** (era só
> uma flag de compilação faltando no meu Makefile, `-fconvert=big-endian`
> — ver `README.md`, seção "Correções do lado do driver"). A
> incompatibilidade do `RADRRTM` do ecRad descrita abaixo é diferente:
> confirmei que o build do próprio ecRad **já** usa
> `-fconvert=big-endian` (está no `Makefile_include.gfortran` dele), e
> mesmo assim a leitura falha — não é o mesmo bug, continua sem solução.

## Onde conseguir

O repositório oficial **é público** e já inclui a pasta `data/`
completa (~22 MB):

```
git clone https://github.com/ecmwf-ifs/ecrad.git
cp -r ecrad/data/* /caminho/para/seu/ecrad-1.5.0_source/data/
```

## Achado real ao validar — dois esquemas de gases, escolha recomendada

O ecRad suporta dois modelos de óptica de gases, escolhidos via
`gas_model_name` no `config.nam`:

| `gas_model_name` | Arquivo(s) | Formato | Status testado aqui |
|---|---|---|---|
| `"RRTMG-IFS"` (padrão do `practical/config.nam`) | `RADRRTM`, `RADSRTM` | binário Fortran não-formatado (registro sequencial) | **Incompatibilidade encontrada** — ver abaixo |
| `"ECCKD"` (desde ecRad 1.5) | `ecckd-*_ckd-definition.nc` | NetCDF | **Validado com sucesso, fim a fim, com resultado numérico fisicamente plausível** |

### O que descobri sobre o `RADRRTM`

Ao copiar o `data/` do repositório oficial (branch `master` atual) para
dentro desta fonte 1.5.0 e tentar usar `gas_model_name = "RRTMG-IFS"`,
a leitura do `RADRRTM` **tecnicamente não dá erro de I/O** (`iostat=0`)
mas retorna valores sem sentido físico (ex.: `rayl=7.5E-226`,
`layreffr=301989888`). Isolei isso com um programa de teste mínimo
lendo só o primeiro registro do arquivo — confirmei que é uma
incompatibilidade de **formato de marcador de registro binário**
(provável: o arquivo do `master` atual foi gravado por um compilador
diferente do gfortran, ou por uma versão do ecRad com layout de
registro levemente diferente desta fonte 1.5.0 específica). Isso não é
um bug no código extraído — é um descompasso de proveniência entre
"dado de uma versão" e "código de outra versão", exatamente o tipo de
armadilha binária que motivou preferir NetCDF (ECCKD) sempre que
possível.

**Recomendação:** use `gas_model_name = "ECCKD"` neste projeto. Além de
mais robusto a diferenças de versão (NetCDF é autodescritivo), é
também o esquema mais moderno do ecRad.

### Configuração mínima testada e validada (ECCKD)

No `config.nam`, a partir do `practical/config.nam` de exemplo, mude:

```
gas_model_name      = "ECCKD",
gas_optics_sw_override_file_name = "ecckd-1.0_sw_climate_rgb-32b_ckd-definition.nc",
gas_optics_lw_override_file_name = "ecckd-1.0_lw_climate_fsck-32b_ckd-definition.nc",
use_general_cloud_optics = true,     ! ECCKD EXIGE isso (senão erro de config)
use_general_aerosol_optics = true,   ! ECCKD EXIGE isso também
use_aerosols = false,                ! a ponte ERA5 não alimenta aerossóis (ver era5_ecrad_bridge.F90)
```

O `practical/config.nam` original vem com `use_general_cloud_optics =
false`, `use_general_aerosol_optics = false` e `use_aerosols = true` —
os três precisam ser trocados exatamente como acima, ou o `setup_radiation`
falha com um erro de configuração claro (`*** Error: ecCKD gas optics
model requires general cloud optics`), ou o `driver_radiation_lw_sw`
tenta usar `aerosol%mixing_ratio` sem estar alocado (segfault).

## Validação real obtida

Com essa configuração, `data/ecckd-1.0_*` copiados e a `namelist.ecrad.test`
de exemplo, o driver rodou de ponta a ponta com física real e produziu:

```
scheme,itime,julian,gmt,level_k,z_m,p_pa,rthratensw_K_s,rthratenlw_K_s
ecrad,1,197.0000,.0000,1, 9198.0, 30000.0, -0.0, -2.10E-05
ecrad,1,197.0000,.0000,2, 5680.0, 50000.0, -0.0, -2.48E-05
...
```

`rthratensw=0` está correto (era noite no horário/longitude do teste);
`rthratenlw` em torno de -2 a -3×10⁻⁵ K/s (≈ -1.7 a -2.9 K/dia de
resfriamento radiativo de onda longa) é fisicamente plausível para a
troposfera livre.

## Bug real encontrado e corrigido durante essa validação

Ao rodar com dados reais, um segfault revelou um bug real em
`era5_ecrad_bridge.F90` (e no arquivo irmão `radiation_driver_monan_mpas.F90`
do projeto de acoplamento MONAN, que tinha o mesmo padrão): chamar
`cloud%allocate(ncol, nlev, 2, .false.)` com o terceiro argumento
(`ntype`) presente faz o ecRad entrar num ramo interno que **anula**
os ponteiros `cloud%q_liq`/`cloud%q_ice` em vez de associá-los —
eles só ficam utilizáveis quando `ntype` é **omitido** da chamada.
Isso nunca tinha sido pego porque nenhuma das duas entregas anteriores
havia rodado contra dados reais até agora — o Fortran não detecta isso
em tempo de compilação. Os dois arquivos já foram corrigidos.

## Bug de "primeira chamada" — workaround aplicado (causa raiz não identificada)

Ao rodar em produção num cluster HPC real (compilador via `ftn`,
Cray/HPE), apareceu um segundo bug real, mais sutil: **a primeira
chamada ao ecRad dentro do processo** retornava taxas de aquecimento
de onda longa fisicamente absurdas perto da superfície (`>0.1 K/s`,
ou seja, dezenas de milhares de K/dia), enquanto **qualquer chamada
seguinte no mesmo processo** — mesmo com dados de entrada quase
idênticos (a mesma coluna, uma hora depois) — retornava valores
corretos.

**Isolamento do problema:** rodando `itime=1` e `itime=2` cada um
sozinho (como única chamada do processo), os dois davam o valor
errado — confirmando que o bug segue "é a primeira chamada do
processo", não "é o dado específico do itime=1".

**Não reproduzi isso no meu ambiente** (gfortran comum), mesmo usando
`-fcheck=all -finit-real=snan -Wuninitialized` no meu código de ponte
(`era5_ecrad_bridge.F90`) — o que aponta fortemente para o problema
estar em algum estado interno do **próprio `libradiation.a`**
(alguma variável `SAVE`/cache que só fica devidamente inicializada
depois de uma primeira chamada completa a `driver_radiation_lw_sw`),
não no código deste projeto. Como o binário do ecRad usado é
pré-compilado, minhas flags de depuração no meu próprio código não
alcançam o que acontece dentro da biblioteca.

**Workaround aplicado** (não é correção da causa raiz): `rad_offline_driver.F90`
agora faz uma chamada de "aquecimento" ao ecRad logo depois da
inicialização — roda uma vez com o primeiro horário disponível e
**descarta o resultado** — antes de entrar no laço de tempo real.
Assim, qualquer efeito de primeira-chamada é absorvido por essa
chamada descartável, e todo valor que vai para o CSV de saída já vem
de uma segunda-ou-mais chamada (a faixa comprovadamente correta).
Validado em produção: depois do workaround, `itime=1` e `itime=2`
passaram a dar exatamente a mesma faixa de valores fisicamente
plausíveis (`-2` a `-3×10⁻⁵` K/s), batendo com o que já era obtido
no ambiente de desenvolvimento.

**Se algum dia quiserem investigar a causa raiz de verdade**: o
próximo passo seria recompilar o `libradiation.a` do próprio ecRad
com `-fcheck=all -finit-real=snan -g` (editando `FFLAGS_ADD` no
`Makefile_include.gfortran` do ecRad) e rodar com `gdb`/`ftn`-debug
no mesmo cluster onde o bug aparece, procurando por uma variável
`SAVE` ou um array alocado sem inicialização explícita em algum lugar
do caminho de `driver_radiation_lw_sw` → `gas_optics`/`cloud_optics`
→ solver de onda longa. Não fiz isso aqui por ser uma investigação
grande dentro de código de terceiros (ECMWF), e o workaround já
resolve o problema prático.
