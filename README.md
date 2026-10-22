# SFINCS practical guide: physics, grids, inputs, and forcing

This guide is intended for coastal and estuarine numerical modelers who are new to SFINCS. It is based on the Fortran source in this checkout and on the Sphinx manual in [`docs/`](docs/). The project overview remains in [`README.rst`](README.rst); this file focuses on how the computational engine interprets grids and forcing files.

> **Version note:** this checkout identifies itself as SFINCS `v2.4.0 Galibier+` in `source/src/sfincs_lib.f90`. Where the source and older manual text disagree, the behavior below follows the source. In particular, `advection = 1` is the current source-code default.

## Contents

- [1. Model concept and physics](#1-model-concept-and-physics)
- [2. Where the model logic lives](#2-where-the-model-logic-lives)
- [3. Preparing the computational grid](#3-preparing-the-computational-grid)
- [4. Main input file](#4-main-input-file)
- [5. Boundary and forcing files](#5-boundary-and-forcing-files)
- [6. Boundary-condition generation workflow](#6-boundary-condition-generation-workflow)
- [7. Running, outputs, and quality control](#7-running-outputs-and-quality-control)
- [8. Common pitfalls](#8-common-pitfalls)
- [9. References and further documentation](#9-references-and-further-documentation)

## 1. Model concept and physics

SFINCS (**Super-Fast INundation of CoastS**) is a reduced-complexity, depth-averaged model for dynamic compound flooding. It is designed for event-scale coastal flood simulations in which tide, storm surge, waves, wind, precipitation, infiltration, and river discharge may all contribute.

### 1.1 Governing equations

The core solver integrates mass continuity and horizontal momentum on a staggered, rectilinear mesh. In simplified form, the continuity equation is

\[
\frac{\partial \zeta}{\partial t}
+ \frac{\partial q_x}{\partial x}
+ \frac{\partial q_y}{\partial y}
= P - I + \frac{Q_s}{A},
\]

where \(\zeta\) is the free-surface elevation, \(q_x\) and \(q_y\) are unit-width fluxes, \(P\) is precipitation, \(I\) is infiltration, \(Q_s\) represents point sources/sinks, and \(A\) is the cell area. In subgrid mode, SFINCS first updates cell volume and then converts volume to water level through precomputed hypsometric tables.

The momentum update follows the computationally efficient Bates et al. (2010) formulation. Conceptually, each face flux is accelerated by water-surface gradient and optional process terms, then damped by bottom friction:

\[
q^{n+1} =
\frac{q_s + \Delta t\,F}
{1 + \Delta t\,g n^2 |q|/h^{7/3}},
\]

with \(F\) containing combinations of the barotropic water-level gradient, optional advection, Manning friction, wind stress, atmospheric-pressure gradients, wave forces, Coriolis acceleration, viscosity, and structure/wavemaker corrections.

Water levels live at cell centers and fluxes at cell faces. Wetting and drying are controlled primarily by `huthresh` (default `0.05 m`). A face is not treated as dynamically wet until its effective depth exceeds that threshold. Flux limiters prevent water from leaving a dry or negative-volume cell.

### 1.2 Numerical integration and timestep

The solver is explicit and adaptive. Each face computes a local timestep from the larger of the shallow-water wave speed \(\sqrt{gh}\) and the local velocity, divided by grid spacing. The global timestep is the minimum value, multiplied by `alpha` and capped by `dtmax`:

```text
dt = alpha * min_i [ 1 / (max(sqrt(g h_i), |u_i|) / dx_i) ]
```

Useful stability controls are `alpha` (default `0.5`), `dtmax` (default `60 s`), `hmin_cfl` (default `0.1 m`), `huthresh` (default `0.05 m`), `theta` (default `1.0`), and optional `viscosity`/`nuvisc`. Set `timestep_analysis = 1` to write diagnostics identifying cells that most often limit the global timestep.


### 1.3 Reduced-physics hydrodynamics in more detail

Let \(z_b(x,y)\) be bed elevation, \(\zeta(x,y,t)\) water-surface elevation, and \(h=\max(\zeta-z_b,0)\) total depth. The depth-integrated shallow-water momentum equations can be written schematically as

\[
\frac{\partial \mathbf{q}}{\partial t}
+\nabla\cdot\left(\frac{\mathbf{q}\mathbf{q}}{h}\right)
+g h\nabla\zeta
+\frac{g n^2 |\mathbf{q}|\mathbf{q}}{h^{7/3}}
= \frac{\boldsymbol{\tau}_s}{\rho_w}
-\frac{h}{\rho_w}\nabla p_a
+\mathbf{F}_{wave}
+\mathbf{F}_{Coriolis}
+\mathbf{F}_{viscosity},
\]

and continuity as

\[
\frac{\partial h}{\partial t}+\nabla\cdot\mathbf{q}=P-I+S.
\]

SFINCS solves these equations in a deliberately reduced way:

- **Flux form:** face flux \(q\) is the primary momentum variable, not cell-centered velocity.
- **Semi-implicit friction:** the Manning denominator in the update stabilizes bottom drag without requiring an iterative nonlinear solve.
- **Optional inertia:** with `advection = 0`, the nonlinear advective term is omitted (often called SFINCS-LIE). With `advection = 1`, a first-order upwind advection approximation is active (SFINCS-SSWE behavior in the documentation). The current source default is enabled.
- **Hydrostatic pressure:** the main barotropic pressure gradient is represented by \(g h \nabla\zeta\). `nonh` adds a pressure-correction method for problems in which vertical acceleration is important; it is not part of the default flood-mode formulation.
- **Wetting/drying:** the model does not solve a complex moving shoreline at every sub-face. It uses depth thresholds, positivity constraints, and in subgrid mode tabulated partial-wetting information.
- **No 3D turbulence closure:** wind stress, bed friction, wave forcing, and vegetation/structures are parameterized as depth-integrated source terms.

This formulation is most defensible for rapidly varying coastal inundation, floodplain flow, and compound flooding where speed and stable wetting/drying are more important than resolving vertical shear, density structure, or spectral wave phase.

### 1.4 Subgrid topography theory

Regular mode represents each computational cell by one bed elevation \(z_b\). Subgrid mode replaces that single value with precomputed relations derived from a high-resolution DEM:

1. **Cell hypsometry:** for each coarse cell, compute cumulative ponded volume \(V(\zeta)\) from subgrid elevations. During the simulation, continuity updates \(V\), then the table is inverted/interpolated to obtain \(\zeta\).
2. **Face conveyance:** for each coarse face, precompute the minimum/maximum subgrid elevation, average conveyance depth, representative \(g n^2\), and wet fraction as functions of water level.
3. **Partial wetting:** the current NetCDF subgrid path stores wet-fraction information (`uv_pwet`). Legacy binary subgrid input assumes wet fractions of one, which can behave differently in partially wet cells.
4. **Fully wet extrapolation:** above the maximum tabulated level, SFINCS extrapolates using the fully wet cell area and the tabulated depth/roughness state.

The physical motivation is that flooding through a coarse cell is controlled by subgrid channels, levees, and depressions that a cell-average DEM cannot represent. The method can retain much of the high-resolution flood extent and conveyance while running on a coarser computational grid. It does not add new degrees of freedom: unresolved hydraulic jumps, detailed levee overtopping jets, and subgrid-scale circulation remain parameterized.

### 1.5 Boundary-condition theory

For `msk = 2`, SFINCS prescribes a time-varying water level at the boundary cell but computes the exchange flux with a weakly reflective characteristic-style relation. In one horizontal dimension, the incoming characteristic combines boundary and interior information through \(u \pm 2\sqrt{gh}\). The source implementation forms an incoming velocity estimate from the difference between the instantaneous and filtered boundary level, then combines it with the interior departure from the filtered level. This is why an open boundary should not be treated as a fixed wall with a prescribed depth only.

Other boundary meanings are simpler:

- `msk = 3` is an outflow boundary with zero water depth at the boundary cell.
- `msk = 5` estimates downstream river level from an interior cell and a user-supplied water-surface slope over a specified or computed distance.
- `msk = 6` imposes a zero normal water-level-gradient (Neumann) condition using the adjacent active cell.

SnapWave boundary theory is different: bulk \(H_s\), \(T_p\), direction, and spreading are converted into directional wave-energy spectra at support points. The stationary solver propagates wave energy over bathymetry, including refraction/shoaling and dissipation, and returns wave forces to SFINCS. It is therefore appropriate for bulk buoy observations. By contrast, `.bzi` is a free-surface time-series boundary perturbation and requires phase-resolved incoming wave information.

### 1.6 Process options

| Process | Main switch/files | Notes |
|---|---|---|
| Tide, surge, mean-water offset | `bndfile` + `bzsfile`, `bcafile`, or `netbndbzsbzifile` | Applied to `msk = 2` cells with a weakly reflective boundary flux (`bndtype = 1`). |
| Incoming wave oscillations | `bzifile` or `zi` in `netbndbzsbzifile` | Specify the incoming component only; SFINCS computes reflection. |
| Stationary wave transformation/setup | `snapwave = 1` plus SnapWave files | Couples wave heights and radiation-stress-derived forces to SFINCS. |
| Dynamic wave generation | `wavemaker_*` files | Internal generating/absorbing wave paddle; more specialized than ordinary flood forcing. |
| Wind | `wndfile`, `amufile`/`amvfile`, `netamuamvfile`, `spwfile`, or `netspwfile` | Uniform, rectilinear gridded, FEWS NetCDF, or tropical-cyclone spiderweb input. |
| Atmospheric pressure | `ampfile`, `netampfile`, `spwfile`, or `netspwfile` | Requires `baro = 1` (default) for pressure effects. |
| Rainfall | `precipfile`/`prcfile`, `amprfile`, `netamprfile`, or spiderweb rainfall | Uniform rates are in mm/hr and converted internally to m/s. |
| Infiltration | `qinf`, `qinffile`, `scsfile`, `sefffile`, `psifile`, Horton files, or NetCDF | Active only when precipitation is active; methods should not be stacked. |
| River/source discharge | `srcfile` + `disfile`, or `netsrcdisfile` | Positive discharge adds water to the snapped cell; negative values extract water. |
| Subgrid topography | `sbgfile` | Uses high-resolution DEM-derived volume and conveyance tables; NetCDF is recommended. |
| Mesh refinement | `qtrfile` | Quadtree NetCDF or binary mesh; a NetCDF quadtree can also carry elevation and mask. |
| Structures | `thdfile`, `weirfile`, `drnfile` | Thin dams, weirs, and drainage links modify fluxes. |
| Fast static inundation | `bathtub = 1` | Disables the full dynamic update; useful only for simplified screening. |

SFINCS is not a replacement for a full baroclinic ocean model. The standard mode is 2D and hydrostatic, with no salinity/temperature transport or density-driven circulation. A non-hydrostatic pressure-correction option (`nonh`) exists for suitable short-wave/tsunami applications, but it should be treated as a specialized option and validated for the target problem.

## 2. Where the model logic lives

SFINCS is implemented in Fortran under [`source/src/`](source/src/). The important control flow is:

```text
sfincs.f90
└── sfincs_lib.f90
    ├── sfincs_initialize()
    │   ├── sfincs_input.f90          read_sfincs_input()
    │   ├── sfincs_meteo.f90          read_meteo_data()
    │   ├── sfincs_domain.f90         initialize_domain()
    │   │   ├── initialize_mesh()
    │   │   ├── initialize_bathymetry()
    │   │   ├── initialize_boundaries()
    │   │   ├── initialize_roughness()
    │   │   └── initialize_hydro()
    │   ├── sfincs_structures.f90     read_structures()
    │   ├── sfincs_boundaries.f90     read_boundary_data(), find_boundary_indices()
    │   ├── sfincs_discharges.f90     read_discharges()
    │   ├── sfincs_snapwave.f90       couple_snapwave(), if snapwave = 1
    │   └── sfincs_output.f90/.F90    initialize_output()
    ├── sfincs_update()
    │   ├── update_meteo_fields(), update_meteo_forcing()
    │   ├── update_infiltration_map()
    │   ├── update_boundaries()
    │   ├── update_discharges()
    │   ├── update_wave_field(), if SnapWave is enabled
    │   ├── sfincs_momentum.f90       compute_fluxes()
    │   ├── compute_fluxes_over_structures()
    │   ├── compute_nonhydrostatic(), if enabled
    │   └── sfincs_continuity.f90     compute_water_levels()
    └── sfincs_finalize()             closes and writes final output
```

Useful source files for a new user:

| File | Main responsibility |
|---|---|
| `source/src/sfincs_input.f90` | Parses every keyword and filename in `sfincs.inp`. This is the authoritative list of supported options for this checkout. |
| `source/src/sfincs_domain.f90` | Builds regular or quadtree meshes, reads masks/elevation, classifies boundary cells, initializes roughness and hydrodynamics. |
| `source/src/sfincs_quadtree.F90` | Reads and searches quadtree meshes. |
| `source/src/sfincs_subgrid.F90` | Reads old binary or current NetCDF subgrid tables. |
| `source/src/sfincs_boundaries.f90` | Reads and updates water-level, astronomical, downstream-river, Neumann, and wave boundary conditions. |
| `source/src/sfincs_meteo.f90` | Reads/interpolates wind, pressure, and rainfall; computes wind stress and net precipitation. |
| `source/src/sfincs_discharges.f90` | Reads and updates point sources/sinks and drainage structures. |
| `source/src/sfincs_momentum.f90` | Computes face fluxes, process accelerations, wet/dry handling, and the CFL timestep. |
| `source/src/sfincs_continuity.f90` | Updates water level or subgrid volume from flux divergence and source terms. |
| `source/src/sfincs_snapwave.f90` and `source/src/snapwave/` | Coupled stationary incident/infragravity wave solver. |
| `source/src/sfincs_ncinput.F90` | FEWS-style NetCDF forcing readers and quadtree NetCDF variable readers. |
| `source/src/sfincs_output.f90`, `sfincs_ncoutput.F90` | ASCII/binary/NetCDF output. |
| `source/src/sfincs_bmi.f90` | Basic Model Interface/XMI coupling functions. |

The repository contains build and cluster scripts in `source/build_scripts/` and `source/sigrid/`, but **it does not contain a general Python preprocessing package for grids or forcing**. The supported external preprocessing route is HydroMT-SFINCS. Manual text and API examples are in `docs/example.rst`, `docs/input.rst`, and `docs/input_forcing.rst`.


## 3. Preparing the computational grid

### 3.1 Coordinate system, datum, and domain extent

1. Prefer a projected Cartesian CRS, such as the appropriate UTM zone, and record it with `epsg`. The current source also has a geographic-coordinate path (`crsgeo = 1`), but most examples and forcing readers assume projected `x`/`y` in metres; use the geographic path only after targeted validation.
2. Use one horizontal CRS for the grid, boundary locations, source locations, observation points, structures, and gridded meteorological input.
3. Use one vertical reference for bed elevation, `zsini`/`inifile`, water-level boundaries, and validation data. SFINCS does not know whether the datum is MSL, NAVD88, NAP, a local chart datum, etc.; datum consistency is the modeler's responsibility.
4. For coastal inundation, make the active domain cover the floodplain and shallow nearshore region that must respond dynamically. The manual recommends placing the seaward boundary around the shallow swash/nearshore transition (approximately 2 m depth) when offshore hydrodynamics are not the target. Do not apply a coastal water-level condition directly at an arbitrary coastline unless that is physically intended.
5. Remove or deactivate areas that do not need to evolve. Masking inactive ocean/land cells is a major reason SFINCS is fast.

### 3.2 Regular mesh

A regular SFINCS grid is staggered, rectilinear, and optionally rotated. The main grid keywords are:

```text
x0       = 318650       # x coordinate of the lower-left grid edge
y0       = 5040000      # y coordinate of the lower-left grid edge
mmax     = 250          # number of cells in x
nmax     = 157          # number of cells in y
dx       = 50           # m
dy       = 50           # m
rotation = 27           # degrees anticlockwise from east
epsg     = 32633
```

Elevation is supplied at cell centers in `depfile`; positive values are upward, so topography is positive and bathymetry is negative relative to the chosen vertical datum.

### 3.3 Active-area and boundary mask

`mskfile` is required. Its cell codes are interpreted by `initialize_boundaries()` as follows:

| Mask value | Meaning | Required companion data |
|---:|---|---|
| `0` | Inactive cell; no fluxes are computed. | None. |
| `1` | Normal active computational cell. | Elevation or subgrid table. |
| `2` | Prescribed water-level boundary cell. | `bndfile` plus `bzsfile` and/or `bcafile`, or `netbndbzsbzifile`. Without forcing, the code warns and uses 0 m. |
| `3` | Outflow boundary; water depth is kept at zero at the boundary cell. | None. |
| `5` | Downstream river boundary using interior level plus a prescribed slope/distance. | `bdrfile`. |
| `6` | Lateral/coastal Neumann boundary, \(\partial\zeta/\partial n = 0\). | An adjacent active `msk = 1` cell. |

For a conventional coastal model, make the seaward active cells `msk = 2`. Keep enclosed or irrelevant areas at `0`. Use `1` for the dynamically computed floodplain and channels. Only use `3`, `5`, and `6` when their physical meaning is intended.

### 3.4 ASCII versus binary regular-grid files

Set `inputformat = asc` while learning or debugging. ASCII `depfile` and `mskfile` contain the complete rectangular grid: `nmax` rows, each with `mmax` values, ordered consistently with the grid definition.

Set `inputformat = bin` for production-sized grids. In binary mode, `indexfile` is also required. The current readers use Fortran unformatted stream files: `indexfile` stores the number of indexed cells followed by their 1-based flattened indices, and field arrays are read in SFINCS' internal cell order. Generating these files by hand is error-prone; use HydroMT-SFINCS or a small writer that is checked against the source reader and a tiny test grid.

### 3.5 Regular, subgrid, and quadtree choices

- **Regular bathymetry:** provide `depfile` and `mskfile`. This is the easiest first model and the best route for learning.
- **Subgrid bathymetry:** provide `sbgfile`. High-resolution DEM data are preprocessed into water-level/volume and representative-depth/conveyance relations. The current NetCDF subgrid format is recommended; old binary subgrid files remain for compatibility. In subgrid mode, `depfile` is not actively used, and friction from `manningfile` is ignored because roughness is already represented in the subgrid tables.
- **Quadtree refinement:** provide `qtrfile` to use locally refined cells. NetCDF quadtree files contain topology and may contain `z` and `mask`; with a NetCDF quadtree, a separate `depfile` can be omitted because elevation can come from the mesh variable `z`.

HydroMT-SFINCS is the recommended preprocessing path. A typical component-based setup is:

```python
from hydromt_sfincs import SfincsModel

sf = SfincsModel(root="my_sfincs_model", mode="w+")
sf.grid.create_from_region(
    region={"geom": "data/region.geojson"},
    res=50.0,
    rotated=True,
    crs=32633,
)
sf.elevation.create(elevation_list=[{"elevtn": "my_dem_or_bathymetry"}])
sf.mask.create_active(...)
sf.mask.create_boundary(...)
# Optional but commonly used:
# sf.subgrid.create(...)
# sf.roughness.create(...)
# sf.infiltration.create_cn(...)
sf.write()
```

The exact method names and arguments are maintained externally; consult the [HydroMT-SFINCS documentation](https://deltares.github.io/hydromt_sfincs/latest/) for the installed version. If preprocessing manually, create a tiny domain first and inspect the log for the expected number of active water-level and flux points.


## 4. Main input file

`sfincs.inp` is the only mandatory control file and is read from the current working directory. It uses `keyword = value` lines. File paths are relative to the run directory unless absolute paths are supplied. All main simulation times use:

```text
yyyymmdd HHMMSS
```

A minimal regular-grid, dynamic coastal model could use:

```text
# ---- horizontal grid ----
x0              = 318650
y0              = 5040000
mmax            = 250
nmax            = 157
dx              = 50.0
dy              = 50.0
rotation        = 0.0
epsg            = 32633
latitude        = 52.0       # used for Coriolis in projected coordinates
crsgeo          = 0

# ---- simulation period ----
tref            = 20230901 000000
tstart          = 20230901 000000
tstop           = 20230903 000000
tspinup         = 0.0

# ---- static input ----
inputformat     = asc
depfile         = sfincs.dep
mskfile         = sfincs.msk
# indexfile is required only for inputformat = bin on a regular grid

# ---- hydrodynamic boundary forcing ----
bndfile         = sfincs.bnd
bzsfile         = sfincs.bzs

# ---- optional point discharge forcing ----
srcfile         = sfincs.src
disfile         = sfincs.dis

# ---- optional uniform meteo ----
# precipfile is accepted as a legacy alias; prcfile is read first
precipfile      = sfincs.prcp
wndfile         = sfincs.wnd

# ---- physics and stability ----
advection       = 1
alpha           = 0.5
theta           = 1.0
huthresh        = 0.05
hmin_cfl        = 0.1
dtmax           = 60.0
manning         = 0.04
rhoa            = 1.25
rhow            = 1024.0
baro            = 1
zsini           = 0.0

# ---- output ----
outputformat    = net
dtout           = 3600
dtmaxout        = 86400
dthisout        = 600
storecumprcp    = 1
storevelmax     = 1
storemeteo      = 1
obsfile         = sfincs.obs
```

For a first run, omit optional files until the basic tide-only model runs. Then add river, rainfall, wind/pressure, and waves one at a time so that any unexpected volume or momentum tendency can be attributed to one process.

### 4.1 Time conventions

- `tref` is the reference time for all ordinary ASCII forcing time columns.
- `tstart` and `tstop` are converted internally to seconds relative to `tref`.
- `tref` should normally equal `tstart` for a hand-built event model, but this is not required.
- `.bzs`, `.dis`, `.wnd`, `.prcp`, `.bzi`, and SnapWave ASCII time-series files use **seconds relative to `tref`**.
- FEWS-style NetCDF forcing uses a `time:units` attribute such as `minutes since 1970-01-01 00:00:00.0 +0000`; the source converts it to seconds relative to `tref`.
- Keep all source data in UTC. If local files use local standard/daylight time, convert them before writing SFINCS files.

### 4.2 Input file inventory

| File keyword | Required? | Contents and units | Format |
|---|---|---|---|
| `sfincs.inp` | Always | Model control, grid, times, process switches, output settings. | ASCII |
| `depfile` | Regular mode; optional if a NetCDF quadtree contains `z`; not used in subgrid mode | Bed elevation at cell centers, m above model vertical datum. | ASCII or binary |
| `mskfile` | Always unless carried in a NetCDF quadtree | Active/boundary mask. | ASCII or binary |
| `indexfile` | Regular mesh with `inputformat = bin` | Number and flattened indices of indexed cells. | Fortran stream binary |
| `qtrfile` | Quadtree only | Quadtree topology; NetCDF may include `z` and `mask`. | Binary or NetCDF |
| `sbgfile` | Subgrid mode | Cell hypsometry and face conveyance/roughness tables. | NetCDF recommended; legacy binary supported |
| `manningfile` | Optional regular mode | Spatial Manning \(n\), s m\(^{-1/3}\); ignored in subgrid mode. | Grid binary |
| `inifile` | Optional | Spatial initial water level, m. Binary; current source also accepts a quadtree NetCDF with variable `zs`. | Binary or NetCDF |
| `rstfile` | Optional | Restart state; takes precedence over `inifile` and `zsini`. | SFINCS restart binary |
| `obsfile` | Optional | Point output locations/names. | ASCII |
| `crsfile` | Optional | Discharge output cross-sections. | Delft3D Tekal-like ASCII |
| `thdfile`, `weirfile`, `drnfile` | Optional structures | Thin dams, weirs, drainage links. | ASCII |


## 5. Boundary and forcing files

### 5.1 Coastal water levels: `.bnd` and `.bzs`

This is the normal route for a measured NOAA water-level time series. The mask selects the model cells that receive the forcing; `.bnd` supplies one or more geographic forcing support points; `.bzs` supplies water level at those points.

`sfincs.bnd` has one `x y` pair per support point, in the model CRS:

```text
# x y
400000.0 1200000.0
480000.0 1250000.0
```

`sfincs.bzs` has time in seconds since `tref`, followed by one water-level column per `.bnd` row. Values are metres above the model vertical datum:

```text
# time_s  zs_at_bnd_1_m  zs_at_bnd_2_m
0.0       0.50           0.75
3600.0    0.60           0.80
7200.0    0.45           0.85
```

Important implementation details from `sfincs_boundaries.f90`:

1. Each `msk = 2` model cell is associated with the two nearest `.bnd` support points.
2. For distances \(d_1\) and \(d_2\), the code computes `fac = d2 / (d1 + d2)` and then uses
   \(\eta_b = fac\,\eta_1 + (1-fac)\,\eta_2\). With one support point, the same value is used everywhere.
3. Values are linearly interpolated in time. Before the first supplied time the first value is held constant; after the last supplied time the last value is held constant. The reader also warns if the forcing does not cover the whole run. Do not rely on endpoint padding: provide a buffer before `tstart` and after `tstop` when practical.
4. The boundary is weakly reflective rather than a rigid Dirichlet clamp at every face. The water-level value is imposed at the boundary cell, while a characteristic/Riemann-style update computes the boundary flux using the prescribed level, filtered mean level, interior level, and \(\sqrt{gh}\). This allows outgoing disturbances to leave with less spurious reflection than a hard wall.
5. Water level is never allowed below the local bed elevation (or below the subgrid minimum elevation in subgrid mode).

For a single NOAA station, use one row in `.bnd` and two columns in `.bzs`: time and station water level. The station should be geographically representative of the model's open-boundary water level. If the station is inside an estuary and the boundary is offshore, a regional tide/surge model or tidal reconstruction may be more appropriate than directly assigning the interior gauge offshore.

### 5.2 Astronomical tide: `.bca`

A Delft3D-style `.bca` file can supply harmonic constituents instead of, or in addition to, a slowly varying `.bzs` record. The reader expects one `[forcing]` block per `.bnd` point, with the same number/order of blocks as `.bnd` rows. Data rows contain constituent name, amplitude, and phase; phase is supplied in degrees and converted internally. `A0` may be included as a mean offset and is inserted as zero if omitted.

Illustrative structure:

```text
[forcing]
Name     = boundary_1
Function = astronomic
Quantity = water level
Unit     = m
A0       0.00   0.0
M2       0.45  132.4
S2       0.11   88.2
N2       0.09  104.0
K1       0.12  201.0
O1       0.10  185.0
```

The nodal correction and astronomical argument are computed at `tref`. The update adds

\[
\eta_{astro}(t)=\sum_c A_c\cos(\omega_c t-\phi_c)
\]

to the time-interpolated `.bzs` value. Therefore:

- If `.bzs` already contains **total measured water level**, do **not** also supply tidal constituents; that would double-count tide.
- Use `.bca` when you want a harmonic/predicted tide and put only a non-tidal residual, storm-surge anomaly, or zero in `.bzs`.
- For NOAA data, a common choice is `bzs = observed total water level` and no `.bca`. Another valid choice is `bzs = observed - predicted` and `.bca` built from local harmonic constituents, but the constituent generation is then a separate tidal-analysis step.

### 5.3 NetCDF water-level input

Set `netbndbzsbzifile` instead of separate ASCII boundary files when using the FEWS-style format. The current reader expects dimensions `time` and `stations`, variables `x`, `y`, `time`, and `zs`, and optional `zi`. The `time:units` attribute must follow the FEWS convention, normally `minutes since 1970-01-01 00:00:00.0 +0000`. Coordinates must be in the same CRS as the SFINCS grid. This route is useful when many stations are needed or when provenance metadata matters, but ASCII is simpler for one NOAA station.


### 5.4 River and point discharges: `.src` and `.dis`

Point discharge forcing is a source term in continuity, not a momentum boundary condition. Each location is snapped to an active SFINCS cell.

`sfincs.src` has one point per row in model coordinates:

```text
# x y
300000.0 1500000.0
380000.0 1650000.0
```

`sfincs.dis` has seconds since `tref` and one discharge column per `.src` point, in \(m^3\,s^{-1}\):

```text
# time_s  Q_src_1_m3_s  Q_src_2_m3_s
0.0       100.0         1000.0
3600.0    300.0         1100.0
7200.0    0.0           1300.0
```

The code applies \(Q\Delta t/A\) to the cell water level in regular mode or \(Q\Delta t\) to cell volume in subgrid mode. Positive values add water. A negative value extracts water, subject to available-volume checks. The reader warns if the series does not cover the whole simulation; provide complete coverage.

If the source data are in cubic feet per second, convert to SI before writing `.dis`:

```text
Q_m3_s = Q_ft3_s * 0.028316846592
```

USGS parameter `00060` is commonly streamflow in cubic feet per second, but Codex must inspect the local file's own metadata rather than inferring units from the filename.

For a tidally influenced river mouth, avoid imposing both an upstream discharge at `srcfile` and an inconsistent downstream water level in a way that creates artificial storage or recirculation. Put `.src` at the upstream inflow location and use `msk = 2` at the coastal/downstream water-level boundary. A `msk = 5` downstream river boundary is a different option: it derives downstream level from an interior point using a slope and distance from `bdrfile`; it does not read a discharge time series.


### 5.5 Wave forcing: choose the right representation

SFINCS has two conceptually different wave inputs:

1. **Bulk spectral boundary conditions for SnapWave.** Use measured significant wave height, peak/dominant period, mean direction, and directional spreading. This is the appropriate route for ordinary buoy or hindcast wave products when `snapwave = 1`.
2. **Incoming free-surface oscillations.** `bzifile` is a high-frequency, approximately zero-mean incoming water-level component added to the coastal boundary. Do not create this by simply treating bulk significant wave height as water level. Generating a defensible incoming wave/IG free-surface series requires spectral or phase-resolving preprocessing and is outside a simple table conversion.

For NOAA/NDBC-style wave data, the usual mapping is:

| Source field | SFINCS/SnapWave use | Units and convention |
|---|---|---|
| `WVHT` or significant wave height | `hs` in SnapWave | metres |
| `DPD` or dominant/peak period | `tp` in SnapWave | seconds |
| `MWD` or mean wave direction | `wd` in SnapWave | degrees clockwise from true north, direction waves come from; source converts to Cartesian going-to radians |
| Directional spread | `ds` in SnapWave | degrees; estimate from spectral data if absent |
| raw free-surface elevation | possible input to a separate incoming-wave analysis for `bzi` | metres; do not confuse with `Hs` |

ASCII SnapWave time-series files are:

```text
# sfincs.inp
snapwave            = 1
snapwave_bndfile    = waves.bnd
snapwave_bhsfile    = waves.bhs
snapwave_btpfile    = waves.btp
snapwave_bwdfile    = waves.bwd
snapwave_bdsfile    = waves.bds
```

`waves.bnd` contains `x y` support points. Each of the other four files contains `time_s` followed by one value per wave support point:

```text
# waves.bhs: time_s Hs_point1 Hs_point2 [m]
0.0     1.20 1.10
3600.0  1.45 1.32

# waves.btp: time_s Tp_point1 Tp_point2 [s]
0.0     9.0  8.6
3600.0  9.4  9.1

# waves.bwd: time_s direction_point1 direction_point2 [deg N, coming-from]
0.0     125.0 130.0
3600.0  132.0 136.0

# waves.bds: time_s spreading_point1 spreading_point2 [deg]
0.0     25.0 25.0
3600.0  28.0 28.0
```

All four time columns must be identical. For one representative buoy, use one support point and one value after time. A single-point alternative is `snapwave_jonswapfile`, whose rows are `time_s Hm0 Tp direction spreading water_level`, but the separate files map more transparently to buoy fields.

The FEWS-style alternative is `netsnapwavefile`, with dimensions `time` and `stations` and variables `x`, `y`, `time`, `hs`, `tp`, `wd`, and `ds`.

SnapWave solves a stationary, directionally distributed wave-energy balance with linear shoaling/refraction, breaking dissipation, bottom friction, optional wind growth, vegetation dissipation, and an infragravity-wave energy balance. SFINCS receives mapped wave forces at the hydrodynamic update interval. Because SnapWave computes wave setup from wave transformation, do not independently add a wave-setup estimate to `.bzs` unless the selected configuration explicitly disables the corresponding force and the double-counting risk has been documented.


### 5.6 Wind, pressure, and rainfall

Uniform ASCII files are the easiest way to start:

```text
# sfincs.wnd: time_s wind_speed_m_s direction_deg_N_coming_from
0.0     5.0  120.0
3600.0  15.0 180.0
7200.0  10.0 165.0

# sfincs.prcp: time_s rainfall_mm_hr
0.0     0.0
3600.0  15.0
7200.0  10.0
```

The wind stress is computed from squared wind speed, the air/water density ratio, and a speed-dependent drag coefficient. The default drag table is \(C_d=0.001\) at low wind, \(0.0025\) at 28 m/s, and \(0.0015\) at 50 m/s.

Spatially varying Delft3D ASCII meteo files are `amufile`/`amvfile` for eastward/northward wind, `amprfile` for precipitation, and `ampfile` for absolute air pressure. FEWS-style NetCDF alternatives are:

- `netamuamvfile`: `time`, `x`, `y`, `eastward_wind`, `northward_wind`;
- `netamprfile`: `time`, `x`, `y`, `Precipitation` or `precipitation`;
- `netampfile`: `time`, `x`, `y`, `barometric_pressure`;
- `netspwfile`: tropical-cyclone spiderweb with `wind_x`, `wind_y`, absolute `pressure`, and optional `precipitation`.

Note the pressure convention: ASCII spiderweb pressure is converted from pressure drop using `gapres`; NetCDF spiderweb pressure is expected to be absolute pressure. Gridded `.amp` and `netampfile` pressure are also absolute. Set `baro = 0` only if pressure-gradient forcing should be disabled deliberately.

### 5.7 Infiltration and initial/restart state

Infiltration is a sink in net precipitation and is activated only when precipitation is active. The supported method families are:

- `qinf`: uniform constant rate, entered as positive mm/hr and applied only above `qinf_zmin` (default 0 m);
- `qinffile` or `infiltration_type = c2d`: spatially varying constant rate in mm/hr;
- `scsfile` or `infiltration_type = cna`: Curve Number without recovery;
- `smaxfile`, `sefffile`, `ksfile` or `infiltration_type = cnb`: Curve Number with recovery;
- `psifile`, `sigmafile`, `ksfile` or `infiltration_type = gai`: Green-Ampt;
- `f0file`, `fcfile`, `kdfile` or `infiltration_type = hor`: Horton decay with recovery.

Do not provide multiple methods and expect them to combine; choose one. For quadtree models, non-uniform infiltration must use the NetCDF `infiltration_file` route.

Initial water level can be uniform (`zsini`), spatially varying (`inifile`), or a full restart (`rstfile`). SFINCS raises the initial level to bed elevation where the prescribed value is below the bed, so dry cells can be represented without negative depths. For event simulations, a short `tspinup` can damp the transition from `zsini` to boundary forcing, but a physically consistent restart is preferable for long tide/river equilibration.


## 6. Boundary-condition generation workflow

This section is written as an implementation checklist for Codex or another local coding agent. The existing field data should be inspected locally; do not assume that every NOAA product has identical column names, flags, datum, or units.

### 6.1 Build a forcing inventory first

Create a manifest that maps every local input to its SFINCS destination before writing model files:

```yaml
model:
  tref_utc: "2023-09-01T00:00:00Z"
  epsg: 32633
  vertical_datum: "NAVD88"       # example only; use the DEM datum
  run_directory: "./sfincs_run"

water_level:
  local_file: "data/noaa_water_level.csv"
  station_id: "..."              # copy from file metadata
  station_lon_lat: [null, null]  # fill after inspecting metadata
  source_time_zone: "UTC"        # verify; convert if not UTC
  source_datum: "..."            # verify: MSL, MLLW, NAVD88, STND, etc.
  source_units: "m"              # verify; convert ft to m if needed
  output_points: "sfincs.bnd"
  output_timeseries: "sfincs.bzs"

river:
  local_file: "data/river_discharge.csv"
  source_units: "ft3/s or m3/s" # verify from metadata
  inflow_lon_lat: [null, null]
  output_points: "sfincs.src"
  output_timeseries: "sfincs.dis"

waves:
  local_file: "data/noaa_wave_data.txt"
  source_units: "SI"            # verify WVHT/DPD/MWD and missing values
  support_lon_lat: [null, null]
  mode: "snapwave_bulk"         # do not select bzi from bulk Hs alone
```

Codex should fail loudly if the time zone, datum, CRS, or units cannot be established. Silent unit/datum conversion is the largest source of forcing errors.

### 6.2 NOAA water level to `sfincs.bnd`/`sfincs.bzs`

Use observed total water level when the goal is to reproduce a measured event. The target is a time series of total still-water level at the coastal support point, in metres above the model datum.

Recommended preprocessing sequence:

1. **Read the local water-level file and print its header/metadata.** Identify station ID, latitude/longitude, product, time zone, datum, unit, sampling interval, and quality columns.
2. **Normalize timestamps to timezone-aware UTC (`pandas.Timestamp(..., tz="UTC")`).** NOAA CO-OPS API requests can specify `time_zone=gmt`; local standard/daylight time must be converted. Drop or resolve duplicated timestamps.
3. **Normalize units to metres.** If source values are feet, multiply by `0.3048`.
4. **Normalize the vertical datum.** Prefer requesting or reading the data directly in the model datum when the station supports it. If that is impossible, apply a station-specific datum offset from authoritative benchmark metadata. Never use a generic MLLW-to-MSL value.
5. **Quality-control values.** Remove missing sentinels and values flagged as bad by the source. Preserve plausible storm surges; do not automatically delete large values as outliers. Investigate isolated spikes and step changes.
6. **Interpolate only small gaps.** Linear interpolation is reasonable for short gaps in a 6-minute tide record. Long gaps should be filled using a defensible tidal prediction/regional model or excluded by shortening the simulation; document the choice.
7. **Place the support point.** Convert station longitude/latitude to model `x y` using the same EPSG as the grid. If the gauge is not dynamically representative of the open boundary, use a regional tide/surge product or document the transfer assumption.
8. **Cover the model period.** Write forcing from before `tstart` to after `tstop`. For a hindcast starting at an arbitrary observed state, also set `zsini` or an `inifile` consistently with the initial boundary level.
9. **Write plain numeric files.** Do not include a header unless the local reader skips it. The Fortran examples use whitespace-separated numeric rows.

Minimal writer logic:

```python
import numpy as np
import pandas as pd
from pyproj import Transformer

def seconds_since_tref(time_utc, tref_utc):
    t = pd.to_datetime(time_utc, utc=True)
    tref = pd.Timestamp(tref_utc).tz_convert("UTC")
    return (t - tref).total_seconds().to_numpy(dtype=float)

def write_bnd_bzs(station_lon, station_lat, epsg, time_utc, eta_m,
                  bnd_path, bzs_path):
    transformer = Transformer.from_crs("EPSG:4326", f"EPSG:{epsg}", always_xy=True)
    x, y = transformer.transform(station_lon, station_lat)
    np.savetxt(bnd_path, np.array([[x, y]]), fmt="%.3f %.3f")
    t = seconds_since_tref(time_utc, "MODEL_TREF_UTC")
    out = np.column_stack([t, eta_m])
    np.savetxt(bzs_path, out, fmt="%.3f %.6f")
```

Replace `"MODEL_TREF_UTC"` with the actual `tref`. The reader does not need evenly spaced data, but strictly increasing times are required for robust interpolation. A measured 6-minute record can be written directly; averaging it to an hourly series may smooth high-water timing and should not be done without a reason.

### 6.3 Choosing between total water level and tidal harmonics

Use this decision rule:

- **Hindcast/validation against the NOAA gauge:** put observed total water level in `.bzs`; omit `bcafile`.
- **Forecast or scenario with a predicted tide:** put astronomical constituents in `.bca` and put storm-surge anomaly or zero in `.bzs`.
- **Forecast with a complete regional ocean forecast:** put total forecast water level in `.bzs`; omit `bcafile`.

If NOAA harmonic constituents are available, convert them to the exact local model datum before writing `.bca`. Constituent phase depends on time reference and convention; verify by reconstructing a known NOAA prediction at the station before attaching the file to SFINCS.


### 6.4 River discharge to `sfincs.src`/`sfincs.dis`

Use discharge as a volumetric source at the upstream edge of the modeled river network:

1. Identify the station coordinates and confirm that the station measures the river entering the SFINCS domain, not a distributary or downstream tidal section already represented by the water-level boundary.
2. Convert timestamps to UTC and discharge to \(m^3\,s^{-1}\).
3. Inspect sign convention. River inflow should normally be positive. Negative values in a source record mean extraction in SFINCS and should only be retained if physically intended.
4. Choose the `.src` point just inside the active domain where the discharge enters. Project the point to model CRS. Ensure that the point snaps to the intended river cell; plotting the snapped cell over the DEM is essential for narrow rivers.
5. Write `.dis` with the same number of value columns as `.src` rows and the same time base as `.bzs`.
6. If multiple tributaries exist, write one `.src` row per tributary and one `.dis` column per row in the same order.

A single-point pair looks like:

```text
# sfincs.src
512345.000  4123456.000

# sfincs.dis
0.0       85.0
3600.0    92.0
7200.0    101.0
```

If the river discharge data are missing during part of the event, reconstruct them using a hydrologic model or interpolate only over short, hydrologically plausible gaps. Do not fill long data gaps with zeros unless zero flow is supported by the station record.

### 6.5 Wave data to SnapWave files

Use bulk wave properties as SnapWave boundary conditions when wave-induced setup, runup contribution, or wave-driven currents matter. The simplest one-buoy implementation is:

```text
# sfincs.inp additions
snapwave         = 1
dtwave           = 3600
snapwave_bndfile = waves.bnd
snapwave_bhsfile = waves.bhs
snapwave_btpfile = waves.btp
snapwave_bwdfile = waves.bwd
snapwave_bdsfile = waves.bds
```

Preprocessing requirements:

1. Convert buoy time to UTC and then seconds since `tref`.
2. Convert wave height to metres if needed. NDBC standard-meteorological data normally report `WVHT` in metres, `DPD` and `APD` in seconds, and `MWD` in degrees clockwise from true north.
3. Treat source-specific missing values such as `MM`, `99`, `999`, or `9999` as missing only after checking the local file format.
4. Use dominant/peak period (`DPD`) for `tp` when available. Average period (`APD`) is a fallback, not an identical quantity.
5. Keep mean wave direction in the coming-from nautical convention. Do not pre-convert to Cartesian components for the ASCII `bwd` file; `snapwave_boundaries.f90` performs the conversion.
6. Provide directional spreading in degrees. If the data set lacks spreading, derive it from directional spectral moments when available. Otherwise use a documented, site-appropriate constant and run a sensitivity test; do not silently treat zero as a realistic spreading.
7. Locate `waves.bnd` at the buoy or at a representative offshore point in the model CRS. The SnapWave grid/domain must cover the propagation path from that support point to the coastline. A deep-water buoy far outside the model should not be copied unchanged to a shallow boundary without a wave-transformation assumption.

For one wave point, each file has exactly two columns:

```text
# waves.bnd
485000.000  1234500.000

# waves.bhs
0.0     1.20
3600.0  1.45

# waves.btp
0.0     9.0
3600.0  9.4

# waves.bwd
0.0     125.0
3600.0  132.0

# waves.bds
0.0     25.0
3600.0  28.0
```

SnapWave boundary time series are linearly interpolated. Keep `Hs` non-negative, `Tp` positive, directions in `[0, 360)`, and spreading positive. Check that wave times cover the full SFINCS simulation.

Do **not** write `bzifile` from `Hs` alone. `bzi` expects an incoming free-surface elevation perturbation in metres at high temporal resolution, normally derived from a wave spectrum or phase-resolving data with the reflected component excluded. If only bulk `Hs/Tp/direction` are available, use SnapWave.


### 6.6 Automated checks Codex should perform after writing files

The local forcing-generation script should end with checks that fail the workflow on structural errors:

```python
def assert_sfincs_ascii_series(path, nvalue_columns, t0_s, t1_s):
    a = np.loadtxt(path, ndmin=2)
    assert a.ndim == 2 and a.shape[1] == nvalue_columns + 1
    assert np.all(np.isfinite(a))
    assert np.all(np.diff(a[:, 0]) > 0), f"non-increasing time in {path}"
    assert a[0, 0] <= t0_s and a[-1, 0] >= t1_s, f"insufficient time coverage in {path}"
    return a

def assert_matching_points(point_path, series_path):
    points = np.loadtxt(point_path, ndmin=2)
    series = np.loadtxt(series_path, ndmin=2)
    assert points.shape[0] == series.shape[1] - 1
```

Additional physical checks:

- `.bzs`: values should be plausible total water levels in metres relative to the selected datum. Compare summary statistics and a plot against the raw gauge record.
- `.bnd`: plot support points, `msk = 2` cells, station location, and coastline in model coordinates.
- `.src`/`.dis`: integrate \(Q\) over time and compare the written-file volume with the source-data volume. Plot the snapped inflow cell on the elevation/mask grid.
- SnapWave files: verify that all four time columns are bit-for-bit identical and that direction/spreading conventions are documented.
- All files: check line counts, monotonic time, no NaN, no tabular header row, and exact column counts.

Generate a forcing report containing maps and time-series overlays. The report should show raw versus converted data so unit, time-zone, and datum mistakes are visible before the model is run.

### 6.7 Recommended files for the current field-data case

For a first dynamically meaningful configuration with the data described by the user:

| Field data | Recommended SFINCS input | `sfincs.inp` entries |
|---|---|---|
| NOAA total water level at one station | One-point `.bnd`/`.bzs` | `bndfile = sfincs.bnd`, `bzsfile = sfincs.bzs`; no `bcafile` |
| River discharge | One or more `.src`/`.dis` points | `srcfile = sfincs.src`, `disfile = sfincs.dis` |
| Bulk wave height/period/direction/spreading | SnapWave boundary files | `snapwave = 1`, `snapwave_bndfile`, `snapwave_bhsfile`, `snapwave_btpfile`, `snapwave_bwdfile`, `snapwave_bdsfile` |
| Uniform rainfall, if available | `.prcp` | `precipfile = sfincs.prcp` or `prcfile = sfincs.prcp` |
| Uniform wind, if available | `.wnd` | `wndfile = sfincs.wnd` |

Start with water level only, then add river discharge, then waves. A tide-only run should reproduce boundary water levels and produce plausible wetting without unbounded depth. A river-only perturbation should conserve the imposed volume approximately, allowing for water crossing open boundaries. A wave-enabled run should be compared against a run with identical water level and discharge but `snapwave = 0`.


## 7. Running, outputs, and quality control

### 7.1 Running

Run from the directory containing `sfincs.inp`:

```bash
sfincs |& tee sfincs_log.txt
```

With Docker on Linux/macOS:

```bash
docker pull deltares/sfincs-cpu
docker run --rm -v "$PWD":/data deltares/sfincs-cpu
```

To build this source checkout with GNU Fortran on Linux:

```bash
cd /workspace/source
autoreconf -vif
./autogen.sh
./configure FC=gfortran FCFLAGS="-fopenmp -O3 -fallow-argument-mismatch -w" \
            FFLAGS="-fopenmp -O3 -fallow-argument-mismatch -w" --disable-shared
make -j4
```

See `source/README.md` and `source/build_scripts/` for platform-specific builds.

### 7.2 Output files

With `outputformat = net`:

- `sfincs_map.nc`: gridded static and time-varying output. Typical variables include `zb`, `msk`, `zs`, `h`, `u`, `v`, `zsmax`, `hmax`, `vmax`, `qmax`, `cumprcp`, and `cuminf`, subject to process/output switches.
- `sfincs_his.nc`: point, cross-section, runup-gauge, and drainage time series when corresponding geometries are supplied.
- `sfincs.rst`: restart file when `dtrstout` or `trstout` requests it.
- Log output: process switches, files read, active-cell counts, boundary/discharge snapping, timestep information, and warnings.

Inspect `zb` and `msk` in `sfincs_map.nc` immediately after the first run. This verifies how SFINCS interpreted the grid before analyzing the hydrodynamics.

### 7.3 Validation hierarchy

1. **Input geometry:** map bed elevation, mask, active cells, boundary points, and discharge points.
2. **Input time series:** plot raw source data and written SFINCS forcing on the same UTC/model-datum axes.
3. **Tide-only run:** compare modeled and observed water level at the gauge/observation point; do not expect exact agreement if the gauge is inland of the boundary, but timing and trend should be plausible.
4. **River-only sensitivity:** compare cumulative inflow volume against \(\int Q\,dt\).
5. **Wave sensitivity:** compare `snapwave = 0` and `snapwave = 1`; inspect wave height and wave-force fields before interpreting flood levels.
6. **Stability:** check maximum water depth and average timestep. Enable `timestep_analysis = 1` if instability appears.
7. **Mass balance:** precipitation plus river input plus boundary inflow should be balanced by storage change, infiltration, and boundary outflow within numerical/output tolerances.


## 8. Common pitfalls

1. **Datum mismatch:** `.dep`, `.bzs`, `zsini`, and observations must use the same vertical datum. A NOAA MLLW record cannot be combined directly with a NAVD88 DEM without a station-specific conversion.
2. **Time-zone mismatch:** all ASCII forcing times are seconds since `tref`; source timestamps should be normalized to UTC before conversion.
3. **Double-counting tide:** total observed water level in `.bzs` plus harmonic tide in `.bca` duplicates the tide.
4. **Double-counting wave setup:** a wave-setup estimate in `.bzs` plus SnapWave wave forces can duplicate setup.
5. **Bulk waves confused with free surface:** `Hs` is not a water-level time series. Use SnapWave files for bulk wave properties; reserve `.bzi` for an incoming free-surface perturbation derived by wave analysis.
6. **Mask/forcing mismatch:** `.bnd`/`.bzs` affect only `msk = 2` cells. If no such cells exist, the log warns that no open boundary points were found.
7. **Discharge snapped to the wrong cell:** visually inspect the `.src` location and active river channel, especially at coarse resolution.
8. **Wrong rain unit:** uniform and Delft3D rainfall inputs are mm/hr; SFINCS converts them to m/s internally.
9. **Wrong wind/wave direction convention:** uniform `.wnd` and SnapWave `bwd` use degrees clockwise from north for the direction the flow/waves come from.
10. **Subgrid roughness:** when `sbgfile` is active, a separate `manningfile` is ignored because subgrid tables carry representative roughness/conveyance.
11. **Old defaults in documentation:** this checkout uses `advection = 1`, `alpha = 0.5`, `theta = 1.0`, and `dtmax = 60 s`; check `source/src/sfincs_input.f90` for the active source defaults.
12. **Whitespace formatting:** use spaces around `=` in `sfincs.inp`. The manual warns that mixed tabs/spaces can prevent a keyword from being recognized.

## 9. References and further documentation

- Leijnse, T., van Ormondt, M., Nederhoff, K., and van Dongeren, A. (2021). *Modeling compound flooding in coastal systems using a computationally efficient reduced-physics solver.* Coastal Engineering, 163, 103796. <https://doi.org/10.1016/j.coastaleng.2020.103796>
- Bates, P. D., Horritt, M. S., and Fewtrell, T. J. (2010). A simple inertial formulation of the shallow water equations for efficient two-dimensional flood inundation modelling. *Journal of Hydrology*. <https://doi.org/10.1016/j.jhydrol.2010.03.027>
- Van Ormondt et al. (2025). Subgrid corrections for the linear inertial equations of a compound flood model. *Geoscientific Model Development*. <https://doi.org/10.5194/gmd-18-843-2025>
- SFINCS manual in this repository: `docs/overview.rst`, `docs/input.rst`, `docs/input_forcing.rst`, `docs/parameters.rst`, and `docs/output.rst`.
- HydroMT-SFINCS: <https://deltares.github.io/hydromt_sfincs/latest/>
- NOAA CO-OPS data API: <https://api.tidesandcurrents.noaa.gov/api/prod>
- NOAA NDBC observation descriptions: <https://www.ndbc.noaa.gov/obsdes.shtml>

## 10. Copy-ready handoff prompt for Codex

Use the following operational instruction with the local data directory:

> Inspect the field-data files and identify their timestamp column/time zone, coordinate metadata, units, vertical datum, missing-value codes, and QC flags. Normalize all times to UTC and choose the SFINCS `tref`. Normalize water level to metres above the DEM datum, river discharge to m3/s, wave height to m, wave period to s, and wave direction to degrees clockwise from north coming-from. Generate `sfincs.bnd/.bzs` from total water level, `sfincs.src/.dis` from discharge, and SnapWave `waves.bnd/.bhs/.btp/.bwd/.bds` from bulk wave data. Do not generate `.bca` from total observed water level and do not generate `.bzi` from bulk `Hs`. Validate monotonic times, finite values, exact column counts, model-period coverage, projected coordinates, and raw-versus-converted plots. Update `sfincs.inp` only after the files pass these checks.

