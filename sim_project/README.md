# sim_project — HPGe gamma-spectroscopy simulation

A Geant4 simulation of a mono-energetic gamma point source illuminating a
high-purity germanium (HPGe) crystal. For every energy-depositing step inside
the crystal it records the **interaction position** and **deposited energy**.
That hit list is the input for the
[SolidStateDetectors.jl](https://github.com/JuliaPhysics/SolidStateDetectors.jl)
signal simulations in `~/ssd_projects` (`ppc`, `ppc_2`), which turn the hits
into detector pulses and, ultimately, an energy spectrum.

```
   122 keV gamma point source (Geant4)
            │
            ▼
   energy deposits in the Ge crystal  ──►  output/hits.csv
            │                              event,x,y,z,edep
            ▼
   SolidStateDetectors.jl (ppc / ppc_2)  ──►  pulses ──► energy spectrum
```

## Detector geometry

| Item            | Value                                              |
|-----------------|----------------------------------------------------|
| Active crystal  | Germanium cylinder, **r = 15 mm, height = 10 mm**  |
| Crystal frame   | axis = z, **N+ face z = 0, P+ point contact z = 10 mm** |
| World           | 20 cm air cube                                     |

The crystal dimensions and coordinate frame match the SSD `config.yaml` files
in `ppc` / `ppc_2` (`tube: r=15, h=10, origin z=5`; N+ contact on the bottom
face, P+ point contact on the top face). A hit at Geant4 world position
`(x, y, z)` in mm therefore maps **directly** onto an SSD
`CartesianPoint(x, y, z) * 1e-3` in metres — no transform needed.

> A real HPGe also has a cryostat, an Al end-cap window and an n+ dead layer.
> These are intentionally omitted for now; add them in `DetectorConstruction.cc`
> as extra volumes when needed (they do not change the crystal frame).

## Source

`G4GeneralParticleSource`, fully macro-controlled. Default (`run.mac`):

- **122 keV** mono-energetic gammas (the dominant Co-57 line);
- **point source** on the detector axis at **z = -30 mm**, i.e. 30 mm outside
  the **N+ contact face** (z = 0), so gammas enter through the N+ side;
- emission restricted to a **cone aimed at the crystal**, half-angle
  `atan(15/30) = 26.57°` — every primary gamma is sent toward the N+ face,
  so no primaries are wasted. The GPS angular frame is rotated with
  `rot1/rot2` so the cone opens toward +z (the crystal).

For a true 4π point source instead, replace the cone block with a single
`/gps/ang/type iso`. Change the energy / position / cone angle by editing the
`/gps/...` lines.

## Physics

`PhysicsList` registers `G4EmStandardPhysics_option4` — the most accurate EM
option, with low-energy models for the photo-electric effect, Compton and
Rayleigh scattering, needed for realistic HPGe spectra. Production cut: 0.1 mm.

No decay physics is included (the source is a plain gamma). To simulate a
radioactive isotope later, register `G4DecayPhysics` +
`G4RadioactiveDecayPhysics` in `PhysicsList.cc`.

## Build

```bash
source /root/geant4/g4_install/bin/geant4.sh        # Geant4 environment
cd /root/geant4/sim_project/build
cmake -DCMAKE_PREFIX_PATH=/root/geant4/g4_install ..
make -j$(nproc)
```

## Run

The program writes `hits.csv` into the current working directory, so run it
from `output/` to collect results there:

```bash
source /root/geant4/g4_install/bin/geant4.sh
cd /root/geant4/sim_project/output
../build/sim_project ../build/run.mac      # batch -> writes output/hits.csv
```

Interactive visualization (Qt OpenGL viewer):

```bash
cd /root/geant4/sim_project/build
./sim_project                              # opens the viewer, fires 300 gammas
```

## Output: `hits.csv`

One row per energy-depositing step inside the crystal:

| Column     | Unit | Meaning                                          |
|------------|------|--------------------------------------------------|
| `event`    | —    | Geant4 event id = one primary gamma              |
| `x_mm`     | mm   | interaction x (SSD crystal frame)                |
| `y_mm`     | mm   | interaction y                                    |
| `z_mm`     | mm   | interaction z (0 = N+ face, 10 = P+ face)        |
| `edep_keV` | keV  | energy deposited in that step                    |

### Feeding it into SSD

Every gamma is prompt, so all hits sharing an `event` id are one detector
pulse. Group the rows by `event` and build one SSD `Event` per Geant4 event:

```julia
using DelimitedFiles, SolidStateDetectors, Unitful

rows  = readdlm("output/hits.csv", ',', skipstart = 1)  # event,x,y,z,edep
ev    = Int.(rows[:, 1])
for e in unique(ev)
    sel      = findall(ev .== e)
    pts      = [CartesianPoint{Float32}(rows[i,2]*1e-3,   # mm -> m
                                        rows[i,3]*1e-3,
                                        rows[i,4]*1e-3) for i in sel]
    energies = [rows[i,5] * u"keV" for i in sel]
    evt = Event(pts, energies)
    simulate!(evt, sim)
    # pulse amplitude -> one count in the energy spectrum
end
```

Building this loader/bridge inside `ppc` and `ppc_2` is the next step.

## Notes

- The run is **serial** (one thread) so `hits.csv` is written by a single
  stream and cannot be corrupted. Switching to multithreaded mode requires a
  different output strategy (per-thread files, a mutex, or `G4AnalysisManager`).
- `hits.csv` is overwritten by every run — rename it to keep results.
- With the cone source every gamma is aimed at the crystal; ~70 % of them
  interact in 10 mm of Ge at 122 keV (the rest pass straight through).
- A photo-electric (full-energy) event produces several sub-0.1 mm hits in one
  `event`. The SSD bridge may want to merge hits closer than the detector's
  position resolution into a single charge cloud.
