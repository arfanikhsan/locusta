# locusta 🦗

My Grasshopper practice: a step-by-step path from **intermediate** to **advanced** parametric design in Rhino + Grasshopper.

Every exercise has its own folder with a brief, learning goals, key components, a step checklist and stretch goals. I add my definition, screenshots and notes as I go.

## Roadmap

Status: ⬜ not started · 🟨 in progress · ✅ done

| # | Exercise | Level | Focus | Plugins | Status |
|---|----------|-------|-------|---------|--------|
| 01 | [Attractor Grid](exercises/1-intermediate/01-attractor-grid) | Intermediate | Distance-based remapping (Remap Numbers, Bounds), Multiple attractors (closest point, weighted sums) | — | ⬜ |
| 02 | [Data Tree Management](exercises/1-intermediate/02-data-trees) | Intermediate | Reading tree paths with Param Viewer, Graft, Flatten, Simplify, Flip Matrix | — | ⬜ |
| 03 | [Surface Paneling](exercises/1-intermediate/03-surface-paneling) | Intermediate | Isotrim / SubSrf workflow, Surface normals and frames | LunchBox (optional) | ⬜ |
| 04 | [Voronoi Patterns](exercises/1-intermediate/04-voronoi-pattern) | Intermediate | Populate / random seeds, Voronoi cells and offsets | Weaverbird (optional) | ⬜ |
| 05 | [Image Sampler Facade](exercises/1-intermediate/05-image-sampler-facade) | Intermediate | Image Sampler setup, Mapping colour channels to geometry | — | ⬜ |
| 06 | [Parametric Tower](exercises/1-intermediate/06-parametric-tower) | Intermediate | Series and transforms per floor, Lofting from section curves | — | ⬜ |
| 07 | [Mesh Fundamentals](exercises/1-intermediate/07-mesh-fundamentals) | Intermediate | Mesh vertices, faces, normals, Mesh colouring by analysis | Weaverbird | ⬜ |
| 08 | [Loops & Recursion](exercises/2-advanced/08-loops-recursion) | Advanced | Loop start/end and data passing, Recursive subdivision rules | Anemone | ⬜ |
| 09 | [Form Finding with Kangaroo](exercises/2-advanced/09-kangaroo-form-finding) | Advanced | Goals: Anchor, Length, Load, Pressure, Mesh relaxation | Kangaroo2 (bundled), Weaverbird | ⬜ |
| 10 | [Optimisation with Galapagos / Wallacei](exercises/2-advanced/10-optimisation) | Advanced | Defining genes and fitness, Single-objective (Galapagos) | Galapagos (bundled), Wallacei or Octopus | ⬜ |
| 11 | [Environmental Analysis](exercises/2-advanced/11-environmental-analysis) | Advanced | EPW weather files, Radiation analysis on facades | Ladybug Tools | ⬜ |
| 12 | [Scripting with Python](exercises/2-advanced/12-scripting-python) | Advanced | RhinoCommon geometry API, Type hints and list/tree access | — | ⬜ |
| 13 | [Agent-Based Systems](exercises/2-advanced/13-agent-systems) | Advanced | Agent state and update rules, Separation, alignment, cohesion | Culebra (optional), Multipipe (bundled R7+) | ⬜ |
| 14 | [Differential Growth](exercises/2-advanced/14-differential-growth) | Advanced | Curve subdivision during simulation, Collision avoidance | Kangaroo2, Anemone | ⬜ |
| 15 | [Digital Fabrication](exercises/2-advanced/15-fabrication) | Advanced | Sectioning / waffle notches, Unrolling developable strips | OpenNest (optional) | ⬜ |
| 16 | [Clusters, User Objects & Clean Definitions](exercises/2-advanced/16-clusters-userobjects) | Advanced | Clusters with inputs/outputs, Publishing user objects (.ghuser) | — | ⬜ |
| 17 | [Capstone: Pavilion](exercises/2-advanced/17-capstone-pavilion) | Advanced | Form finding + paneling + analysis + fabrication in one pipeline, Managing a large definition | As needed | ⬜ |

### Intermediate (01–07)
Data trees, attractors, paneling, patterns, and meshes. These are the core skills that every later exercise builds on.

### Advanced (08–17)
Loops, physics simulation, optimisation, environmental analysis, scripting, agent systems, fabrication, and keeping large definitions clean. Ends with a capstone pavilion that combines everything.

## Structure

```
locusta/
├── exercises/
│   ├── 1-intermediate/     # 01–07
│   ├── 2-advanced/         # 08–17
│   └── _template/          # Blank template for extra exercises
├── scripts/
│   ├── python/             # Reusable Python script components
│   └── csharp/             # Reusable C# script components
└── references/             # Tutorials, papers, inspiration
```

Each exercise folder:

```
NN-name/
├── README.md        # Brief, goals, steps, notes
├── definition.gh    # Grasshopper definition
├── model.3dm        # Rhino model (if needed)
└── images/          # Canvas + viewport screenshots (result.png = cover)
```

## Workflow

1. Open the exercise README and read the brief.
2. Build the definition and save it as `definition.gh` in the folder.
3. Tick off steps (`[ ]` → `[x]`) and update the status in the table above.
4. Save a screenshot as `images/result.png`.
5. Write down what you learned.

To add an extra exercise:

```bash
cp -r exercises/_template exercises/2-advanced/18-my-exercise
```

## Plugins used

Install these from [Food4Rhino](https://www.food4rhino.com/) or the Rhino 8 Package Manager (`_PackageManager`):

- **Kangaroo2** (bundled): physics and form finding
- **Galapagos** (bundled): optimisation
- **Weaverbird**: mesh subdivision
- **Anemone**: loops
- **Ladybug Tools**: environmental analysis
- **Wallacei** or **Octopus**: multi-objective optimisation
- Optional: LunchBox, Culebra, OpenNest

## Tips

- Save as `.ghx` if you want readable diffs; `.gh` is smaller but binary.
- Internalise referenced geometry, or commit the `.3dm` next to the `.gh`, so the file opens cleanly.
- For large `.3dm` files, use [Git LFS](https://git-lfs.com/): `git lfs track "*.3dm"`.

## Environment

- Rhino 8 with Grasshopper
