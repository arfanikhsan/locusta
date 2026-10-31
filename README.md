# locusta

A workspace for Rhino + Grasshopper exercises: parametric definitions, scripts, and notes.

## Structure

```
locusta/
├── exercises/          # One folder per exercise
│   └── _template/      # Copy this to start a new exercise
│       ├── README.md   # Brief, approach, parameters, learnings
│       └── images/     # Screenshots / renders of the result
├── scripts/            # Reusable script components
│   ├── python/         # GhPython / Python 3 script components
│   └── csharp/         # C# script components
└── references/         # Papers, links, inspiration
```

## Starting a new exercise

```bash
cp -r exercises/_template exercises/01-my-exercise
```

Then save your files inside that folder:

- `definition.gh` — the Grasshopper definition
- `model.3dm` — the Rhino model (if referenced geometry is needed)
- `images/` — screenshots of the canvas and viewport

Naming convention: `NN-short-name` (e.g. `01-attractor-grid`, `02-voronoi-panel`).

## Exercises

| # | Name | Topics | Plugins |
|---|------|--------|---------|
| — | — | — | — |

## Tips

- Prefer saving definitions as `.ghx` (XML) if you want readable diffs; `.gh` is smaller but binary.
- Internalise referenced geometry, or commit the `.3dm` alongside the `.gh`, so the definition opens cleanly.
- List any third-party plugins (Kangaroo, Ladybug, Weaverbird, Pufferfish…) in the exercise README.
- Large `.3dm` files: consider [Git LFS](https://git-lfs.com/) (`git lfs track "*.3dm"`).

## Environment

- Rhino: 8
- Grasshopper: bundled with Rhino
