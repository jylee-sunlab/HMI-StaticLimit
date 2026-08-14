# HMI-StaticLimit

Lightweight MATLAB demonstration code for static and harmonic mechanical observability in connected solid geometries.

## Files

- `HMI_demo.m` — universal analysis code for a compatible tetrahedral mesh MAT file.
- `prepare_mesh_sphere.m` — sphere example.
- `prepare_mesh_cube.m` — cube example.
- `prepare_mesh_torus.m` — torus example.
- `prepare_mesh_spiky_virus.m` — symmetric spiky-particle example.
- `prepare_mesh_nanostar.m` — five-arm nanostar example.

Each mesh-preparation script writes a MAT file containing a scalar structure named `model`.

## Requirements

- MATLAB.
- Gmsh command-line executable for the mesh-preparation scripts.
- No additional MATLAB toolboxes are required by the supplied analysis code.

Gmsh can be found in one of the following ways.

1. Put `gmsh` on the system PATH.
2. Set the environment variable `GMSH_EXE`.
3. Pass the full executable path to a mesh-preparation function.

Example:

```matlab
prepare_mesh_sphere('C:\path\to\gmsh.exe');
```

## Quick start

Generate a mesh.

```matlab
prepare_mesh_sphere;
```

Run the analysis.

```matlab
results = HMI_demo('mesh_sphere.mat');
```

Other supplied examples are run in the same way.

```matlab
prepare_mesh_cube;
results = HMI_demo('mesh_cube.mat');

prepare_mesh_torus;
results = HMI_demo('mesh_torus.mat');

prepare_mesh_spiky_virus;
results = HMI_demo('mesh_spiky_virus.mat');

prepare_mesh_nanostar;
results = HMI_demo('mesh_nanostar.mat');
```

## Model data contract

`HMI_demo.m` expects a MAT file containing one scalar structure named `model`.

Required fields are

```text
model.schemaVersion
model.name
model.nodes
model.elements.tet4
model.region.solidMask
model.region.fluidMask
model.region.portIndex
model.boundary.outerFaces
model.scale.a
```

Optional fields are

```text
model.ports(k).elements
model.ports(k).nominalCenter
model.meta
```

The supplied mesh-preparation scripts provide examples of this data contract.

## Output

For an input such as

```text
mesh_sphere.mat
```

the analysis creates

```text
mesh_sphere_HMI_results/
```

containing the numerical summary, MAT results, and separate figures.

Generated meshes and result folders are excluded from Git tracking by `.gitignore`.

## Scope

The supplied meshes are intentionally lightweight demonstration meshes.
They are intended for rapid reproduction of the numerical mechanism rather than high-accuracy production calculations.
