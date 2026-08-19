# Hidden Mechanical Information beyond the Static Limit

## Demonstration code

- `HMI_demo.m` — evaluates static and harmonic observability, hidden activation, and low-frequency opening from a compatible tetrahedral mesh.
- `prepare_mesh_sphere.m` — generates the sphere example mesh.
- `prepare_mesh_cube.m` — generates the cube example mesh.
- `prepare_mesh_torus.m` — generates the torus example mesh.
- `prepare_mesh_spiky_virus.m` — generates the symmetric spiky-particle example mesh.
- `prepare_mesh_nanostar.m` — generates the five-arm nanostar example mesh.

## Usage

The code follows a two-step **mesh preparation → analysis** pattern.

```matlab
prepare_mesh_sphere;
results = HMI_demo('mesh_sphere.mat');
```

The other supplied geometries are used in the same way.

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

Each mesh-preparation routine writes a MAT file containing a scalar structure named `model`.
`HMI_demo` reads this structure, assembles the mixed finite-element system, evaluates the observable and hidden dimensions, computes the hidden activation spectrum, and writes numerical results and figures to a separate output folder.

## Requirements

- MATLAB with `decomposition` and `exportgraphics` support.
- Gmsh command-line executable for the mesh-preparation routines.
- No additional MATLAB toolboxes are required.
- The analysis code uses only base MATLAB sparse matrix operations, factorization, and plotting utilities.

Gmsh can be supplied through the system PATH, the `GMSH_EXE` environment variable, or directly as a f