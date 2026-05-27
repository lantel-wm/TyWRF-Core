# CUDA-Ready Layout

The CPU reference implementation must use data structures that can map cleanly
to future CUDA kernels.

## State Layout Rules

- Use structure-of-arrays.
- Give each physical variable its own contiguous buffer.
- Keep explicit halo cells.
- Keep `i` as the contiguous dimension.
- Do not use nested `std::vector` in hot paths.
- Do not use arrays of structs for grid fields.

Canonical 3D indexing:

```text
idx = ((j * nz) + k) * nx + i
```

where `nx` includes the active span plus halo cells for the field view.

`include/tywrf/field_view.hpp` defines the shared kernel layout contract:

- `FieldView2D<T>` and `FieldView3D<T>` are standard-layout, trivially
  copyable kernel views. They contain only a pointer, total extents, strides,
  and halo metadata.
- `FieldLayout3D::index(i, j, k)` follows
  `idx = ((j * nz) + k) * nx + i`; `stride_i` is always `1`, `stride_k` is
  `nx`, and `stride_j` is `nx * nz`.
- Active-domain loop bounds are explicit: `i_begin()/i_end()`,
  `j_begin()/j_end()`, and `k_begin()/k_end()`.

`include/tywrf/grid.hpp` keeps the v1 grid shapes minimal but WRF-aware:

- mass fields use `mass_nx, mass_ny, mass_nz`;
- `U` is staggered in `i`;
- `V` is staggered in `j`;
- `W`, `PH`, and `PHB` use the full vertical level count;
- surface fields use the horizontal mass grid and horizontal halo.

`include/tywrf/state.hpp` owns one buffer per v1 core variable. This is the
structure-of-arrays state container; host-side `FieldStorage2D<T>` and
`FieldStorage3D<T>` wrap the owning vectors, while `State::view()` returns a
bundle of POD views suitable for CPU scalar, OpenMP, and future CUDA kernel
entry points.

## Kernel Interface Rules

Kernels receive plain views only:

```text
pointer + shape + stride + halo
```

Kernels must not perform NetCDF I/O, logging, virtual dispatch, or hidden
allocation. CPU scalar, OpenMP, and future CUDA kernels should share the same
view shape and indexing contracts.

`include/tywrf/dynamics/tendency.hpp` provides the first state tendency apply
skeleton. It accepts only `FieldView2D/3D` or `StateView` POD bundles and
applies `field += dt * tendency` across active cells. Halo cells are skipped
explicitly, and the loop order is `j-k-i` for 3D and `j-i` for 2D so `i`
remains the innermost contiguous dimension. The zero-tendency variant validates
and counts the active region without writing field storage.

## Loop Rules

- Make `i` the innermost loop.
- Keep boundary and halo updates explicit.
- Prefer tiled loops when kernels grow beyond simple smoke checks.
- Keep WRF physics staging buffers explicit and outside core kernels.
