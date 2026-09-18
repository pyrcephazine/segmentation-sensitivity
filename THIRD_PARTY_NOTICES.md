# Source provenance

The study replaces the physical model and workflow of multihepatic
(https://github.com/pyrcephazine/multihepatic), starting revision
2283b8eda631fe3d62931abeed35c0c8e2817f74. The unused original application is
not needed to run this repository.

`src/segmentation_sensitivity/network_mesh.py` is incorporated from networks_fenicsx
(https://github.com/scientificcomputing/networks_fenicsx), revision
4d0dd39d6a788185c589c9f4c5017701cd8d27db. It includes this study's changes:
explicit MPI cell partitioner, optional omission of multiplier submeshes, and
orientation-independent edge-color lookup. These changes are part of the
local source file; no upstream checkout, installation, or patch is required.

FEniCSx, PETSc, MPI, and fenicsx_ii are external, unmodified dependencies.

## License for network_mesh.py

Copyright 2025 Jørgen S. Dokken, Cécile Daversin-Catty, Joseph P. Dean

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
