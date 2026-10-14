---
type: entity
projects: [godt]
created: 2026-07-18
updated: 2026-07-18
status: current
---

# (3D) Shepp-Logan Phantom

Standard synthetic test object for tomographic reconstruction algorithms: an ellipsoid-composed digital phantom originally from X-ray CT literature, extended to 3D (MATLAB File Exchange implementation by M. Schabel, used by [[2411-08423-godt]]). In [[2411-08423-godt]]'s simulation, real-valued RI 1.33–1.335, 40 μm max size, used to generate 180 synthetic object waves via the [[brenner1993|wave propagation method]] and validate that GODT reconstruction matches the finite-difference derivative of the known ground truth.

## Relevance

Per [[godt]]: the standard synthetic-phantom choice for any forward-model/reconstruction validation in this project's simulation thread — likely the first phantom to reuse for multi-shear extension testing before moving to more realistic phantoms.
