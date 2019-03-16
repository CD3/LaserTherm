---
title : LaserTherm
---

# Overview

# Purpose

- High Fidelity Laser-Tissue interaction Simulations
- Knowledge Base
- Teaching
- Experimentation

# Goals

- Build on top of small, simple, independent, swappable pieces.
    - A student should be able to implement a new heat solver without knowing how the configuration engine works.
    - A simulation should be able to use an implicit or explicit heat solver.
- Favor simple over complex at the bottom.
    - Low-level components can *assume* certain conditions are true so simplify their implementation.
    - High-level components deal with cases were low-level component assumptions are not true.
- Provide a library for outside clients.
- Write unit tests.
- Write derivations.


# Design

- Use Standard Library, Boost, and Eigen
- Callbacks implemented with Boost.Signals2
-

