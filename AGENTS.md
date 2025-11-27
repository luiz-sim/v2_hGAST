Purpose

This repository contains the hGAST-based modelling of a floating wind turbine with catenary moorings. Automated agents (e.g. Copilot / Codex) may be used to extend and maintain the code. This document explains the key design choices and how agents should interact with the codebase, with special focus on the mooring repositioning control feature.

Relevant files and modules

truss.f90
Implements the mooring line model as 3D truss elements (catenary chains) and the mooring time integration (MOORINGS_tr).

truss.inp
Input file defining mooring line elements, nodes, rest lengths, diameters, etc.

repositioning.inp (pending)
New input file describing the desired floater repositioning trajectory in sway (time vs target offset) and controller gains.

Other standard hGAST input files (paths.inp, dfile_el.inp, etc.) should not be changed by agents unless explicitly requested.

Repositioning control – overview

The repository implements an optional mooring-based repositioning controller. When enabled, it:

Reads a time series of desired sway offsets from repositioning.inp.

Uses a PI-type controller to adjust the rest length of the top few truss elements in one mooring line.

Moves the floater laterally (in sway) by exploiting mooring asymmetry, without introducing artificial external forces.

Leaves core behaviour unchanged when repositioning.inp is not present.

Key aspects:

A boolean repos_active indicates whether repositioning.inp was successfully read.

Arrays t_repos(:) and y_repos(:) store the time series (time, y_target).

Controller parameters repos_Kp, repos_Ki, repos_vwinch define a PI law and maximum rate of rest-length change.

Arrays ALENG0_tr, ALENG_ctrl, ALENG_min, ALENG_max track original and current rest lengths and bounds for each truss element.

winch_elem(:) identifies the small set of elements (top segments of one mooring line) that are actuated.

The actual mooring element formulations (stiffness, forces) are unchanged; they simply use the possibly updated body_tr(e)%ALENG_tr values.

Control logic

At each time step in the mooring integrator:

update_repositioning(TTIME, DT) is called.

It interpolates the current target sway y_target_curr from the time series.

It computes position error e = y_target_curr - (y_float - y_ref) using:

y_float: current floater sway DOF in global coordinates (set externally).

y_ref: sway position at t=0 (reference).

A PI law computes a desired total rest-length rate dL_total_dt, rate-limited by repos_vwinch.

The incremental change dL_step = dL_total_dt * DT is distributed across active winch elements, respecting per-element bounds.

For each winch element, ALENG_ctrl and body_tr(e)%ALENG_tr are updated accordingly.

If at any point no elements can be shortened/lengthened further (all at bounds), update_repositioning does nothing; the mooring system simply holds the floater where it can.

Expectations for future agents

When modifying or extending this feature:

Do not change default behaviour when repositioning.inp is absent.

Prefer to extend the existing control interfaces (init_repositioning, update_repositioning, etc.) rather than introducing ad-hoc logic elsewhere.

Keep controller-related parameters and arrays within the mooring/truss module; avoid scattering repositioning logic across unrelated modules.

If adding 2D control (surge+sway) or more advanced controllers:

Extend repositioning.inp format carefully (e.g. add an optional x_target column).

Maintain backward compatibility with existing files that only specify y_target.

When changing the node or element indexing scheme, update winch_elem selection logic accordingly and ensure that the chosen elements correspond to the top segments of a specific mooring line.

When possible, add or update small test cases (numerical, not unit tests) that:

Run with and without repositioning.inp.

Exercise a simple repositioning manoeuvre (e.g. 0 → 50 m → 0 in sway) and check:

that the floater moves in the correct direction,

that rest lengths remain within [ALENG_min, ALENG_max],

and that the solver does not diverge.

Coding style

Preserve existing Fortran style (indentation, naming conventions).

Prefer small, clearly named subroutines for new logic (init_repositioning, init_winch_elements, update_repositioning, set_floater_sway) rather than large monolithic blocks.

Avoid duplicating logic (e.g. interpolation routines); reuse helper functions where available.

Add brief comments near new code explaining:

inputs/outputs,

how it interacts with existing modules,

and any assumptions (e.g., hard-coded choice of mooring line for winching).