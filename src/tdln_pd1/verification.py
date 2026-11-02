"""Executable numerical and directional verification for the TDLN twin."""

from __future__ import annotations

import numpy as np
import pandas as pd

from .simulator import Dose, STATE_NAMES, TwinParameters, simulate
from .spatial_dynamics import (SPECIES, ZONE_NAMES, SpatialParameters, simulate_spatial,
                               simulate_spatial_coupled)


def _relative_error(a: float, b: float) -> float:
    return float(abs(a - b) / max(abs(b), 1e-12))


def verify_twin() -> dict:
    """Run prespecified software invariants and return a machine-readable report."""
    checks: list[dict] = []

    baseline = simulate(days=21)
    array = baseline.loc[:, STATE_NAMES].to_numpy(dtype=float)
    finite_nonnegative = bool(np.isfinite(array).all() and (array >= 0).all())
    checks.append({"name": "systemic_finite_nonnegative", "passed": finite_nonnegative,
                   "minimum_state": float(array.min())})

    rejected = []
    for layer, call in (
        ("systemic", lambda: simulate(days=1, parameters=TwinParameters(tumour_growth=-0.1))),
        ("spatial", lambda: simulate_spatial(
            days=1, parameters=SpatialParameters(antigen_decay=-0.1),
        )),
    ):
        try:
            call()
        except ValueError:
            rejected.append(layer)
    checks.append({"name": "input_parameter_domain_guards",
                   "passed": set(rejected) == {"systemic", "spatial"},
                   "rejected_invalid_layers": rejected,
                   "expected": "negative parameters rejected before integration"})

    zero = simulate(days=21, doses=(Dose(0, 0.0), Dose(14, 0.0)))
    zero_difference = np.abs(zero.loc[:, STATE_NAMES].to_numpy() - array)
    zero_scale = np.maximum(np.abs(array), 1.0)
    zero_error = float(np.max(zero_difference / zero_scale))
    checks.append({"name": "zero_dose_limit", "passed": zero_error <= 1e-6,
                   "maximum_scaled_relative_error": zero_error,
                   "maximum_absolute_error": float(zero_difference.max()),
                   "relative_tolerance": 1e-6,
                   "note": "Zero-dose event segmentation may restart the adaptive solver."})

    dose_amounts = (0.0, 0.25, 0.5, 1.0, 2.0)
    dose_rows = []
    for amount in dose_amounts:
        trajectory = baseline if amount == 0 else simulate(
            days=21, doses=(Dose(0, amount), Dose(14, amount)),
        )
        dose_rows.append({"dose": amount,
                          "final_tumour_cells": float(trajectory.tumour_cells.iloc[-1]),
                          "peak_pd1_occupancy": float(trajectory.pd1_occupancy.max())})
    tumour = np.array([row["final_tumour_cells"] for row in dose_rows])
    occupancy = np.array([row["peak_pd1_occupancy"] for row in dose_rows])
    monotone = bool((np.diff(tumour) <= 0).all() and (np.diff(occupancy) >= 0).all())
    checks.append({"name": "dose_response_direction", "passed": monotone,
                   "expected": "tumour non-increasing and occupancy non-decreasing", "series": dose_rows})

    strict = simulate(days=21, doses=(Dose(0, 1.0), Dose(14, 1.0)), rtol=1e-9, atol=1e-11)
    standard = simulate(days=21, doses=(Dose(0, 1.0), Dose(14, 1.0)))
    convergence = {
        name: _relative_error(float(standard[name].iloc[-1]), float(strict[name].iloc[-1]))
        for name in ("tumour_cells", "tpex_tdln", "effector_tumour", "pd1_occupancy")
    }
    max_convergence_error = max(convergence.values())
    checks.append({"name": "systemic_solver_tolerance_convergence",
                   "passed": max_convergence_error <= 1e-4,
                   "relative_endpoint_errors": convergence,
                   "maximum_relative_error": max_convergence_error, "tolerance": 1e-4})

    transport = simulate_spatial(days=7, reactions=False)
    totals = transport.groupby(["day", "species"]).amount.sum().unstack()
    drift = {
        species: float(np.max(np.abs(totals[species] / totals[species].iloc[0] - 1.0)))
        for species in SPECIES
    }
    max_drift = max(drift.values())
    checks.append({"name": "spatial_transport_mass_conservation", "passed": max_drift <= 1e-6,
                   "relative_drift_by_species": drift, "maximum_relative_drift": max_drift,
                   "tolerance": 1e-6})

    capacity_parameters = SpatialParameters(
        antigen_decay=0.0, antigen_uptake=0.0, dc_death=0.0,
        priming_rate=0.0, tpex_growth=2.0, tpex_carrying_capacity=1.0e5,
        tpex_differentiation=0.0, tpex_death=0.0,
        tpex_to_medulla=0.0, tpex_to_follicle=0.0,
    )
    capacity_initial = pd.DataFrame([
        {"species": species, "zone": zone, "amount": (
            1.0 if species == "antigen" and zone == "paracortex" else
            8.0e4 if species == "tpex" and zone == "paracortex" else 0.0
        )}
        for species in SPECIES
        for zone in ZONE_NAMES
    ])
    capacity_trajectory = simulate_spatial(
        days=20, step=0.2, parameters=capacity_parameters, initial=capacity_initial,
    )
    capacity_tpex = capacity_trajectory[
        capacity_trajectory.species.eq("tpex")
        & capacity_trajectory.zone.eq("paracortex")
    ].amount
    maximum_capacity_ratio = float(
        capacity_tpex.max() / capacity_parameters.tpex_carrying_capacity
    )
    final_capacity_ratio = float(
        capacity_tpex.iloc[-1] / capacity_parameters.tpex_carrying_capacity
    )
    capacity_passed = maximum_capacity_ratio <= 1.0002 and abs(final_capacity_ratio - 1) <= 2e-4
    checks.append({
        "name": "spatial_tpex_niche_capacity",
        "passed": bool(capacity_passed),
        "maximum_capacity_ratio": maximum_capacity_ratio,
        "final_capacity_ratio": final_capacity_ratio,
        "tolerance": 2e-4,
        "note": "Numerical saturation invariant under a high-antigen stress scenario.",
    })

    control_system = simulate(days=14)
    treated_system = simulate(days=14, doses=(Dose(0, 1.0),))
    control_spatial = simulate_spatial_coupled(control_system, days=14)
    treated_spatial = simulate_spatial_coupled(treated_system, days=14)
    def endpoint(frame, species, zone):
        return float(frame[(frame.species == species) & (frame.zone == zone)].amount.iloc[-1])
    spatial_direction = {
        "control_efferent_effector": endpoint(control_spatial, "effector", "efferent_lymph"),
        "treated_efferent_effector": endpoint(treated_spatial, "effector", "efferent_lymph"),
        "control_paracortex_tpex": endpoint(control_spatial, "tpex", "paracortex"),
        "treated_paracortex_tpex": endpoint(treated_spatial, "tpex", "paracortex"),
    }
    coupled_passed = (spatial_direction["treated_efferent_effector"] >
                      spatial_direction["control_efferent_effector"] and
                      spatial_direction["treated_paracortex_tpex"] >
                      spatial_direction["control_paracortex_tpex"])
    checks.append({"name": "coupled_checkpoint_direction", "passed": bool(coupled_passed),
                   **spatial_direction})

    spatial_dose_rows = []
    zero_spatial = None
    for amount in dose_amounts:
        system = simulate(days=14, doses=(Dose(0, amount),))
        spatial = simulate_spatial_coupled(system, days=14)
        if amount == 0:
            zero_spatial = spatial
        spatial_dose_rows.append({
            "dose": amount,
            "peak_pd1_occupancy": float(system.pd1_occupancy.max()),
            "final_efferent_effector": endpoint(
                spatial, "effector", "efferent_lymph"
            ),
            "final_paracortex_tpex": endpoint(spatial, "tpex", "paracortex"),
        })
    assert zero_spatial is not None
    spatial_zero_difference = np.abs(
        zero_spatial.amount.to_numpy(dtype=float)
        - control_spatial.amount.to_numpy(dtype=float)
    )
    spatial_zero_scale = np.maximum(
        np.abs(control_spatial.amount.to_numpy(dtype=float)), 1.0
    )
    spatial_zero_error = float(np.max(spatial_zero_difference / spatial_zero_scale))
    spatial_occupancy = np.array([
        row["peak_pd1_occupancy"] for row in spatial_dose_rows
    ])
    spatial_effector = np.array([
        row["final_efferent_effector"] for row in spatial_dose_rows
    ])
    spatial_tpex = np.array([
        row["final_paracortex_tpex"] for row in spatial_dose_rows
    ])
    spatial_dose_passed = bool(
        spatial_zero_error <= 1e-6
        and (np.diff(spatial_occupancy) >= 0).all()
        and (np.diff(spatial_effector) >= 0).all()
        and (np.diff(spatial_tpex) >= 0).all()
    )
    checks.append({
        "name": "spatial_zero_dose_and_dose_direction",
        "passed": spatial_dose_passed,
        "maximum_zero_dose_scaled_relative_error": spatial_zero_error,
        "relative_tolerance": 1e-6,
        "expected": (
            "zero dose equivalent to control; occupancy, efferent effectors and "
            "paracortical Tpex non-decreasing over the prespecified dose series"
        ),
        "series": spatial_dose_rows,
        "note": "Software direction check under default hypotheses, not a clinical dose response.",
    })

    return {"schema_version": 1, "research_only": True,
            "passed": all(check["passed"] for check in checks), "checks": checks}
