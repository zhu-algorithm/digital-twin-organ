"""Six-zone TDLN response driven by the closed-loop liver tumour model."""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np
import pandas as pd

from .liver_tdln_coupling import simulate_liver_tdln
from .simulator import Dose, TwinParameters
from .spatial_dynamics import (SpatialInitialState, SpatialParameters,
                               simulate_spatial)


@dataclass(frozen=True)
class LiverTDLNSpatialBridge:
    antigen_drainage_gain: float = 0.45
    dc_migration_gain: float = 100.0
    # Nominal healthy-tissue flow used to form a dimensionless delivery term.
    perfusion_reference: float = 0.8657193939964243


def _validate_bridge(bridge: LiverTDLNSpatialBridge) -> None:
    values = np.asarray(list(asdict(bridge).values()), dtype=float)
    if not np.isfinite(values).all() or np.any(values < 0):
        raise ValueError("liver–spatial bridge parameters must be finite and non-negative")
    if bridge.perfusion_reference <= 0:
        raise ValueError("bridge perfusion_reference must be strictly positive")


def simulate_liver_tdln_spatial(*, days: float = 28.0, step: float = 0.25,
                                doses: tuple[Dose, ...] = (),
                                bridge: LiverTDLNSpatialBridge | None = None,
                                spatial_parameters: SpatialParameters | None = None,
                                spatial_initial: SpatialInitialState | None = None,
                                spatial_rtol: float = 1e-8,
                                spatial_atol: float = 1e-10,
                                **coupled_kwargs) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return simultaneous organ/system trajectory and its six-zone TDLN field."""
    b = bridge or LiverTDLNSpatialBridge()
    _validate_bridge(b)
    systemic = simulate_liver_tdln(days=days, step=step, doses=doses, **coupled_kwargs)
    timeline = systemic.day.to_numpy(dtype=float)
    antigen = systemic.tumour_antigen.to_numpy(dtype=float)
    perfusion_factor = (
        systemic.liver_perfusion.to_numpy(dtype=float) / b.perfusion_reference
    )
    immune_parameters = coupled_kwargs.get("immune_parameters") or TwinParameters()
    # This is the explicit influx term in the systemic mature-DC equation.
    dc_migration_flux = immune_parameters.dc_migration * antigen
    occupancy_values = systemic.pd1_occupancy.to_numpy(dtype=float)
    antigen_source = lambda time: b.antigen_drainage_gain * float(
        np.interp(time, timeline, antigen * perfusion_factor))
    dc_source = lambda time: b.dc_migration_gain * float(
        np.interp(time, timeline, dc_migration_flux * perfusion_factor))
    occupancy = lambda time: float(np.interp(time, timeline, occupancy_values))
    initial = spatial_initial or SpatialInitialState(afferent_antigen=0.0, afferent_dc=0.0)
    spatial = simulate_spatial(
        days=days, step=step, parameters=spatial_parameters, initial=initial,
        reactions=True, checkpoint_occupancy=occupancy,
        afferent_sources={"antigen": antigen_source, "mature_dc": dc_source},
        rtol=spatial_rtol, atol=spatial_atol,
    )
    spatial["pd1_occupancy"] = np.interp(spatial.day, timeline, occupancy_values)
    spatial["liver_perfusion_delivery_factor"] = np.interp(
        spatial.day, timeline, perfusion_factor)
    spatial["systemic_dc_migration_flux"] = np.interp(
        spatial.day, timeline, dc_migration_flux)
    spatial["afferent_antigen_source"] = b.antigen_drainage_gain * np.interp(
        spatial.day, timeline, antigen * perfusion_factor)
    spatial["afferent_dc_source"] = b.dc_migration_gain * np.interp(
        spatial.day, timeline, dc_migration_flux * perfusion_factor)
    spatial["liver_tumour_fraction"] = np.interp(
        spatial.day, timeline, systemic.liver_tumour_fraction)
    return systemic, spatial


def liver_tdln_spatial_counterfactual(*, doses: tuple[Dose, ...], days: float = 28.0,
                                      step: float = 0.25,
                                      **kwargs) -> pd.DataFrame:
    """Matched untreated/anti-PD-1 spatial fields with identical bridge assumptions."""
    if not doses:
        raise ValueError("spatial liver counterfactual requires an explicit dose schedule")
    frames = []
    for arm, schedule in (("control", ()), ("anti_pd1", doses)):
        _, spatial = simulate_liver_tdln_spatial(
            days=days, step=step, doses=schedule, **kwargs)
        spatial["arm"] = arm
        frames.append(spatial)
    return pd.concat(frames, ignore_index=True)

