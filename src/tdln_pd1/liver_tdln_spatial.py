"""Six-zone TDLN response driven by the closed-loop liver tumour model."""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np
import pandas as pd

from .liver_tdln_coupling import simulate_liver_tdln
from .simulator import Dose
from .spatial_dynamics import (SpatialInitialState, SpatialParameters,
                               simulate_spatial)


@dataclass(frozen=True)
class LiverTDLNSpatialBridge:
    antigen_drainage_gain: float = 0.45


def _validate_bridge(bridge: LiverTDLNSpatialBridge) -> None:
    values = np.asarray(list(asdict(bridge).values()), dtype=float)
    if not np.isfinite(values).all() or np.any(values < 0):
        raise ValueError("liver–spatial bridge parameters must be finite and non-negative")


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
    occupancy_values = systemic.pd1_occupancy.to_numpy(dtype=float)
    source = lambda time: b.antigen_drainage_gain * float(np.interp(time, timeline, antigen))
    occupancy = lambda time: float(np.interp(time, timeline, occupancy_values))
    initial = spatial_initial or SpatialInitialState(afferent_antigen=0.0, afferent_dc=0.0)
    spatial = simulate_spatial(
        days=days, step=step, parameters=spatial_parameters, initial=initial,
        reactions=True, checkpoint_occupancy=occupancy,
        afferent_sources={"antigen": source}, rtol=spatial_rtol, atol=spatial_atol,
    )
    spatial["pd1_occupancy"] = np.interp(spatial.day, timeline, occupancy_values)
    spatial["afferent_antigen_source"] = b.antigen_drainage_gain * np.interp(
        spatial.day, timeline, antigen)
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

