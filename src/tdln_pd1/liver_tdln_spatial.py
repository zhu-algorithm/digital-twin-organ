"""Six-zone TDLN response driven by the closed-loop liver tumour model."""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

from .liver_tdln_coupling import (
    COMBINED_STATE_NAMES, IMMUNE_DYNAMIC_STATES, LiverTDLNCouplingParameters,
    _combined_rhs, _decorate_coupled_output, _validate_coupling,
    _validated_combined_states,
)
from .liver_tumor import (
    STATE_NAMES as LIVER_STATES, LiverTumourInitialState,
    LiverTumourParameters, _inflow as _liver_inflow,
    _validate as _validate_liver, _validate_parameters as _validate_liver_parameters,
)
from .simulator import (STATE_NAMES as IMMUNE_STATES, Dose, InitialState,
                        TwinParameters, _validate_model_inputs)
from .spatial_dynamics import (
    SPECIES, ZONE_NAMES, SpatialInitialState, SpatialParameters,
    _index as _spatial_index, _rhs as _spatial_rhs,
    _validate_parameters as _validate_spatial_parameters,
    _validated_spatial_solution, initial_vector,
)


@dataclass(frozen=True)
class LiverTDLNSpatialBridge:
    antigen_drainage_gain: float = 0.45
    dc_migration_gain: float = 100.0
    # Nominal healthy-tissue flow used to form a dimensionless delivery term.
    perfusion_reference: float = 0.8657193939964243
    # Fractional viable-rim killing per efferent effector cell and day.
    efferent_effector_killing_gain: float = 2.0e-6


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
    """Jointly solve liver, systemic immunity and the six-zone TDLN field.

    The spatial efferent-effector state feeds back into viable-rim killing in
    the same ODE solve.  A zero feedback gain recovers the one-way bridge.
    """
    b = bridge or LiverTDLNSpatialBridge()
    _validate_bridge(b)
    supported = {
        "liver_parameters", "immune_parameters", "coupling_parameters",
        "liver_initial", "immune_initial", "flow_multiplier",
        "portal_flow_multiplier", "arterial_flow_multiplier",
        "angiogenesis_inhibition", "rtol", "atol",
    }
    unknown = sorted(set(coupled_kwargs) - supported)
    if unknown:
        raise TypeError(f"unexpected coupled simulation arguments: {unknown}")
    lp = coupled_kwargs.get("liver_parameters") or LiverTumourParameters()
    ip = coupled_kwargs.get("immune_parameters") or TwinParameters()
    cp = coupled_kwargs.get("coupling_parameters") or LiverTDLNCouplingParameters()
    li = coupled_kwargs.get("liver_initial") or LiverTumourInitialState()
    ii = coupled_kwargs.get("immune_initial") or InitialState()
    flow = coupled_kwargs.get("flow_multiplier", 1.0)
    portal_flow = coupled_kwargs.get("portal_flow_multiplier")
    arterial_flow = coupled_kwargs.get("arterial_flow_multiplier")
    inhibition = coupled_kwargs.get("angiogenesis_inhibition", 0.0)
    systemic_rtol = coupled_kwargs.get("rtol", 1e-7)
    systemic_atol = coupled_kwargs.get("atol", 1e-9)
    sp = spatial_parameters or SpatialParameters()
    initial = spatial_initial or SpatialInitialState(afferent_antigen=0.0, afferent_dc=0.0)

    _validate_liver_parameters(lp)
    _validate_liver(li, "liver initial state")
    _validate_model_inputs(ip, ii)
    _validate_coupling(cp)
    _validate_spatial_parameters(sp)
    if not np.isfinite(days) or not np.isfinite(step) or days <= 0 or step <= 0:
        raise ValueError("days and step must be finite and strictly positive")
    for name, value in (("systemic rtol", systemic_rtol), ("systemic atol", systemic_atol),
                        ("spatial rtol", spatial_rtol), ("spatial atol", spatial_atol)):
        if not np.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be finite and strictly positive")
    if not np.isfinite(flow) or flow <= 0:
        raise ValueError("flow_multiplier must be finite and strictly positive")
    for name, value in (("portal_flow_multiplier", portal_flow),
                        ("arterial_flow_multiplier", arterial_flow)):
        if value is not None and (not np.isfinite(value) or value <= 0):
            raise ValueError(f"{name} must be finite and strictly positive")
    if not np.isfinite(inhibition) or not 0 <= inhibition <= 1:
        raise ValueError("angiogenesis_inhibition must lie in [0, 1]")
    if any(not np.isfinite(d.day) or not np.isfinite(d.amount) or d.day < 0
           or d.day >= days or d.amount < 0 for d in doses):
        raise ValueError("doses must have non-negative amount and occur before the simulation endpoint")

    n_system = len(COMBINED_STATE_NAMES)
    spatial_y0 = initial_vector(initial)
    y = np.concatenate((np.asarray(list(asdict(li).values())),
                        np.asarray([getattr(ii, name) for name in IMMUNE_DYNAMIC_STATES]),
                        spatial_y0))
    portal, arterial, _, _ = _liver_inflow(lp, flow, portal_flow, arterial_flow)

    def joint_rhs(time, state):
        system_state = np.maximum(state[:n_system], 0.0)
        spatial_state = np.maximum(state[n_system:], 0.0)
        system_dydt = _combined_rhs(
            time, system_state, lp, ip, cp, flow, inhibition,
            portal_flow, arterial_flow)
        rim_index = LIVER_STATES.index("tumour_rim")
        necrotic_index = LIVER_STATES.index("necrotic_tumour")
        vessels = system_state[LIVER_STATES.index("functional_vessels")]
        viable = (system_state[LIVER_STATES.index("tumour_rim")]
                  + system_state[LIVER_STATES.index("tumour_core")])
        perfusion = ((portal + arterial) * vessels
                     / (1.0 + lp.vessel_compression * viable))
        perfusion_factor = perfusion / b.perfusion_reference
        antigen_index = len(LIVER_STATES) + IMMUNE_DYNAMIC_STATES.index("tumour_antigen")
        occupancy_index = len(LIVER_STATES) + IMMUNE_DYNAMIC_STATES.index("pd1_occupancy")
        antigen = system_state[antigen_index]
        sources = {
            "antigen": b.antigen_drainage_gain * antigen * perfusion_factor,
            "mature_dc": b.dc_migration_gain * ip.dc_migration * antigen * perfusion_factor,
        }
        spatial_dydt = _spatial_rhs(
            time, spatial_state, sp, True, system_state[occupancy_index], sources)
        efferent = spatial_state[_spatial_index("effector", "efferent_lymph")]
        spatial_killing = b.efferent_effector_killing_gain * efferent * system_state[rim_index]
        spatial_killing = min(spatial_killing, system_state[rim_index])
        system_dydt[rim_index] -= spatial_killing
        system_dydt[necrotic_index] += spatial_killing
        return np.concatenate((system_dydt, spatial_dydt))

    events: dict[float, float] = {}
    for dose in doses:
        events[dose.day] = events.get(dose.day, 0.0) + dose.amount
    boundaries = sorted({0.0, days, *events})
    plasma_index = len(LIVER_STATES) + IMMUNE_DYNAMIC_STATES.index("anti_pd1_plasma")
    frames = []
    for start, end in zip(boundaries[:-1], boundaries[1:]):
        if start in events:
            y[plasma_index] += events[start]
        times = np.unique(np.append(np.arange(start, end, step), end))
        solved = solve_ivp(joint_rhs, (start, end), y, t_eval=times,
                           rtol=min(systemic_rtol, spatial_rtol),
                           atol=min(systemic_atol, spatial_atol))
        if not solved.success:
            raise RuntimeError(f"joint liver–TDLN integration failed: {solved.message}")
        system_values = _validated_combined_states(
            solved.y[:n_system], systemic_rtol, systemic_atol)
        spatial_values = _validated_spatial_solution(
            solved.y[n_system:], spatial_rtol, spatial_atol)
        values = np.vstack((system_values, spatial_values)).T
        if frames:
            frames[-1] = frames[-1][:-1] if start in events else frames[-1]
            if start not in events:
                values, times = values[1:], solved.t[1:]
        frames.append([(float(time), value) for time, value in zip(times, values)])
        y = values[-1]
    records = [item for frame in frames for item in frame]
    timeline = np.asarray([item[0] for item in records])
    values = np.asarray([item[1] for item in records])
    systemic = pd.DataFrame(values[:, :n_system], columns=COMBINED_STATE_NAMES)
    systemic.insert(0, "day", timeline)
    systemic.insert(1, "administered_dose", systemic.day.map(events).fillna(0.0))
    systemic = _decorate_coupled_output(
        systemic, cp, lp, flow, portal_flow, arterial_flow)

    rows = []
    spatial_values = values[:, n_system:]
    for time_index, day in enumerate(timeline):
        for species in SPECIES:
            for zone in ZONE_NAMES:
                rows.append({"day": day, "species": species, "zone": zone,
                             "amount": spatial_values[time_index,
                                                       _spatial_index(species, zone)]})
    spatial = pd.DataFrame(rows)
    antigen = systemic.tumour_antigen.to_numpy(dtype=float)
    perfusion_factor = systemic.liver_perfusion.to_numpy(dtype=float) / b.perfusion_reference
    dc_migration_flux = ip.dc_migration * antigen
    occupancy_values = systemic.pd1_occupancy.to_numpy(dtype=float)
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
    efferent = spatial[
        spatial.species.eq("effector") & spatial.zone.eq("efferent_lymph")
    ]
    systemic["spatial_efferent_effector"] = np.interp(
        systemic.day, efferent.day, efferent.amount)
    systemic["spatial_effector_killing_flux"] = np.minimum(
        b.efferent_effector_killing_gain * systemic.spatial_efferent_effector
        * systemic.liver_tumour_rim,
        systemic.liver_tumour_rim)
    spatial["spatial_effector_killing_flux"] = np.interp(
        spatial.day, systemic.day, systemic.spatial_effector_killing_flux)
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

