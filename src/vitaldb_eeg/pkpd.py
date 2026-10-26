"""Reconstruct propofol and remifentanil plasma and effect-site concentrations."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.linalg import expm


@dataclass(frozen=True)
class ThreeCompartmentEffectSiteModel:
    """Linear three-compartment PK model with an effect-site compartment."""

    v1: float
    v2: float
    v3: float
    cl1: float
    cl2: float
    cl3: float
    ke0: float
    amount_unit: str
    concentration_unit: str

    def system_matrix_per_minute(self) -> np.ndarray:
        values = (self.v1, self.v2, self.v3, self.cl1, self.cl2, self.cl3, self.ke0)
        if not all(np.isfinite(value) and value > 0 for value in values):
            raise ValueError(f"PK/PD parameters must be positive and finite: {values}")
        return np.asarray(
            [
                [
                    -(self.cl1 + self.cl2 + self.cl3) / self.v1,
                    self.cl2 / self.v2,
                    self.cl3 / self.v3,
                    0.0,
                ],
                [self.cl2 / self.v1, -self.cl2 / self.v2, 0.0, 0.0],
                [self.cl3 / self.v1, 0.0, -self.cl3 / self.v3, 0.0],
                [self.ke0 / self.v1, 0.0, 0.0, -self.ke0],
            ],
            dtype=float,
        )


def marsh_propofol_model(
    weight_kg: float, *, ke0_per_minute: float = 0.26
) -> ThreeCompartmentEffectSiteModel:
    """Marsh propofol model used by the original Diprifusor implementation."""

    weight = float(weight_kg)
    if not np.isfinite(weight) or weight <= 0:
        raise ValueError("weight_kg must be positive and finite")
    v1 = 0.228 * weight
    return ThreeCompartmentEffectSiteModel(
        v1=v1,
        v2=0.463 * weight,
        v3=2.893 * weight,
        # Marsh published k10, k12 and k13 as rate constants (/min).
        # Convert them to clearances by multiplying by the central volume.
        cl1=0.119 * v1,
        cl2=0.112 * v1,
        cl3=0.0419 * v1,
        ke0=float(ke0_per_minute),
        amount_unit="mg",
        concentration_unit="mg/L (= ug/mL)",
    )


def james_lean_body_mass_kg(*, sex: str, weight_kg: float, height_cm: float) -> float:
    weight = float(weight_kg)
    height = float(height_cm)
    if not np.isfinite(weight) or weight <= 0 or not np.isfinite(height) or height <= 0:
        raise ValueError("weight_kg and height_cm must be positive and finite")
    normalized = str(sex).strip().upper()
    if normalized in {"M", "MALE", "1"}:
        return float(1.10 * weight - 128.0 * (weight / height) ** 2)
    if normalized in {"F", "FEMALE", "0"}:
        return float(1.07 * weight - 148.0 * (weight / height) ** 2)
    raise ValueError(f"Unsupported sex value for the Minto model: {sex!r}")


def minto_remifentanil_model(
    *, age_years: float, sex: str, weight_kg: float, height_cm: float
) -> ThreeCompartmentEffectSiteModel:
    """Minto remifentanil population model with the published age-dependent ke0."""

    age = float(age_years)
    if not np.isfinite(age) or age <= 0:
        raise ValueError("age_years must be positive and finite")
    lbm = james_lean_body_mass_kg(sex=sex, weight_kg=weight_kg, height_cm=height_cm)
    age_delta = age - 40.0
    lbm_delta = lbm - 55.0
    return ThreeCompartmentEffectSiteModel(
        v1=5.10 - 0.0201 * age_delta + 0.0720 * lbm_delta,
        v2=9.82 - 0.0811 * age_delta + 0.1080 * lbm_delta,
        v3=5.42,
        cl1=2.60 - 0.0162 * age_delta + 0.0191 * lbm_delta,
        cl2=2.05 - 0.0301 * age_delta,
        cl3=0.076 - 0.00113 * age_delta,
        ke0=0.595 - 0.007 * age_delta,
        amount_unit="ug",
        concentration_unit="ug/L (= ng/mL)",
    )


def _transition(
    model: ThreeCompartmentEffectSiteModel, step_seconds: float
) -> tuple[np.ndarray, np.ndarray]:
    step_minutes = float(step_seconds) / 60.0
    if not np.isfinite(step_minutes) or step_minutes <= 0:
        raise ValueError("step_seconds must be positive and finite")
    augmented = np.zeros((5, 5), dtype=float)
    augmented[:4, :4] = model.system_matrix_per_minute()
    augmented[0, 4] = 1.0
    transition = expm(augmented * step_minutes)
    return transition[:4, :4], transition[:4, 4]


def effect_site_concentration_from_rate(
    infusion_rate_per_minute: np.ndarray,
    *,
    model: ThreeCompartmentEffectSiteModel,
    step_seconds: float = 1.0,
    initial_state: np.ndarray | None = None,
) -> np.ndarray:
    """Reconstruct Ce using zero-order-hold infusion rates on a regular grid.

    The value at index ``i`` is the end-of-step effect-site concentration after
    applying ``infusion_rate_per_minute[i]`` during one grid interval.
    """

    rates = np.asarray(infusion_rate_per_minute, dtype=float).reshape(-1)
    if np.any(np.isfinite(rates) & (rates < 0)):
        raise ValueError("Infusion rates cannot be negative")
    phi, gamma = _transition(model, step_seconds)
    state = (
        np.zeros(4, dtype=float)
        if initial_state is None
        else np.asarray(initial_state, dtype=float).copy()
    )
    if state.shape != (4,):
        raise ValueError("initial_state must have four elements")
    output = np.full(rates.shape, np.nan, dtype=float)
    for index, rate in enumerate(rates):
        if not np.isfinite(rate):
            continue
        state = phi @ state + gamma * float(rate)
        output[index] = state[3]
    return output
