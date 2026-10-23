"""Shared model parameters."""

from dataclasses import dataclass


@dataclass(frozen=True)
class Parameters:
    """Physical and numerical settings."""

    tissue_resolution: int = 12
    vessel_resolution: int = 8
    viscosity: float = 0.0035
    tissue_mobility: float = 2.0e-9
    drainage: float = 5.0e-6
    wall_conductivity: float = 5.0e-7
    terminal_conductance: float = 2.0e-11
    healthy_terminal_factor: float = 1.0
    terminal_radius: float = 0.0025
    spread_radius_factor: float = 5.0
    inlet_pressure: float = 80.0 * 133.322
    venous_pressure: float = 10.0 * 133.322
    conductance_scale: float = 2.0e-11
    angles: int = 12
    disk_rings: int = 3
    terminal_samples: int = 256
    multiplier_spacing: float = 0.004


MODELS = ("poiseuille", "darcy", "lateral", "terminal", "contact", "spread", "line_source")

MODEL_LABELS = {
    "poiseuille": "Vessel flow",
    "darcy": "Fixed permeability",
    "lateral": "Wall exchange",
    "terminal": "Terminal exchange",
    "contact": "Equal pressures",
    "spread": "Spread exchange",
    "line_source": "Line source",
}

MODEL_EQUATIONS = {
    "poiseuille": r"$q=-\frac{\pi r^4}{8\mu}P'$",
    "darcy": r"$q=-\frac{\pi r^2K_v}{\mu}P'$",
    "lateral": r"$j=2\pi r L_p(P-\overline p)$",
    "terminal": r"$Q_a=G_a(P_a-\overline p_a)$",
    "contact": r"$P=\overline p$",
    "spread": r"$j=\beta_\epsilon(P-\overline p_\epsilon)$",
    "line_source": r"$-\nabla\cdot(M\nabla p)+dp=j\delta_\Lambda$",
}
