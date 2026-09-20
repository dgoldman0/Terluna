"""Uniform-motion horizon diagnostic, not a lunar twilight radiative model."""
from __future__ import annotations
import json

def diagnostic(period_days: float = 29.53, interval_deg: float = 6.0, solar_diameter_arcmin: float = 32.0) -> dict[str, float]:
    if period_days <= 0 or interval_deg < 0 or solar_diameter_arcmin <= 0:
        raise ValueError("Period and solar diameter must be positive; angular interval must be nonnegative.")
    rate = 360.0 / (period_days * 24.0)
    return {"period_days": period_days, "angular_rate_deg_per_hour": rate,
            "interval_degrees": interval_deg, "angular_interval_hours": interval_deg / rate,
            "solar_diameter_arcmin": solar_diameter_arcmin,
            "solar_disk_crossing_hours": (solar_diameter_arcmin / 60.0) / rate}

if __name__ == "__main__":
    print(json.dumps(diagnostic(), indent=2))
