"""Project configuration loading with a simple YAML fallback."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from lanm.paths import CONFIG_PATH

try:
    import yaml  # type: ignore
except ModuleNotFoundError:  # pragma: no cover - exercised when pyyaml is absent
    yaml = None


@dataclass(frozen=True, slots=True)
class ProjectConfig:
    project_name: str
    target_metal: str
    competitors: tuple[str, ...]
    temperature_K: int
    first_shell_cutoff_A: float
    second_sphere_cutoff_A: float
    ddg_target_kcal_mol: float
    ore_weights_status: str


def _coerce_scalar(value: str) -> object:
    stripped = value.strip()
    if stripped.startswith("[") and stripped.endswith("]"):
        items = [item.strip() for item in stripped[1:-1].split(",") if item.strip()]
        return [item.strip("'\"") for item in items]
    if stripped.lower() in {"true", "false"}:
        return stripped.lower() == "true"
    try:
        if "." in stripped or "e" in stripped.lower():
            return float(stripped)
        return int(stripped)
    except ValueError:
        return stripped.strip("'\"")


def _load_yaml_with_fallback(path: Path) -> dict[str, object]:
    if yaml is not None:
        payload = yaml.safe_load(path.read_text(encoding="utf-8"))
        if isinstance(payload, dict):
            return payload
        raise ValueError(f"Unexpected YAML payload in {path}")
    parsed: dict[str, object] = {}
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        key, value = line.split(":", 1)
        parsed[key.strip()] = _coerce_scalar(value)
    return parsed


def load_project_config(path: Path = CONFIG_PATH) -> ProjectConfig:
    payload = _load_yaml_with_fallback(path)
    return ProjectConfig(
        project_name=str(payload["project_name"]),
        target_metal=str(payload["target_metal"]),
        competitors=tuple(str(item) for item in payload["competitors"]),
        temperature_K=int(payload["temperature_K"]),
        first_shell_cutoff_A=float(payload["first_shell_cutoff_A"]),
        second_sphere_cutoff_A=float(payload["second_sphere_cutoff_A"]),
        ddg_target_kcal_mol=float(payload["ddg_target_kcal_mol"]),
        ore_weights_status=str(payload["ore_weights_status"]),
    )
