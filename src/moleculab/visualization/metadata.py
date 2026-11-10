"""Metadata for reproducible scientific figures."""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

FIGURE_METADATA_SCHEMA_VERSION = "1.1"


def _as_mapping(value: Any, name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be a mapping.")
    return value


def _as_list(value: Any, name: str) -> list[Any]:
    if value is None:
        return []
    if not isinstance(value, (list, tuple)):
        raise ValueError(f"{name} must be a list.")
    return list(value)


def _required_string(data: Mapping[str, Any], key: str) -> str:
    value = data.get(key)
    if value is None or str(value) == "":
        raise ValueError(f"{key} is required.")
    return str(value)


def _required_int(data: Mapping[str, Any], key: str) -> int:
    value = data.get(key)
    if value is None:
        raise ValueError(f"{key} is required.")
    return int(value)


def _string_tuple(data: Mapping[str, Any], key: str) -> tuple[str, ...]:
    return tuple(str(item) for item in _as_list(data.get(key), key))


def _float_tuple(data: Mapping[str, Any], key: str) -> tuple[float, ...] | None:
    value = data.get(key)
    if value is None:
        return None
    return tuple(float(item) for item in _as_list(value, key))


def _pair_tuple(
    data: Mapping[str, Any],
    key: str,
) -> tuple[tuple[str, str], ...]:
    value = data.get(key)
    if value is None:
        return ()
    if isinstance(value, Mapping):
        return tuple(sorted((str(k), str(v)) for k, v in value.items()))
    pairs = []
    for item in _as_list(value, key):
        if isinstance(item, (list, tuple)) and len(item) == 2:
            pairs.append((str(item[0]), str(item[1])))
        else:
            raise ValueError(f"{key} must contain two-item pairs.")
    return tuple(sorted(pairs))


def rendering_config_sha256(
    config: Sequence[tuple[str, str]],
) -> str:
    """Return a stable fingerprint for a rendering configuration."""

    payload = json.dumps(
        {key: value for key, value in sorted(config)},
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


@dataclass(frozen=True)
class FigureContext:
    source_structure: str
    structure_sha256: str
    protein_chain: str
    ligand_instance: str
    model_number: int
    coordinate_units: str = "angstrom"

    def to_dict(self) -> dict[str, Any]:
        return {
            "source_structure": self.source_structure,
            "structure_sha256": self.structure_sha256,
            "protein_chain": self.protein_chain,
            "ligand_instance": self.ligand_instance,
            "model_number": self.model_number,
            "coordinate_units": self.coordinate_units,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> FigureContext:
        data = _as_mapping(data, "figure context")
        return cls(
            source_structure=_required_string(data, "source_structure"),
            structure_sha256=_required_string(data, "structure_sha256"),
            protein_chain=_required_string(data, "protein_chain"),
            ligand_instance=_required_string(data, "ligand_instance"),
            model_number=_required_int(data, "model_number"),
            coordinate_units=str(data.get("coordinate_units", "angstrom")),
        )


@dataclass(frozen=True)
class FigureMetadata:
    figure_id: str
    title: str
    image_path: str
    context: FigureContext
    visualization_tool: str
    tool_version: str
    source_records: tuple[str, ...]
    source_checksums: tuple[tuple[str, str], ...]
    rendering_config: tuple[tuple[str, str], ...]
    rendering_config_sha256: str | None
    scientific_limitations: tuple[str, ...]
    alt_text: str
    caption: str
    camera_view: tuple[float, ...] | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "source_records", tuple(sorted(self.source_records)))
        object.__setattr__(
            self,
            "source_checksums",
            tuple(sorted(self.source_checksums)),
        )
        object.__setattr__(
            self,
            "rendering_config",
            tuple(sorted(self.rendering_config)),
        )
        if self.rendering_config_sha256 is None:
            object.__setattr__(
                self,
                "rendering_config_sha256",
                rendering_config_sha256(self.rendering_config),
            )

    def to_dict(self) -> dict[str, Any]:
        return {
            "figure_id": self.figure_id,
            "title": self.title,
            "image_path": self.image_path,
            "context": self.context.to_dict(),
            "visualization_tool": self.visualization_tool,
            "tool_version": self.tool_version,
            "source_records": list(self.source_records),
            "source_checksums": {
                key: value for key, value in self.source_checksums
            },
            "rendering_config": {
                key: value for key, value in self.rendering_config
            },
            "rendering_config_sha256": self.rendering_config_sha256,
            "scientific_limitations": list(self.scientific_limitations),
            "alt_text": self.alt_text,
            "caption": self.caption,
            "camera_view": (
                None if self.camera_view is None else list(self.camera_view)
            ),
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> FigureMetadata:
        data = _as_mapping(data, "figure metadata")
        return cls(
            figure_id=_required_string(data, "figure_id"),
            title=_required_string(data, "title"),
            image_path=_required_string(data, "image_path"),
            context=FigureContext.from_dict(
                _as_mapping(data.get("context"), "context")
            ),
            visualization_tool=_required_string(data, "visualization_tool"),
            tool_version=_required_string(data, "tool_version"),
            source_records=_string_tuple(data, "source_records"),
            source_checksums=_pair_tuple(data, "source_checksums"),
            rendering_config=_pair_tuple(data, "rendering_config"),
            rendering_config_sha256=(
                None
                if data.get("rendering_config_sha256") is None
                else str(data.get("rendering_config_sha256"))
            ),
            scientific_limitations=_string_tuple(data, "scientific_limitations"),
            alt_text=_required_string(data, "alt_text"),
            caption=_required_string(data, "caption"),
            camera_view=_float_tuple(data, "camera_view"),
        )


def write_figure_metadata(
    path: str | Path,
    figures: Sequence[FigureMetadata],
) -> Path:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(
            {
                "schema_version": FIGURE_METADATA_SCHEMA_VERSION,
                "figures": [figure.to_dict() for figure in figures],
            },
            indent=2,
            sort_keys=False,
        )
        + "\n",
        encoding="utf-8",
    )
    return destination


def load_figure_metadata(path: str | Path) -> tuple[FigureMetadata, ...]:
    data = _as_mapping(
        json.loads(Path(path).read_text(encoding="utf-8")),
        "figure metadata collection",
    )
    return tuple(
        FigureMetadata.from_dict(_as_mapping(item, "figure metadata"))
        for item in _as_list(data.get("figures"), "figures")
    )


def validate_figure_metadata(figures: Sequence[FigureMetadata]) -> list[str]:
    errors: list[str] = []
    figure_ids = set()
    for index, figure in enumerate(figures):
        if not figure.figure_id or not figure.title or not figure.image_path:
            errors.append(
                f"figures[{index}] requires figure_id, title, and image_path."
            )
        if figure.figure_id in figure_ids:
            errors.append(f"figures[{index}].figure_id is duplicated.")
        figure_ids.add(figure.figure_id)
        if not figure.visualization_tool or not figure.tool_version:
            errors.append(
                f"figures[{index}] requires visualization_tool and tool_version."
            )
        if not figure.context.structure_sha256:
            errors.append(f"figures[{index}] requires a structure checksum.")
        if not figure.source_records:
            errors.append(f"figures[{index}] requires source records.")
        record_names = set(figure.source_records)
        checksum_names = {key for key, _ in figure.source_checksums}
        if record_names != checksum_names:
            errors.append(
                f"figures[{index}] source_checksums must cover every source record."
            )
        if not figure.rendering_config:
            errors.append(f"figures[{index}] requires a rendering configuration.")
        expected_fingerprint = rendering_config_sha256(figure.rendering_config)
        if figure.rendering_config_sha256 != expected_fingerprint:
            errors.append(
                f"figures[{index}] rendering_config_sha256 does not match the config."
            )
        if not figure.scientific_limitations:
            errors.append(f"figures[{index}] requires scientific limitations.")
        if not figure.alt_text:
            errors.append(f"figures[{index}] requires alt text.")
    return errors
