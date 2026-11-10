"""Metadata for reproducible scientific figures."""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

FIGURE_METADATA_SCHEMA_VERSION = "1.0"


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


@dataclass(frozen=True)
class FigureContext:
    source_structure: str
    structure_sha256: str
    protein_chain: str
    ligand_instance: str
    model_number: int
    coordinate_units: str = "angstrom"
    generated_from: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return {
            "source_structure": self.source_structure,
            "structure_sha256": self.structure_sha256,
            "protein_chain": self.protein_chain,
            "ligand_instance": self.ligand_instance,
            "model_number": self.model_number,
            "coordinate_units": self.coordinate_units,
            "generated_from": list(self.generated_from),
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
            generated_from=_string_tuple(data, "generated_from"),
        )


@dataclass(frozen=True)
class FigureMetadata:
    figure_id: str
    title: str
    image_path: str
    context: FigureContext
    visualization_tool: str
    tool_version: str
    scientific_limitations: tuple[str, ...]
    alt_text: str
    caption: str
    camera_view: tuple[float, ...] | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "figure_id": self.figure_id,
            "title": self.title,
            "image_path": self.image_path,
            "context": self.context.to_dict(),
            "visualization_tool": self.visualization_tool,
            "tool_version": self.tool_version,
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
        if not figure.scientific_limitations:
            errors.append(f"figures[{index}] requires scientific limitations.")
        if not figure.alt_text:
            errors.append(f"figures[{index}] requires alt text.")
    return errors
