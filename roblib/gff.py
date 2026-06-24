"""
Small GFF/GFF3 reader and writer helpers.

The writer supports Biopython SeqRecord/SeqFeature objects without depending on
BCBBio-Gff. Coordinates are written as GFF3 1-based inclusive positions.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List, TextIO, Tuple, Union
from urllib.parse import quote, unquote


@dataclass
class GffRecord:
    """
    One GFF/GFF3 feature row.
    """

    seqid: str
    source: str
    feature_type: str
    start: int
    end: int
    score: str = "."
    strand: str = "."
    phase: str = "."
    attributes: Dict[str, List[str]] = field(default_factory=dict)


def _open_text(path_or_handle: Union[str, Path, TextIO], mode: str) -> Tuple[TextIO, bool]:
    if hasattr(path_or_handle, "read") or hasattr(path_or_handle, "write"):
        return path_or_handle, False
    return open(path_or_handle, mode), True


def parse_gff_attributes(text: str) -> Dict[str, List[str]]:
    """
    Parse GFF3 key=value attributes, with a small fallback for GFF2 key value
    attributes.
    """

    attributes: Dict[str, List[str]] = {}
    text = text.strip()
    if not text or text == ".":
        return attributes

    for item in text.split(";"):
        item = item.strip()
        if not item:
            continue

        if "=" in item:
            key, value = item.split("=", 1)
        elif " " in item:
            key, value = item.split(None, 1)
            value = value.strip().strip('"')
        else:
            key, value = item, ""

        key = unquote(key.strip())
        values = [unquote(v.strip()) for v in value.split(",") if v.strip()]
        attributes[key] = values or [""]

    return attributes


def format_gff_attributes(attributes: Dict[str, Iterable[Any]]) -> str:
    """
    Format attributes as GFF3 key=value pairs.
    """

    parts = []
    for key, values in attributes.items():
        if values is None:
            continue

        if isinstance(values, (str, bytes)):
            value_list = [values]
        else:
            value_list = list(values)

        if len(value_list) == 0:
            continue

        escaped_key = quote(str(key), safe="")
        escaped_values = [
            quote(str(value), safe=f"~!$'()*+./:-@_")
            for value in value_list
            if value is not None
        ]
        if len(escaped_values) == 0:
            continue
        parts.append(f"{escaped_key}={','.join(escaped_values)}")

    return ";".join(parts) if parts else "."


def read_gff(path_or_handle: Union[str, Path, TextIO]) -> Iterator[GffRecord]:
    """
    Yield GffRecord objects from a GFF or GFF3 file.
    """

    handle, close_handle = _open_text(path_or_handle, "r")
    try:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\n")
            if line == "##FASTA":
                break
            if not line or line.startswith("#"):
                continue

            columns = line.split("\t")
            if len(columns) != 9:
                raise ValueError(f"Line {line_number} is not a valid GFF row: {line}")

            yield GffRecord(
                seqid=columns[0],
                source=columns[1],
                feature_type=columns[2],
                start=int(columns[3]),
                end=int(columns[4]),
                score=columns[5],
                strand=columns[6],
                phase=columns[7],
                attributes=parse_gff_attributes(columns[8]),
            )
    finally:
        if close_handle:
            handle.close()


read_gff3 = read_gff


def write_gff3(records: Iterable[Any], path_or_handle: Union[str, Path, TextIO]) -> None:
    """
    Write GFF3 for an iterable of GffRecord or Biopython SeqRecord objects.
    """

    handle, close_handle = _open_text(path_or_handle, "w")
    try:
        handle.write("##gff-version 3\n")
        for record in records:
            if isinstance(record, GffRecord):
                _write_gff_record(record, handle)
            else:
                _write_seqrecord(record, handle)
    finally:
        if close_handle:
            handle.close()


def _write_gff_record(record: GffRecord, handle: TextIO) -> None:
    handle.write(
        "\t".join(
            [
                record.seqid,
                record.source,
                record.feature_type,
                str(record.start),
                str(record.end),
                record.score,
                record.strand,
                record.phase,
                format_gff_attributes(record.attributes),
            ]
        )
        + "\n"
    )


def _write_seqrecord(record: Any, handle: TextIO) -> None:
    seqid = str(getattr(record, "id", "."))
    seq_length = len(getattr(record, "seq", ""))
    if seq_length > 0:
        handle.write(f"##sequence-region {seqid} 1 {seq_length}\n")

    features = getattr(record, "features", [])
    for index, feature in enumerate(features, start=1):
        for gff_record in _feature_to_gff_records(record, feature, index):
            _write_gff_record(gff_record, handle)


def _feature_to_gff_records(record: Any, feature: Any, index: int) -> Iterator[GffRecord]:
    location = getattr(feature, "location", None)
    if location is None:
        return

    qualifiers = getattr(feature, "qualifiers", {})
    source = _first_qualifier(qualifiers, "source") or "GenBank"
    feature_type = str(getattr(feature, "type", "feature"))
    attributes = _feature_attributes(record, feature, index)
    phase = _feature_phase(feature)

    parent_id = attributes.get("ID", [None])[0]
    parts = list(getattr(location, "parts", []) or [])
    if len(parts) <= 1:
        yield GffRecord(
            seqid=str(getattr(record, "id", ".")),
            source=source,
            feature_type=feature_type,
            start=_location_start(location),
            end=_location_end(location),
            strand=_location_strand(location),
            phase=phase,
            attributes=attributes,
        )
        return

    yield GffRecord(
        seqid=str(getattr(record, "id", ".")),
        source=source,
        feature_type=feature_type,
        start=_location_start(location),
        end=_location_end(location),
        strand=_location_strand(location),
        phase=phase,
        attributes=attributes,
    )

    for part_index, part in enumerate(parts, start=1):
        part_attributes: Dict[str, List[str]] = {"Parent": [parent_id]} if parent_id else {}
        if parent_id:
            part_attributes["ID"] = [f"{parent_id}.part{part_index}"]

        yield GffRecord(
            seqid=str(getattr(record, "id", ".")),
            source=source,
            feature_type=f"{feature_type}_part",
            start=_location_start(part),
            end=_location_end(part),
            strand=_location_strand(part),
            phase=phase,
            attributes=part_attributes,
        )


def _feature_attributes(record: Any, feature: Any, index: int) -> Dict[str, List[str]]:
    qualifiers = getattr(feature, "qualifiers", {})
    feature_type = str(getattr(feature, "type", "feature"))
    feature_id = _first_qualifier(qualifiers, "ID")
    if feature_id is None:
        feature_id = (
            _first_qualifier(qualifiers, "protein_id")
            or _first_qualifier(qualifiers, "locus_tag")
            or _first_qualifier(qualifiers, "gene")
            or f"{getattr(record, 'id', 'record')}.{feature_type}.{index}"
        )

    attributes: Dict[str, List[str]] = {"ID": [feature_id]}

    name = (
        _first_qualifier(qualifiers, "Name")
        or _first_qualifier(qualifiers, "gene")
        or _first_qualifier(qualifiers, "locus_tag")
        or _first_qualifier(qualifiers, "product")
    )
    if name:
        attributes["Name"] = [name]

    for key, values in qualifiers.items():
        if key in {"ID", "Name"}:
            continue
        if isinstance(values, (str, bytes)):
            attributes[key] = [str(values)]
        else:
            attributes[key] = [str(value) for value in values]

    return attributes


def _first_qualifier(qualifiers: Dict[str, Any], key: str) -> Union[str, None]:
    values = qualifiers.get(key)
    if values is None:
        return None
    if isinstance(values, (str, bytes)):
        return str(values)
    if len(values) == 0:
        return None
    return str(values[0])


def _location_start(location: Any) -> int:
    return int(getattr(location, "start")) + 1


def _location_end(location: Any) -> int:
    return int(getattr(location, "end"))


def _location_strand(location: Any) -> str:
    strand = getattr(location, "strand", None)
    if strand == 1:
        return "+"
    if strand == -1:
        return "-"
    return "."


def _feature_phase(feature: Any) -> str:
    if getattr(feature, "type", None) != "CDS":
        return "."

    codon_start = _first_qualifier(getattr(feature, "qualifiers", {}), "codon_start")
    if codon_start is None:
        return "0"

    try:
        return str((int(codon_start) - 1) % 3)
    except ValueError:
        return "0"
