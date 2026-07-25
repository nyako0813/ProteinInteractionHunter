"""TSV writer for the flat manual structure-review queue."""

import csv
import io
from collections.abc import Sequence
from pathlib import Path

from protein_interaction_hunter.models.structure_queue import StructurePredictionQueueEntry

STRUCTURE_QUEUE_COLUMNS = tuple(StructurePredictionQueueEntry.model_fields)


class StructureQueueTsvWriter:
    def write(self, entries: Sequence[StructurePredictionQueueEntry], path: Path) -> Path:
        output_path = path.expanduser().resolve()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        buffer = io.StringIO(newline="")
        writer = csv.DictWriter(
            buffer,
            fieldnames=STRUCTURE_QUEUE_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for entry in sorted(
            entries, key=lambda item: (item.rank, item.query_id, item.candidate_id)
        ):
            writer.writerow(entry.model_dump(mode="json"))
        output_path.write_text(buffer.getvalue(), encoding="utf-8", newline="\n")
        return output_path
