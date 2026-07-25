"""Side-effect-free logging configuration."""

import logging
from pathlib import Path


def configure_logging(level: str, directory: Path | None = None) -> logging.Logger:
    """Create an application logger only when explicitly called."""
    logger = logging.getLogger("protein_interaction_hunter")
    logger.setLevel(level)
    logger.handlers.clear()
    logger.propagate = False
    stream = logging.StreamHandler()
    stream.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
    logger.addHandler(stream)
    if directory is not None:
        directory.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(directory / "run.log", encoding="utf-8")
        file_handler.setFormatter(
            logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")
        )
        logger.addHandler(file_handler)
    return logger
