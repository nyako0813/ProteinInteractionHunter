"""Typed errors used by local validation and serialization."""


class ProteinInteractionHunterError(Exception):
    """Base application error."""


class ConfigurationError(ProteinInteractionHunterError):
    """Configuration loading or validation failed."""


class InputValidationError(ProteinInteractionHunterError):
    """A required local input could not be validated."""


class SerializationError(ProteinInteractionHunterError):
    """A validated model could not be serialized safely."""
