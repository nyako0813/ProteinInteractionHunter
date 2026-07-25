"""Shared strict model behavior."""

from pydantic import BaseModel, ConfigDict


class StrictModel(BaseModel):
    """Reject unknown fields and validate assignments."""

    model_config = ConfigDict(
        extra="forbid",
        validate_assignment=True,
        use_enum_values=False,
    )
