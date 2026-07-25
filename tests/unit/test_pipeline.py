"""Unimplemented analyses must fail explicitly."""

import pytest

from protein_interaction_hunter.application.pipeline import InteractionCandidatePipeline


def test_pipeline_does_not_fake_mvp0_success() -> None:
    with pytest.raises(NotImplementedError, match="not implemented in MVP-0"):
        InteractionCandidatePipeline().run()
