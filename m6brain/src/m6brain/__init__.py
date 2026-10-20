"""Zero-shot clock of brain proteome age from m6Brain, an m6A-aware model of brain proteome aging."""

from .core import gene, score, weights, weights_version

__version__ = "0.1.0"
__all__ = ["__version__", "gene", "score", "weights", "weights_version"]
