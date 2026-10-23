"""m6Brain clock of brain proteome age: an elastic net guided by m6Brain, an m6A-aware model of brain proteome aging."""

from .core import gene, predict_age, score, weights, weights_version

__version__ = "0.2.0"
__all__ = ["__version__", "gene", "predict_age", "score", "weights", "weights_version"]
