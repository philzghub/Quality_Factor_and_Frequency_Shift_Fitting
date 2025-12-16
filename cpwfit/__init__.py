"""CPW resonator fitting: Qi & frequency shift models."""
from importlib.metadata import version, PackageNotFoundError

# Public re-exports (nice UX)
from . import plotting
from .models import qint_singlecurve, qint_allcurves_lin, qint_allcurves_semilog , qint_loss_channels, freqshift_loss_channels

# Optional version string
try:
    __version__ = version("cpwfit")
except PackageNotFoundError:
    __version__ = "0.0.0"

__all__ = [
    "plotting",
    "qint_singlecurve",
    "qint_allcurves",
    "qint_loss_channels",
    "freqshift_loss_channels",
    "__version__",
]
