from dataclasses import dataclass
import mpmath as mp
import random


# One numerical phase-space point for the double-collinear kernel.
@dataclass(frozen=True)
class Kinematics:
    x: mp.mpf
    y: mp.mpf
    p2: mp.mpf
    ml: mp.mpf
    Ml: mp.mpf
    mk: mp.mpf
    Mk: mp.mpf

    @property
    def X(self):
        """LaTeX: eq:X-kappa-def."""
        return 1 - self.x

    @property
    def lam(self):
        """LaTeX: eq:lambda."""
        return self.y / self.X

    @property
    def kappa(self):
        """LaTeX: eq:X-kappa-def."""
        return self.lam * (1 - self.lam)

    def in_support(self):
        """LaTeX: eq:support-theta."""
        # Support of the longitudinal fractions.
        return (0 <= self.x <= 1) and (0 <= self.y <= 1 - self.x)

    def mu_l_x(self):
        """LaTeX: eq:mu-l-def."""
        # Effective mass scale of the outer l-collinear pair at fixed x.
        return (1 - self.x) * self.ml**2 + self.x * self.Ml**2

    def mu_l_u(self, u):
        """LaTeX: eq:mu-l-def."""
        # Same outer-pair mass combination, evaluated at Feynman parameter u.
        return u * self.ml**2 + (1 - u) * self.Ml**2

    def mu_k_lam(self):
        """LaTeX: eq:lambda."""
        # Effective mass scale of the nested k-subloop pair.
        return (1 - self.lam) * self.mk**2 + self.lam * self.Mk**2
