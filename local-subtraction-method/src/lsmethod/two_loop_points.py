from dataclasses import dataclass
from fractions import Fraction


@dataclass(frozen=True)
class TwoLoopPoint:
    name: str
    p_plus: Fraction
    p_minus: Fraction
    p_perp2: Fraction
    x: Fraction
    y: Fraction
    ml2: Fraction
    Ml2: Fraction
    mk2: Fraction
    Mk2: Fraction

    @property
    def p2(self):
        return self.p_plus * self.p_minus - self.p_perp2

    @property
    def M2(self):
        if self.ml2 == self.Ml2 == self.mk2 == self.Mk2:
            return self.ml2
        raise ValueError("Method 3 currently implements the equal-mass setup.")

    @property
    def lam(self):
        return self.y / (self.p_plus - self.x)

    def in_support(self):
        return self.x > 0 and self.y > 0 and self.x + self.y < self.p_plus


POINTS = {
    "branch_safe_rational": TwoLoopPoint(
        name="branch_safe_rational",
        p_plus=Fraction(1),
        p_minus=Fraction(0),
        p_perp2=Fraction(0),
        x=Fraction(3, 10),
        y=Fraction(1, 4),
        ml2=Fraction(49, 100),
        Ml2=Fraction(4),
        mk2=Fraction(1, 4),
        Mk2=Fraction(81, 100),
    ),
    "equal_mass_onshell_branch": TwoLoopPoint(
        name="equal_mass_onshell_branch",
        p_plus=Fraction(1),
        p_minus=Fraction(0),
        p_perp2=Fraction(0),
        x=Fraction(1, 4),
        y=Fraction(1, 4),
        ml2=Fraction(1),
        Ml2=Fraction(1),
        mk2=Fraction(1),
        Mk2=Fraction(1),
    ),
    "equal_mass_offshell_positive": TwoLoopPoint(
        name="equal_mass_offshell_positive",
        p_plus=Fraction(1),
        p_minus=Fraction(5),
        p_perp2=Fraction(0),
        x=Fraction(1, 4),
        y=Fraction(1, 4),
        ml2=Fraction(1),
        Ml2=Fraction(1),
        mk2=Fraction(1),
        Mk2=Fraction(1),
    ),
}


def get_point(name):
    try:
        return POINTS[name]
    except KeyError as exc:
        known = ", ".join(sorted(POINTS))
        raise ValueError(f"Unknown point {name!r}. Known points: {known}") from exc
