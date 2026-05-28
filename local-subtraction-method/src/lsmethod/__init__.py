__all__ = [
    "Kinematics",
    "closed_form",
    "delta0",
    "delta1",
    "method1_w_integral",
    "method2_u_integral",
]


def __getattr__(name):
    if name == "Kinematics":
        from .kinematics import Kinematics

        return Kinematics

    if name in {"closed_form", "delta0", "delta1", "method1_w_integral", "method2_u_integral"}:
        from . import kernels

        return getattr(kernels, name)

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
