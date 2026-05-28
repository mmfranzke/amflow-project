from lsmethod.equations import EQUATIONS


def test_core_equation_labels_exist():
    required = [
        "delta0",
        "delta1",
        "closed_form",
        "method1_w_integral",
        "C_u",
        "method2_u_integral",
        "convolution",
    ]

    for key in required:
        assert key in EQUATIONS
        assert EQUATIONS[key].startswith("eq:")
