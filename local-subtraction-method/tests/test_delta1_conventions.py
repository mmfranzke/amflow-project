import mpmath as mp

from lsmethod.kinematics import Kinematics
from lsmethod.kernels import delta1


def delta1_current_pdf(kin):
    lam = kin.lam
    return (
        -kin.mu_k_lam()
        + lam * (1 - lam) * kin.Ml**2
        + lam * (1 - lam * kin.x) * kin.p2
    )


def delta1_old_compare(kin):
    lam = kin.lam
    return (
        -kin.mu_k_lam()
        + lam * (1 - lam) * kin.Ml**2
        + lam * (1 - kin.x) * kin.p2
    )


def test_current_pdf_delta1_agrees_with_kernel():
    kin = Kinematics(
        x=mp.mpf("0.25"),
        y=mp.mpf("0.25"),
        p2=mp.mpf("5.0"),
        ml=mp.mpf("1.0"),
        Ml=mp.mpf("1.0"),
        mk=mp.mpf("1.0"),
        Mk=mp.mpf("1.0"),
    )

    assert abs(delta1(kin) - delta1_current_pdf(kin)) < mp.mpf("1e-40")


def test_current_pdf_and_old_compare_differ_offshell():
    kin = Kinematics(
        x=mp.mpf("0.25"),
        y=mp.mpf("0.25"),
        p2=mp.mpf("5.0"),
        ml=mp.mpf("1.0"),
        Ml=mp.mpf("1.0"),
        mk=mp.mpf("1.0"),
        Mk=mp.mpf("1.0"),
    )

    assert abs(delta1_current_pdf(kin) - delta1_old_compare(kin)) > mp.mpf("1e-20")


def test_current_pdf_and_old_compare_agree_onshell():
    kin = Kinematics(
        x=mp.mpf("0.25"),
        y=mp.mpf("0.25"),
        p2=mp.mpf("0.0"),
        ml=mp.mpf("1.0"),
        Ml=mp.mpf("1.0"),
        mk=mp.mpf("1.0"),
        Mk=mp.mpf("1.0"),
    )

    assert abs(delta1_current_pdf(kin) - delta1_old_compare(kin)) < mp.mpf("1e-40")
