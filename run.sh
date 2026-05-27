#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if command -v WolframKernel >/dev/null 2>&1; then
  KERNEL="$(command -v WolframKernel)"
elif command -v MathKernel >/dev/null 2>&1; then
  KERNEL="$(command -v MathKernel)"
elif command -v wolfram >/dev/null 2>&1; then
  KERNEL="$(command -v wolfram)"
elif command -v math >/dev/null 2>&1; then
  KERNEL="$(command -v math)"
elif [ -x "/Applications/Mathematica.app/Contents/MacOS/MathKernel" ]; then
  KERNEL="/Applications/Mathematica.app/Contents/MacOS/MathKernel"
elif [ -x "/Applications/Wolfram.app/Contents/MacOS/WolframKernel" ]; then
  KERNEL="/Applications/Wolfram.app/Contents/MacOS/WolframKernel"
else
  echo "Could not find Mathematica/Wolfram kernel."
  echo "On Euler, load the Mathematica module first, e.g. module load <mathematica-module>."
  exit 1
fi

TARGET="${1:-bubble}"

case "$TARGET" in
  bubble)
    # Sanity check for AMFlow, reducer, and dynamic-library setup.
    SCRIPT="$PROJECT_DIR/targets/RunBubble.wl"
    ;;
  oneloop-kernel-direct-cut)
    # Diagnostic direct cut-linear attempt. Not the recommended route.
    SCRIPT="$PROJECT_DIR/targets/RunOneLoopKernelDirectCut.wl"
    ;;
  oneloop-kernel-uncut-plus)
    # Computes F+ with the uncut L+i0 linear denominator.
    SCRIPT="$PROJECT_DIR/targets/RunOneLoopKernelUncutPlus.wl"
    ;;
  oneloop-kernel-uncut-minus)
    # Computes J- with the uncut -L+i0 linear denominator.
    SCRIPT="$PROJECT_DIR/targets/RunOneLoopKernelUncutMinus.wl"
    ;;
  compare-oneloop)
    # Fresh one-loop check: computes F+, J-, and compares directly.
    SCRIPT="$PROJECT_DIR/targets/CompareOneLoop.wl"
    ;;
  compare-oneloop-direct-cut)
    # Legacy comparison for the direct cut-linear attempt.
    SCRIPT="$PROJECT_DIR/checks/CompareOneLoopDirectCut.wl"
    ;;
  compare-oneloop-from-files)
    # Post-processes exported F+ and J- results to compare to the analytic kernel.
    SCRIPT="$PROJECT_DIR/checks/CompareOneLoopFromFiles.wl"
    ;;
  oneloop-kernel-aux-plus)
    # Computes F+ with an extra auxiliary denominator set to exponent 0.
    SCRIPT="$PROJECT_DIR/targets/RunOneLoopKernelAuxDenPlus.wl"
    ;;
  oneloop-kernel-aux-minus)
    # Computes J- with an extra auxiliary denominator set to exponent 0.
    SCRIPT="$PROJECT_DIR/targets/RunOneLoopKernelAuxDenMinus.wl"
    ;;
  compare-oneloop-aux)
    # Fresh zero-power denominator test: computes F+, J-, and compares directly.
    SCRIPT="$PROJECT_DIR/targets/CompareOneLoopAuxDen.wl"
    ;;
  compare-oneloop-aux-from-files)
    # Post-processes exported auxiliary F+ and J- results to compare to the analytic kernel.
    SCRIPT="$PROJECT_DIR/checks/CompareOneLoopAuxDenFromFiles.wl"
    ;;
  twoloop-kernel-uncut-pp)
    # Computes the PP two-loop uncut GaugeLink piece.
    SCRIPT="$PROJECT_DIR/targets/RunTwoLoopKernelUncutPP.wl"
    ;;
  
  twoloop-kernel-uncut-pm)
    # Computes the PM two-loop uncut GaugeLink piece.
    SCRIPT="$PROJECT_DIR/targets/RunTwoLoopKernelUncutPM.wl"
    ;;
  twoloop-kernel-uncut-mp)
    # Computes the MP two-loop uncut GaugeLink piece.
    SCRIPT="$PROJECT_DIR/targets/RunTwoLoopKernelUncutMP.wl"
    ;;
  twoloop-kernel-uncut-mm)
    # Computes the MM two-loop uncut GaugeLink piece.
    SCRIPT="$PROJECT_DIR/targets/RunTwoLoopKernelUncutMM.wl"
    ;;
  compare-twoloop)
    # Fresh two-loop check: computes all four pieces and compares directly.
    SCRIPT="$PROJECT_DIR/targets/CompareTwoLoop.wl"
    ;;
  compare-twoloop-quick)
    # Cheap full smoke test: all four pieces, same eps order, low precision.
    export AMFLOW_PRECISION_GOAL="${AMFLOW_PRECISION_GOAL:-6}"
    export AMFLOW_EPS_ORDER="${AMFLOW_EPS_ORDER:-4}"
    SCRIPT="$PROJECT_DIR/targets/CompareTwoLoop.wl"
    ;;
  compare-twoloop-from-files)
    # Post-processes exported two-loop sign pieces to compare to the closed form.
    SCRIPT="$PROJECT_DIR/checks/CompareTwoLoopFromFiles.wl"
    ;;
  compare-ksubloop)
    # Cheap full smoke test: all four pieces, same eps order, low precision.
    export AMFLOW_PRECISION_GOAL="${AMFLOW_PRECISION_GOAL:-6}"
    export AMFLOW_EPS_ORDER="${AMFLOW_EPS_ORDER:-4}"
    SCRIPT="$PROJECT_DIR/targets/CompareKSubloop.wl"
    ;;
  outer-eq29-plus)
    # Computes the +L GaugeLink piece for the direct Eq. (29) outer integral.
    export AMFLOW_PRECISION_GOAL="${AMFLOW_PRECISION_GOAL:-6}"
    export AMFLOW_EPS_ORDER="${AMFLOW_EPS_ORDER:-4}"
    SCRIPT="$PROJECT_DIR/targets/RunOuterEq29Plus.wl"
    ;;
  outer-eq29-minus)
    # Computes the -L GaugeLink piece for the direct Eq. (29) outer integral.
    export AMFLOW_PRECISION_GOAL="${AMFLOW_PRECISION_GOAL:-6}"
    export AMFLOW_EPS_ORDER="${AMFLOW_EPS_ORDER:-4}"
    SCRIPT="$PROJECT_DIR/targets/RunOuterEq29Minus.wl"
    ;;
  compare-outer-eq29)
    # Fresh direct Eq. (29) check: computes +L, -L, and compares to Eq. (51)/(53).
    export AMFLOW_PRECISION_GOAL="${AMFLOW_PRECISION_GOAL:-6}"
    export AMFLOW_EPS_ORDER="${AMFLOW_EPS_ORDER:-4}"
    SCRIPT="$PROJECT_DIR/targets/CompareOuterEq29.wl"
    ;;
  *)
    echo "Unknown target: $TARGET"
    echo "Allowed: bubble, oneloop-kernel-direct-cut, oneloop-kernel-uncut-plus, oneloop-kernel-uncut-minus, compare-oneloop, compare-oneloop-direct-cut, compare-oneloop-from-files, oneloop-kernel-aux-plus, oneloop-kernel-aux-minus, compare-oneloop-aux, compare-oneloop-aux-from-files, twoloop-kernel-uncut-pp, twoloop-kernel-uncut-pm, twoloop-kernel-uncut-mp, twoloop-kernel-uncut-mm, compare-twoloop, compare-twoloop-quick, compare-twoloop-from-files, compare-ksubloop, outer-eq29-plus, outer-eq29-minus, compare-outer-eq29"
    exit 1
    ;;
esac

mkdir -p "$PROJECT_DIR/logs"

echo "Running $SCRIPT"
"$KERNEL" -script "$SCRIPT" 2>&1 | tee "$PROJECT_DIR/logs/${TARGET}.log"
