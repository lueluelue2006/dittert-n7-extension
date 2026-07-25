from __future__ import annotations

import contextlib
import io
import subprocess
import sys
import unittest
from fractions import Fraction as F
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import dittert_n6_certificates as n6  # noqa: E402
import dittert_n7_n15_certificates as n7  # noqa: E402


def bernstein_basis(degree: int, index: int, t: F) -> F:
    return F(comb(degree, index)) * t**index * (1 - t) ** (degree - index)


def eval_bivariate_grid(grid: list[list[F]], u: F, v: F) -> F:
    dx = len(grid) - 1
    dy = len(grid[0]) - 1
    return sum(
        grid[i][j]
        * bernstein_basis(dx, i, u)
        * bernstein_basis(dy, j, v)
        for i in range(dx + 1)
        for j in range(dy + 1)
    )


def eval_sparse_bivariate(poly: n7.Poly, x: F, y: F) -> F:
    return sum(c * x**i * y**j for (i, j), c in poly.items())


def eval_tensor(arr, parameters: tuple[F, ...]) -> F:
    shape = arr.shape
    total = F(0)
    import itertools

    for index in itertools.product(*[range(size) for size in shape]):
        weight = F(1)
        for axis, i in enumerate(index):
            degree = shape[axis] - 1
            weight *= bernstein_basis(degree, i, parameters[axis])
        total += arr[index] * weight
    return total


class ImportAndStateTests(unittest.TestCase):
    def test_n7_import_has_no_output(self) -> None:
        code = "import dittert_n7_n15_certificates"
        completed = subprocess.run(
            [sys.executable, "-c", code],
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=True,
        )
        self.assertEqual(completed.stdout, "")
        self.assertEqual(completed.stderr, "")

    def test_n7_gamma_state_is_immutable(self) -> None:
        before = n7.single_bridge_polynomials(1, 5)
        n7.verify_n15()
        after = n7.single_bridge_polynomials(1, 5)
        self.assertEqual(before, after)
        self.assertEqual(n7.GAMMA7, n7.gamma(7))

    def test_n6_import_has_no_verification_output(self) -> None:
        code = "import dittert_n6_certificates"
        completed = subprocess.run(
            [sys.executable, "-c", code],
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=True,
        )
        self.assertEqual(completed.stdout, "")


class BivariateBernsteinTests(unittest.TestCase):
    def setUp(self) -> None:
        self.poly = {
            (0, 0): F(1),
            (1, 0): F(2),
            (0, 1): F(3),
            (1, 1): F(4),
            (2, 0): F(5),
        }
        self.x0, self.x1 = F(1, 5), F(3, 5)
        self.y0, self.y1 = F(1, 7), F(4, 7)
        self.grid = n7.power_to_bernstein_on_box(
            self.poly, self.x0, self.x1, self.y0, self.y1
        )

    def test_conversion_preserves_exact_values(self) -> None:
        u, v = F(2, 5), F(3, 5)
        x = self.x0 + (self.x1 - self.x0) * u
        y = self.y0 + (self.y1 - self.y0) * v
        self.assertEqual(
            eval_bivariate_grid(self.grid, u, v),
            eval_sparse_bivariate(self.poly, x, y),
        )

    def test_de_casteljau_children_cover_parent(self) -> None:
        u, v = F(2, 7), F(4, 9)
        left, right = n7.split_x(self.grid)
        self.assertEqual(
            eval_bivariate_grid(left, u, v),
            eval_bivariate_grid(self.grid, u / 2, v),
        )
        self.assertEqual(
            eval_bivariate_grid(right, u, v),
            eval_bivariate_grid(self.grid, (1 + u) / 2, v),
        )
        low, high = n7.split_y(self.grid)
        self.assertEqual(
            eval_bivariate_grid(low, u, v),
            eval_bivariate_grid(self.grid, u, v / 2),
        )
        self.assertEqual(
            eval_bivariate_grid(high, u, v),
            eval_bivariate_grid(self.grid, u, (1 + v) / 2),
        )

    def test_constant_and_degenerate_degree(self) -> None:
        grid = n7.power_to_bernstein_on_box(
            n7.pconst(7), F(0), F(1), F(0), F(1)
        )
        self.assertEqual(grid, [[F(7)]])


class DimensionSixBernsteinTests(unittest.TestCase):
    def test_univariate_conversion_and_split(self) -> None:
        t = n6.sp.symbols("t")
        poly = n6.sp.Poly(2 + 3 * t + 5 * t**2, t, domain=n6.sp.QQ)
        lo, hi = n6.sp.Rational(1, 4), n6.sp.Rational(3, 4)
        control = n6.univariate_bernstein(poly, t, lo, hi)

        def evaluate(control_points: list[F], u: F) -> F:
            degree = len(control_points) - 1
            return sum(
                c * bernstein_basis(degree, i, u)
                for i, c in enumerate(control_points)
            )

        u = F(2, 5)
        x = F(1, 4) + F(1, 2) * u
        expected = F(2) + F(3) * x + F(5) * x**2
        self.assertEqual(evaluate(control, u), expected)
        left, right = n6.split_univariate(control)
        self.assertEqual(evaluate(left, u), evaluate(control, u / 2))
        self.assertEqual(evaluate(right, u), evaluate(control, (1 + u) / 2))

    def test_tensor_conversion_with_zero_degree_axis(self) -> None:
        expr = 1 + 2 * n6.X + 3 * n6.Z + 4 * n6.X * n6.Z + 5 * n6.U**2
        poly = n6.sp.Poly(expr, *n6.VARS, domain=n6.sp.QQ)
        arr, degrees = n6.power_to_bernstein(poly)
        self.assertEqual(degrees, (1, 1, 2, 0))
        parameters = (F(2, 5), F(1, 3), F(3, 7), F(5, 9))
        x, z, u, _v = parameters
        expected = F(1) + 2 * x + 3 * z + 4 * x * z + 5 * u**2
        self.assertEqual(eval_tensor(arr, parameters), expected)

        left, right = n6.split_axis(arr, 2)
        child_parameters = (parameters[0], parameters[1], F(2, 5), parameters[3])
        parent_left = (parameters[0], parameters[1], F(1, 5), parameters[3])
        parent_right = (parameters[0], parameters[1], F(7, 10), parameters[3])
        self.assertEqual(eval_tensor(left, child_parameters), eval_tensor(arr, parent_left))
        self.assertEqual(eval_tensor(right, child_parameters), eval_tensor(arr, parent_right))


if __name__ == "__main__":
    unittest.main()
