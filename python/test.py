import unittest

from .example_alkyla import solve_alkyla
from .example_box import solve_box
from .example_rosen_suzuki import solve_rozen_suzuki


class TestExtension(unittest.TestCase):
    def test_can_solve(self):
        # Verify it does  not raise any exception
        solve_box()
        solve_alkyla()
        solve_rozen_suzuki()


if __name__ == '__main__':
    unittest.main()
