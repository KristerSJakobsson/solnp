import os
import sys
import unittest

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from examples.example_alkyla import solve_alkyla
from examples.example_box import solve_box
from examples.example_rosen_suzuki import solve_rozen_suzuki


class TestExtension(unittest.TestCase):

    def test_can_solve(self):
        # Verify it does  not raise any exception
        solve_box()
        solve_alkyla()
        solve_rozen_suzuki()
