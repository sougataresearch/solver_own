import numpy as np

from sougata_solver.geometry import Circle, Rectangle
from sougata_solver.materials import Material


def test_circle_dc_value_equals_area():
    mat = Material("core", 4.0)
    circle = Circle(center=(0.0, 0.0), radius=0.3, material=mat)
    assert np.isclose(circle.fourier_transform(0.0, 0.0), circle.area)


def test_rectangle_dc_value_equals_area():
    mat = Material("core", 4.0)
    rect = Rectangle(center=(0.0, 0.0), halfwidth=(0.2, 0.15), material=mat)
    assert np.isclose(rect.fourier_transform(0.0, 0.0), rect.area)


def test_circle_off_center_dc_value_still_equals_area():
    mat = Material("core", 4.0)
    circle = Circle(center=(0.3, -0.1), radius=0.2, material=mat)
    assert np.isclose(circle.fourier_transform(0.0, 0.0), circle.area)


def test_rectangle_contains_matches_geometry():
    mat = Material("core", 4.0)
    rect = Rectangle(center=(0.0, 0.0), halfwidth=(0.2, 0.1), material=mat)
    assert rect.contains(0.0, 0.0)
    assert rect.contains(0.19, 0.09)
    assert not rect.contains(0.21, 0.0)
    assert not rect.contains(0.0, 0.11)


def test_circle_contains_matches_geometry():
    mat = Material("core", 4.0)
    circle = Circle(center=(0.0, 0.0), radius=0.25, material=mat)
    assert circle.contains(0.0, 0.0)
    assert circle.contains(0.24, 0.0)
    assert not circle.contains(0.26, 0.0)
