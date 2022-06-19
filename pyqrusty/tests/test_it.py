import pytest

from pyqrusty import py_sum

# content of test_sample.py
def inc(x):
    return x + 1


def test_answer():
    assert inc(3) == 4

def test_sum():
    assert py_sum(3,4) == 7
