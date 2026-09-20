import numpy as np
import pytest

from pygcgc import IonExpressionError, fold, parse_ion_expr
from tests.conftest import make_run


@pytest.mark.parametrize("text, num, den, out", [
    ("", None, None, "TIC"),
    ("tic", None, None, "TIC"),
    ("178", (178,), None, "178"),
    ("178 + 252", (178, 252), None, "178+252"),
    ("200:203", (200, 201, 202, 203), None, "200:203"),
    ("203:200+57", (200, 201, 202, 203, 57), None, "200:203+57"),
    ("57/71", (57,), (71,), "57/71"),
    ("57+71/85", (57, 71), (85,), "57+71/85"),
    ("TIC/57", None, (57,), "TIC/57"),
    ("57+57", (57,), None, "57"),
])
def test_parse(text, num, den, out):
    e = parse_ion_expr(text)
    assert e.numerator == num and e.denominator == den
    assert e.text() == out
    assert parse_ion_expr(e.text()) == e     # text round-trips


@pytest.mark.parametrize("bad", ["abc", "57//71", "57/71/85", "5.5", "TIC/TIC", "57+"])
def test_bad(bad):
    with pytest.raises(IonExpressionError):
        parse_ion_expr(bad)


def test_images():
    run = make_run([(62.0, 1.2, 50, 1e5), (62.0, 1.2, 51, 5e4)])
    f = fold(run, 4.0)
    np.testing.assert_allclose(parse_ion_expr("50+51").image(run, f),
                               f.ion_image(run, 50) + f.ion_image(run, 51))
    ratio = parse_ion_expr("51/50").image(run, f)
    peak = np.nanargmax(f.ion_image(run, 50))
    assert ratio.ravel()[peak] == pytest.approx(0.5, rel=1e-6)
    assert np.isnan(ratio).mean() > 0.9          # empty pixels aren't turned into ratios
    with pytest.raises(KeyError):
        parse_ion_expr("999").validate(run)
