from cmcbrowser import CMCBrowser
import pytest


def test_gzip_snapshot():

    browser = CMCBrowser()

    browser.load_snapshot(
        model_name="N4e5_rv1_rg8_Z0.02",
        ss_name="initial.snap0147.dat.gz",
        distance=5.0,
        mode="dat.gz",
    )


def test_h5_snapshot():

    browser = CMCBrowser()

    browser.load_snapshot(
        model_name="N2e5_rv0.5_rg20_Z0.02",
        ss_name="king.window.snapshots.h5",
        distance=5,
        h5_key="8(t=12.000067Gyr)",
    )


def test_list_models():
    browser = CMCBrowser()
    browser.list_models()


def test_list_snapshots():
    browser = CMCBrowser()
    browser.list_snapshots(model_name="N4e5_rv1_rg8_Z0.02")
