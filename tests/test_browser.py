from cmcbrowser import CMCBrowser
import pytest


@pytest.fixture
def cmc_browser():
    browser = CMCBrowser()
    return browser


def test_load_snapshot(cmc_browser):
    cmc_browser.load_snapshot(
        model_name="N4e5_rv1_rg8_Z0.02", ss_name="initial.snap0158.dat.gz"
    )


def test_list_models(cmc_browser):
    cmc_browser.list_models()


def test_list_snapshots(cmc_browser):
    cmc_browser.list_snapshots(model_name="N4e5_rv1_rg8_Z0.02")
