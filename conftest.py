"""Repository-wide test behaviour: a test that opens an absent frozen product is skipped, with how to restore it.

Frozen products (shared/frozen_products.json) are large results of closed studies that a checkout links to, restores
from Git's history, or goes without. Where one is absent, the tests that read it cannot run; they are reported as
skipped with the restoring command, as tests elsewhere skip when an ignored input is not at hand. A test that fails
for any other reason, or that reads a frozen product that is present, is reported as it is.
"""
import pytest

from shared import frozen_products


def _reason(path):
    return f'frozen product absent: {path} ({frozen_products.RESTORE})'


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(item, call):
    outcome = yield
    report = outcome.get_result()
    if report.failed and call.excinfo is not None:
        path = frozen_products.absent_frozen_product(call.excinfo.value)
        if path is not None:
            report.outcome = 'skipped'
            report.longrepr = (str(item.path), item.location[1] or 0, 'Skipped: ' + _reason(path))


@pytest.hookimpl(hookwrapper=True)
def pytest_make_collect_report(collector):
    outcome = yield
    report = outcome.get_result()
    call = getattr(report, 'call', None)
    excinfo = getattr(call, 'excinfo', None)
    if report.failed and excinfo is not None:
        path = frozen_products.absent_frozen_product(excinfo.value)
        if path is not None:
            report.outcome = 'skipped'
            report.longrepr = (str(collector.path), 0, 'Skipped: ' + _reason(path))
