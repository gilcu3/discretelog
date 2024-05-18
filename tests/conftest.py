#! /usr/bin/env python

import pytest
import random


def pytest_addoption(parser):
    parser.addoption(
        "--slow", action="store_true", default=False, help="run slow tests"
    )
    parser.addoption(
        "--perf", action="store_true", default=False, help="run performance tests"
    )


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: mark test as slow to run")
    config.addinivalue_line("markers", "perf: mark test as perf to run")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--slow"):
        skip_slow = None
    else:
        skip_slow = pytest.mark.skip(reason="need --slow option to run")
    if config.getoption("--perf"):
        skip_perf = None
    else:
        skip_perf = pytest.mark.skip(reason="need --perf option to run")

    for item in items:
        if "slow" in item.keywords and skip_slow is not None:
            item.add_marker(skip_slow)
        if "perf" in item.keywords and skip_perf is not None:
            item.add_marker(skip_perf)


@pytest.fixture
def frandom():
    return random.Random(0)
