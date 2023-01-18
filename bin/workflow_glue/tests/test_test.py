"""A dummy test."""

import argparse

from workflow_glue import cnv_plot


def test():
    """Just showing that we can import using the workflow-glue."""
    assert isinstance(cnv_plot.argparser(), argparse.ArgumentParser)
