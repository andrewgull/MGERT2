import logging
import sys


def setup_logging(log_file=None, logger_name=None):
    """Set up logging to a file when provided, otherwise fall back to stdout."""
    if log_file:
        handlers = [logging.FileHandler(log_file, mode="a")]
    else:
        handlers = [logging.StreamHandler(sys.stdout)]

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=handlers,
    )
    return logging.getLogger(logger_name)
