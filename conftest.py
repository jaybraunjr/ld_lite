"""Put the repo root on sys.path so `pytest` works without installing first."""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
