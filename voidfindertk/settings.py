import datetime as dt
import json
import os
import pathlib

from . import __version__ as VERSION
from .utils.bunch import Bunch

# =============================================================================
# DEFAULT CONF
# =============================================================================
_EMPTY_CONF = {
    "void_finder_tk_version": VERSION,
    "created_at": None,
    "paths": {
        "zobov_path": "",
        "popcorn_path": "",
    },
}

# =============================================================================
# API
# =============================================================================


def create_empty_conf(path):
    conf = _EMPTY_CONF.copy()
    conf["created_at"] = dt.datetime.now(dt.timezone.utc).isoformat()
    with open(path, "w") as fp:
        json.dump(conf, fp, indent=2)


def read_conf(path):
    with open(path) as fp:
        return Bunch(path, json.load(fp))


# =============================================================================
# CONSTANTS
# =============================================================================

USER_HOME_PATH = pathlib.Path(os.path.expanduser("~"))

DEFAULT_CONF_PATH = USER_HOME_PATH / ".voidfindertk.json"

if not DEFAULT_CONF_PATH.exists():
    print("Creating new configurarion...")
    create_empty_conf(DEFAULT_CONF_PATH)
    print(f"Please configure {DEFAULT_CONF_PATH}")

SETTINGS = read_conf(DEFAULT_CONF_PATH)
