__version__ = "1.0"

import contextlib
from typing import Optional, List, TextIO


def main(argv: Optional[List[str]] = None, stream: Optional[TextIO] = None) -> int:

    from topology_cmd.topology_cmd import main_app

    #with contextlib.ExitStack() as ctx:
    return main_app(__version__)