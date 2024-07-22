import contextlib
import os


try:
    chdir = contextlib.chdir
except AttributeError:

    class chdir(contextlib.AbstractContextManager):
        """Non thread-safe context manager to change the current working \
        directory.

        """

        def __init__(self, path):
            self.path = path
            self._old_cwd = []

        def __enter__(self):
            self._old_cwd.append(os.getcwd())
            os.chdir(self.path)

        def __exit__(self, *excinfo):
            os.chdir(self._old_cwd.pop())
