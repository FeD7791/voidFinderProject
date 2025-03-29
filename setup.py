from setuptools import setup
from setuptools.command.build import build
import subprocess

import os


class MakeRun(build):
    def run(self):
        # Run `make` before the build process
        makefile_dir = os.path.join(os.path.dirname(__file__), "voidfindertk/zobov")
        try:
            subprocess.check_call(["make"], cwd=makefile_dir)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Make failed with error: {e}")
        
        # Continue with the normal build process
        super().run()



setup(
    cmdclass={"build": MakeRun}
)