import os
import sys
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext

# Define the three C modules
extensions = [
    Extension(
        name="voidfindertk.zobov.zobov_loader",
        sources=["voidfindertk/zobov/zobov_loader.c"],
        extra_compile_args=["-fPIC", "-O3"],
    ),
    Extension(
        name="voidfindertk.zobov.zones_in_void",
        sources=["voidfindertk/zobov/zones_in_void.c"],
        extra_compile_args=["-fPIC", "-O3"],
    ),
    Extension(
        name="voidfindertk.zobov.tracers_in_zones",
        sources=["voidfindertk/zobov/tracers_in_zones.c"],
        extra_compile_args=["-fPIC", "-O3"],
    ),
]

class BuildFailed(Exception):
    pass

class CustomBuildExt(build_ext):
    def run(self):
        try:
            super().run()
        except Exception as e:
            raise BuildFailed(f"Compilation failed: {str(e)}")

    def build_extension(self, ext):
        try:
            super().build_extension(ext)
            print(f"Successfully built {ext.name}")
        except Exception as e:
            print(f"Error building {ext.name}: {str(e)}", file=sys.stderr)
            raise