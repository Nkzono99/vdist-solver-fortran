import subprocess

from setuptools import Command, find_packages, setup
from setuptools.command.build import build


class CustomCommnad(Command):
    def initialize_options(self) -> None:
        pass

    def finalize_options(self) -> None:
        pass

    def run(self) -> None:
        command = "make"
        subprocess.run(command, shell=True)


class CustomBuild(build):
    sub_commands = [("build_custom", None)] + build.sub_commands


setup(
    cmdclass={"build": CustomBuild, "build_custom": CustomCommnad},
    packages=find_packages(),
)
