# cx_Freeze setup file

import sys
from cx_Freeze import setup, Executable

# Dependencies are automatically detected, but it might need fine tuning.
build_exe_options = {"packages": ["pysam", "numpy", "matplotlib"]}

# GUI applications require a different base on Windows (the default is for a
# console application).
base = None
if sys.platform == "win32":
    base = "Win32GUI"

setup(  name = "DiscoPlot",
        version = "1.0.4",
        description = "Discordant read visualisation.",
        options = {"build_exe": build_exe_options},
        executables = [Executable("DiscoPlot.py", base=base)])