import sys
from pathlib import Path

# this allows us to import files without the .py extension
# I don't entirely understand the need for types, but it seems
# like that is needed since the migration from imp to
# importlib
import types
import importlib.machinery
loader = importlib.machinery.SourceFileLoader('common', 'workflow/rules/common.smk')
mod = types.ModuleType(loader.name)
loader.exec_module(mod)

sys.modules["common"] = mod
from common import *
