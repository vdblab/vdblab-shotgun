import sys
from pathlib import Path

# this allows us to import files without the .py extension
# I don't entirely understand the need for types, but it seems
# like that is needed since the migration from imp to
# importlib
import types
import importlib.machinery

loader = importlib.machinery.SourceFileLoader("common", "workflow/rules/common.smk")
mod = types.ModuleType(loader.name)
loader.exec_module(mod)

sys.modules["common"] = mod
from common import *

# TODO: many of the snakemake input/param functions rely on accessing `config`
# from the global namespace. That won't work defining that globally in this
# test module because the imported modules have a different namespace.

def test_make_assembly_split_names():
    split_names = make_assembly_split_names(10000)
    assert len(str(split_names[0])) == 3, "bad padding!"
    assert len(str(split_names[1001])) == 4, "bad padding!"

def test_bbmap_dedup_params_flags():
    test_conditions = [
        [{"dedup_platform": "SRA"}, "dedupe"],
        [{"dedup_platform": "MiniSeq"}, "dedupe optical dupedist=40"],
        [{"dedup_platform": "NovaSeq"}, "dedupe optical dupedist=12000"],
        [{"dedup_platform": "NextSeq"}, "dedupe optical spany adjacent dupedist=40"],
        [{"dedup_platform": "HiSeq-3000"}, "dedupe optical dupedist=2500"],
        [{"dedup_platform": "HiSeq-4000"}, "dedupe optical dupedist=2500"],
    ]
    config = {}
    for i in test_conditions:
        config = i[0]
        outcome = i[1]
        assert bbmap_dedup_params_flags(None, config) == outcome
        print(f"passed dedup flags test for {config['dedup_platform']}")
