import os
import importlib
import sys

def import_path(path):
    module_name = os.path.basename(path).replace('-', '_')
    spec = importlib.util.spec_from_loader(
        module_name,
        importlib.machinery.SourceFileLoader(module_name, path)
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    sys.modules[module_name] = module
    return module

common = import_path('workflow/rules/common.smk')

def test_make_assembly_split_names():
    assert common.make_assembly_split_names(4) == ["001", "002", "003", "004"]
    assert common.make_assembly_split_names(100)[99] == "100"
