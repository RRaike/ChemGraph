import pkgutil
import importlib

# === Automatically import all modules in this package === #
for loader, module_name, is_pkg in pkgutil.iter_modules(__path__):
    importlib.import_module(f".{module_name}", package=__name__)