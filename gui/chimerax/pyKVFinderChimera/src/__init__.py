# vim: set expandtab shiftwidth=4 softtabstop=4:

from chimerax.core.toolshed import BundleAPI
from chimerax.core.commands import run

# Subclass from chimerax.core.toolshed.BundleAPI and
# override the method for registering commands,
# inheriting all other methods from the base class.
class _MyAPI(BundleAPI):
    api_version = 1     # start_tool called with BundleInfo and
                        # ToolInfo instance (vs. BundleInfo and
                        # tool name when api_version==0 [the default])

    # Override method
    @staticmethod
    def start_tool(session, bi, ti):
        # session is an instance of chimerax.core.session.Session
        # bi is an instance of chimerax.core.toolshed.BundleInfo
        # ti is an instance of chimerax.core.toolshed.ToolInfo

        # This method is called once for each time the tool is invoked.

        # We check the name of the tool, which should match one of the
        # ones listed in bundle_info.xml (without the leading and
        # trailing whitespace), and create and return an instance of the
        # appropriate class from the ``tool`` module.
        if ti.name == "Cavities":
            resolveImports(session)
            from . import kvfinder
            return kvfinder.KVFinder(session, ti.name)
        raise ValueError("trying to start unknown tool: %s" % ti.name)

    @staticmethod
    def get_class(class_name):
        # class_name will be a string
        if class_name == "KVFinder":
            from . import kvfinder
            return kvfinder.KVFinder
        raise ValueError("Unknown class name '%s'" % class_name)

def resolveImports(session):
    erros = False
    try:
        import pyKVFinder
    except:
        session.logger.info("pyKVFinder isn't installed.\nInstalling...")
        run(session, "pip install pyKVFinder")
        erros = True
    try:
        import PyQt5
    except:
        session.logger.info("PyQt5 isn't installed.\nInstalling...")
        run(session, "pip install PyQt5")  
        erros = True

    if erros:
        session.logger.info("Please restart your ChimeraX")

# Create the ``bundle_api`` object that ChimeraX expects.
bundle_api = _MyAPI()

