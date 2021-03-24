# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SOURCEDIR     = source
BUILDDIR      = build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@mv source/ Makefile INSTRUCTIONS.rst ..
	@if [ -f ".nojekyll" ]; then rm .nojekyll; fi
	@if [ -f ".buildinfo" ]; then rm .buildinfo; fi
	@if [ -d ".doctrees" ]; then rm -r .doctrees; fi
	@if [ -f "index.html" ]; then rm -r *; fi
	@mv ../source/ ../Makefile ../INSTRUCTIONS.rst .
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
	@mv $(BUILDDIR)/$@/* $(BUILDDIR)/$@/.nojekyll $(BUILDDIR)/$@/.buildinfo .
	@mv $(BUILDDIR)/doctrees .doctrees
	@rm -r $(BUILDDIR)
