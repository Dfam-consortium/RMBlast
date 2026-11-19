# Makefile — build + versioned, structured install for the rmblastn Rust port.
#
# Layout produced by `make install`:
#   $(PREFIX)/$(NAME)-$(VERSION)/bin/rmblastn      the compiled binary
#   $(PREFIX)/$(NAME)-$(VERSION)/bin/dustmasker    NCBI-dustmasker-compatible CLI
#   $(PREFIX)/$(NAME)-$(VERSION)/wrappers/*        RepeatMasker-facing wrapper scripts
#   $(PREFIX)/$(NAME)-$(VERSION)/verify_wrappers/* same, but A/B-checking rmblastn
#   $(PREFIX)/$(NAME)-$(VERSION)/README.md         documentation (see DOCS)
#   $(PREFIX)/$(NAME)-$(VERSION)/matrices/*        scoring matrices (opt-in, see MATRIX_SRC)
#
# Examples:
#   make                          # build (release)
#   make install PREFIX=/opt      # -> /opt/rmblast-<version>/...
#   make install PREFIX=/opt MATRIX_SRC=../   # also install ../*.matrix
#   make dist                     # build a versioned tarball (staged under ./stage)
#   make version                  # print the workspace version from Cargo.toml

PREFIX     ?= /usr/local
DESTDIR    ?=
NAME       ?= rmblast
CARGO      ?= cargo
INSTALL    ?= install
BIN         = rmblastn
# Additional binaries installed alongside rmblastn.  `dustmasker` is a drop-in
# replacement for the NCBI application of the same name (FASTA in, text out);
# installing it shadows NCBI's copy for anything that finds it on PATH first.
EXTRA_BINS  = dustmasker
# The wrapper scripts that ship, by name: a file dropped into wrappers/ or
# verify_wrappers/ does not end up in a release by accident.
WRAPPERS    = rmblastn makeblastdb blastdbcmd blastdb_aliastool dustmasker
MANIFEST    = Cargo.toml
# Documentation shipped in a binary release.  Listed by name rather than globbed
# as *.md, so that a maintainer note dropped at the top level does not silently
# end up in a release tarball.
DOCS        = README.md LICENSE
# Directory holding *.matrix files to bundle.  Empty by default (opt-in) since the
# canonical matrices live outside the source tree; set on the command line.
MATRIX_SRC ?=

# Read the workspace version from the root Cargo.toml at parse time.  It lives
# in [workspace.package] and every crate inherits it, so rmblastn and dustmasker
# always report the same version.  The sed range confines the match to that table.
VERSION := $(shell sed -n '/^\[workspace\.package\]/,/^\[/s/^version *= *"\([^"]*\)".*/\1/p' $(MANIFEST) | head -1)
$(if $(VERSION),,$(error could not read version from $(MANIFEST)))

# Versioned package root (DESTDIR-aware for staged packaging).
PKGDIR = $(DESTDIR)$(PREFIX)/$(NAME)-$(VERSION)

.PHONY: all build test test-full version check-release install uninstall dist clean

all: build

version:
	@echo $(VERSION)

build:
	$(CARGO) build --release

# Before tagging a release (see the script for what it checks):
#   make check-release TAG=X.Y.Z
check-release: build
	bash .github/scripts/check_release.sh $(VERSION) "$(TAG)" $(BIN) $(EXTRA_BINS)

test:
	$(CARGO) test --release

# Also runs the #[ignore]d chr22 seed and ungapped-hit regression tests
# (rmblast-lib/tests/), which compare against NCBI-derived references and add
# about 40 s, the CLI tests that need the `diagnostics` test hooks (built in
# their own target dir so target/release keeps the normal binary), and the
# wrapper tests, which use target/release.
test-full:
	$(CARGO) test --release -- --include-ignored
	CARGO_TARGET_DIR=target/diagnostics $(CARGO) test --release -p rmblastn --features diagnostics --test cli
	bash tests/wrappers.sh

# Install into a top-level <prefix>/<name>-<version>/ hierarchy.
#
# verify_wrappers/ ships alongside wrappers/: it is only useful where NCBI
# rmblastn is also installed, but it costs nothing and it has to stay in sync
# with wrappers/.  Both directories locate the binaries relative to themselves
# (../bin here, ../target/release in a source checkout), so neither needs
# editing after install.
install: build
	$(INSTALL) -d $(PKGDIR)/bin
	$(INSTALL) -m 0755 target/release/$(BIN) $(PKGDIR)/bin/$(BIN)
	@for b in $(EXTRA_BINS); do \
	  $(INSTALL) -m 0755 target/release/$$b $(PKGDIR)/bin/$$b && echo "  bin       <- $$b"; \
	done
	$(INSTALL) -d $(PKGDIR)/wrappers
	$(INSTALL) -m 0755 $(addprefix wrappers/,$(WRAPPERS)) $(PKGDIR)/wrappers/
	$(INSTALL) -d $(PKGDIR)/verify_wrappers
	$(INSTALL) -m 0755 $(addprefix verify_wrappers/,$(WRAPPERS)) $(PKGDIR)/verify_wrappers/
	$(INSTALL) -m 0644 verify_wrappers/README.md $(PKGDIR)/verify_wrappers/
	@for d in $(DOCS); do \
	  $(INSTALL) -m 0644 $$d $(PKGDIR)/ && echo "  doc       <- $$d"; \
	done
	@if [ -n "$(MATRIX_SRC)" ] && [ -d "$(MATRIX_SRC)" ]; then \
	  $(INSTALL) -d $(PKGDIR)/matrices; \
	  $(INSTALL) -m 0644 $(MATRIX_SRC)/*.matrix $(PKGDIR)/matrices/ && \
	  echo "  matrices  <- $(MATRIX_SRC)/*.matrix"; \
	else \
	  echo "  matrices  (skipped: set MATRIX_SRC=<dir> to bundle *.matrix files)"; \
	fi
	@echo "Installed $(NAME) $(VERSION) -> $(PKGDIR)"
	@for b in $(BIN) $(EXTRA_BINS); do echo "  $$($(PKGDIR)/bin/$$b --version)"; done

uninstall:
	rm -rf $(PKGDIR)

# Versioned tarball <name>-<version>.tar.gz of the installed hierarchy, built
# from a staging install: under DESTDIR if given (the release workflow passes
# one), else ./stage.  It never installs into PREFIX itself, which defaults
# to /usr/local.
DIST_DESTDIR = $(if $(DESTDIR),$(DESTDIR),$(CURDIR)/stage)
dist: build
	$(MAKE) install DESTDIR=$(DIST_DESTDIR)
	tar -C $(DIST_DESTDIR)$(PREFIX) -czf $(NAME)-$(VERSION).tar.gz $(NAME)-$(VERSION)
	@echo "Wrote $(NAME)-$(VERSION).tar.gz (staged in $(DIST_DESTDIR)$(PREFIX))"

clean:
	$(CARGO) clean
