# Makefile — build + versioned, structured install for the rmblastn Rust port.
#
# Layout produced by `make install`:
#   $(PREFIX)/$(NAME)-$(VERSION)/bin/rmblastn      the compiled binary
#   $(PREFIX)/$(NAME)-$(VERSION)/bin/dustmasker    NCBI-dustmasker-compatible CLI
#   $(PREFIX)/$(NAME)-$(VERSION)/wrappers/*        RepeatMasker-facing wrapper scripts
#   $(PREFIX)/$(NAME)-$(VERSION)/*.md              documentation
#   $(PREFIX)/$(NAME)-$(VERSION)/matrices/*        scoring matrices (opt-in, see MATRIX_SRC)
#
# Examples:
#   make                          # build (release)
#   make install PREFIX=/opt      # -> /opt/rmblast-0.1.0/...
#   make install PREFIX=/opt MATRIX_SRC=../   # also install ../*.matrix
#   make dist PREFIX=/tmp/stage   # build a versioned tarball
#   make version                  # print the version from Cargo.toml

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
MANIFEST    = rmblastn/Cargo.toml
# Directory holding *.matrix files to bundle.  Empty by default (opt-in) since the
# canonical matrices live outside the source tree; set on the command line.
MATRIX_SRC ?=

# Read the package version from Cargo.toml at parse time (the only top-level
# `version =` in the binary crate's manifest; deps use workspace/path).
VERSION := $(shell sed -n 's/^version *= *"\([^"]*\)".*/\1/p' $(MANIFEST) | head -1)
$(if $(VERSION),,$(error could not read version from $(MANIFEST)))

# Versioned package root (DESTDIR-aware for staged packaging).
PKGDIR = $(DESTDIR)$(PREFIX)/$(NAME)-$(VERSION)

.PHONY: all build test version install uninstall dist clean

all: build

version:
	@echo $(VERSION)

build:
	$(CARGO) build --release

test:
	$(CARGO) test --release

# Install into a top-level <prefix>/<name>-<version>/ hierarchy.
install: build
	$(INSTALL) -d $(PKGDIR)/bin
	$(INSTALL) -m 0755 target/release/$(BIN) $(PKGDIR)/bin/$(BIN)
	@for b in $(EXTRA_BINS); do \
	  $(INSTALL) -m 0755 target/release/$$b $(PKGDIR)/bin/$$b && echo "  bin       <- $$b"; \
	done
	$(INSTALL) -d $(PKGDIR)/wrappers
	$(INSTALL) -m 0755 wrappers/* $(PKGDIR)/wrappers/
	$(INSTALL) -m 0644 *.md $(PKGDIR)/
	@if [ -n "$(MATRIX_SRC)" ] && [ -d "$(MATRIX_SRC)" ]; then \
	  $(INSTALL) -d $(PKGDIR)/matrices; \
	  $(INSTALL) -m 0644 $(MATRIX_SRC)/*.matrix $(PKGDIR)/matrices/ && \
	  echo "  matrices  <- $(MATRIX_SRC)/*.matrix"; \
	else \
	  echo "  matrices  (skipped: set MATRIX_SRC=<dir> to bundle *.matrix files)"; \
	fi
	@echo "Installed $(NAME) $(VERSION) -> $(PKGDIR)"

uninstall:
	rm -rf $(PKGDIR)

# Versioned tarball of the installed hierarchy: <name>-<version>.tar.gz
dist: install
	tar -C $(DESTDIR)$(PREFIX) -czf $(NAME)-$(VERSION).tar.gz $(NAME)-$(VERSION)
	@echo "Wrote $(NAME)-$(VERSION).tar.gz"

clean:
	$(CARGO) clean
