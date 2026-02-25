# ============================
#  Makefile (BK-like layout)
#  - forces gfortran
#  - puts .o/.mod in build/
# ============================

# ---- toolchain (FORCED) ----
FC := gfortran

# ---- build dirs ----
BUILDDIR := build
MODDIR   := $(BUILDDIR)

# ---- target ----
TARGET := phi4_mmc

# ---- mode: release | debug ----
MODE ?= release

# ---- flags ----
FSTD   := -std=f2008
FBASE  := -Wall -frecursive $(FSTD)

ifeq ($(MODE),debug)
  OPT   := -O0 -g
  CHECK := -fcheck=all -fbacktrace
else
  OPT   := -O3
  CHECK :=
endif

FFLAGS  := $(OPT) $(FBASE) $(CHECK)
LDFLAGS :=

# ---- sources (compile order matters due to modules) ----
SRCS := \
  mod_kinds.f90 \
  mod_types_phi4.f90 \
  mod_cli_phi4.f90 \
  mod_input_phi4.f90 \
  mod_random.f90 \
  mod_phi4_model.f90 \
  mod_metropolis_phi4.f90 \
  mod_phi4_init.f90 \
  mod_ensemble_dispatch.f90 \
  mod_io_append_dispatch.f90 \
  mod_io_observables_phi4.f90 \
  mod_io_restart_phi4.f90 \
  main_phi4.f90

# ---- objects ----
OBJS := $(patsubst %.f90,$(BUILDDIR)/%.o,$(SRCS))

# ---- default ----
all: $(TARGET)

# ---- link ----
$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJS)

# ---- ensure build dir ----
$(BUILDDIR):
	@mkdir -p $(BUILDDIR)

# ---- compile rule ----
# -J$(MODDIR): write .mod to build/
# -I$(MODDIR): search .mod in build/
$(BUILDDIR)/%.o: %.f90 | $(BUILDDIR)
	$(FC) $(FFLAGS) -J$(MODDIR) -I$(MODDIR) -c $< -o $@

# ---- utilities ----
clean:
	@rm -f $(BUILDDIR)/*.o $(BUILDDIR)/*.mod *.mod

distclean: clean
	@rm -rf $(BUILDDIR) $(TARGET)

rebuild: distclean all

run: $(TARGET)
	./$(TARGET)

debug:
	@$(MAKE) MODE=debug rebuild

release:
	@$(MAKE) MODE=release rebuild

print:
	@echo "FC=$(FC)"
	@echo "MODE=$(MODE)"
	@echo "FFLAGS=$(FFLAGS)"
	@echo "SRCS=$(SRCS)"

.PHONY: all clean distclean rebuild run debug release print
