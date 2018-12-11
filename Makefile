CXX = g++
CXXFLAGS += -O2 -Wall -Wextra -Werror -I.
ROOTFLAGS := `root-config --cflags --glibs`

SRCDIR = ./src
BINDIR = ./bin
BUILDDIR = ./build
PDFDIR = ./pdfDir

SRCS = $(wildcard $(SRCDIR)/*.C)
EXES = $(patsubst $(SRCDIR)/%.C,$(BINDIR)/%.exe,$(SRCS))
DEPS = $(patsubst $(SRCDIR)/%.C,$(BUILDDIR)/%.d,$(SRCS))

.PHONY: all clean mkdir

all: $(EXES) mkdir

$(BINDIR)/%.exe: $(SRCDIR)/%.C
	@mkdir -p $(BINDIR) $(BUILDDIR)
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -MMD -MF $(BUILDDIR)/$(*F).d $< -o $@

mkdir:
	@mkdir -p $(PDFDIR)

clean:
	@$(RM) $(EXES) $(DEPS)
	@rm -rf $(BINDIR)
	@rm -rf $(BUILDDIR)
	rm -f *~
	rm -f \#*.*#
	rm -f $(PWD)/include/#*.*#
	rm -f $(PWD)/include/*~
	rm -f $(PWD)/src/#*.*#
	rm -f $(PWD)/src/*~

-include $(DEPS)
