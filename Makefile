# The platform file must define the following:
# CXX : the C++ compiler
# CPPFLAGS : the preprocessor flags
# CXXFLAGS : the c++ compile flags
# LDFLAGS : the link flags
CXX := g++
CPPFLAGS :=
CXXFLAGS := -pedantic -Wall -Wextra -Winline -Wfatal-errors -O3 -funroll-loops
LDFLAGS :=

CPPFLAGS += -I./lib

LIB_SRCS = $(wildcard lib/*.cc)
LIB_OBJS = $(patsubst %.cc,%.o,$(LIB_SRCS))
LIB_TARGET = libcosmospec.a

TEST_SRCS = $(wildcard test/*.cc)
TEST_OBJS = $(patsubst %.cc,%.o,$(TEST_SRCS))
TEST_TARGET = test/run_tests.ex

# Targets
all: $(LIB_TARGET) $(APP_TARGET) tests

$(LIB_TARGET): $(LIB_OBJS)
	@echo [MAKE] Archiving static library.
	ar cr $(LIB_TARGET) $(LIB_OBJS)
	ranlib $(LIB_TARGET)

$(TEST_TARGET): $(LIB_TARGET) $(TEST_OBJS)
	@echo [MAKE] Linking test runner.
	$(CXX) $(LDFLAGS) $(TEST_OBJS) $(LIB_TARGET) -o $@

.PHONY: tests
tests: $(TEST_TARGET)
	@echo [MAKE] Running tests.
	@./$(TEST_TARGET)

.PHONY: test
test: tests

.PHONY: clean
clean:
	rm -rf $(LIB_OBJS) $(LIB_TARGET) $(TEST_OBJS) $(TEST_TARGET)
