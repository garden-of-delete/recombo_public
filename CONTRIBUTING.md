# Contributing to RECOMBO

## Build System

RECOMBO uses CMake (minimum version 3.16). To build:

```bash
cmake -B build
cmake --build build
```

Production code is compiled into a static library (`recombo_core`) which all executables and tests link against.

## Testing

### Frameworks

Tests use [Google Test](https://github.com/google/googletest) (v1.14.0), fetched automatically via CMake's FetchContent.

### Test Directory Structure

```
test/
  unit/          # Unit tests for individual classes and functions
  regression/    # Regression tests protecting behavioral stability
src/tests/       # Legacy tests (being migrated to test/)
```

### Naming Conventions

- **Test files:** `snake_case_test.cpp`, named after the class or module under test (e.g. `pseudorandom_test.cpp`)
- **Test suites:** PascalCase with `Test` suffix, named after the class or module (e.g. `PseudorandomTest`)
- **Test names:** PascalCase, describe the behavior being verified (e.g. `ProducesDeterministicSequenceFromSeed`)

Example:

```cpp
// pseudorandom_test.cpp
TEST(PseudorandomTest, ProducesDeterministicSequenceFromSeed)
{
    // ...
}
```

### Running Tests

```bash
# All legacy unit tests
./build/src/bin/unitTest

# Regression tests
./build/test/regression/src/bin/regressionTest
```

## Docker

A Dockerfile is provided for reproducible builds on `debian:bookworm-slim`. It builds all targets and runs the unit test suite.

```bash
docker build -t recombo .
```

## Code Style

The codebase contains legacy C code (`cross.c`, `homfly1_21.c`, `polynomesAZ.c`) which uses K&R-style function definitions. These files receive relaxed compiler settings in CMake.

## Pull Requests

- Keep PRs focused on a single concern
- Ensure all existing tests pass before submitting
- Add tests for new functionality
