# Makefile for pion photoproduction Fortran 77 implementation

# Compiler and flags
FC = gfortran
FFLAGS = -O2

# Objects for minimal test
OBJS_TEST = constants.o physics_funcs.o test_getrawm.o

# Objects for full implementation
OBJS_FULL = constants.o vector_ops.o legendre.o block_data.o said_parser.o physics_funcs.o getrawm.o test_getrawm.o

# Targets
TEST_TARGET = test_minimal
FULL_TARGET = test_getrawm

# Default target
all: $(FULL_TARGET)

# Build the minimal test executable
$(TEST_TARGET): $(OBJS_TEST)
	$(FC) $(FFLAGS) -o $@ $(OBJS_TEST)

# Build the full executable
$(FULL_TARGET): $(OBJS_FULL)
	$(FC) $(FFLAGS) -o $@ $(OBJS_FULL)

# Object file rules
%.o: %.f
	$(FC) $(FFLAGS) -c $<

# Dependencies
constants.o: constants.f
vector_ops.o: vector_ops.f constants.o
legendre.o: legendre.f
said_parser.o: said_parser.f parseConstants.def constants.o
physics_funcs.o: physics_funcs.f constants.o vector_ops.o legendre.o said_parser.o
getrawm.o: getrawm.f constants.o vector_ops.o legendre.o said_parser.o physics_funcs.o
test_getrawm.o: test_getrawm.f constants.o

# Clean rule
clean:
	rm -f $(OBJS_TEST) $(OBJS_FULL) $(TEST_TARGET) $(FULL_TARGET)