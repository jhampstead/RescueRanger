# Define the compiler
CC = gcc

# Define the directory containing the source files
SRC_DIR = src

# Define the directory for the output files (optional)
OBJ_DIR = obj

# Collect all .c files in the source directory
SRC_FILES = $(shell find $(SRC_DIR) -name '*.c')

# Define the names of the object files
OBJ_FILES = $(patsubst $(SRC_DIR)/%.c,$(OBJ_DIR)/%.o,$(SRC_FILES))

# Define the name of the executable
TARGET = annotate_ptv_rescue

# Default rule
all: $(TARGET)

# Rule to link the executable
$(TARGET): $(OBJ_FILES)
	$(CC) -g -o $@ $^ bcftools/gff.o bcftools/regidx.o -Ibcftools -lhts -lz

# Rule to compile each source file into an object file
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c | $(OBJ_DIR)
	@mkdir -p $(dir $@)
	$(CC) -g -c -o $@ $< -DLOG_USE_COLOR 

# Rule to create the object directory if it doesn't exist
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

# Rule to clean up the build artifacts
clean:
	rm -rf $(OBJ_DIR) $(TARGET)

.PHONY: all clean

