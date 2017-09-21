# -----------------------------------------------
# Makefile for compiling GLS3D
#
# Author : Aniruddha Panda
# Date   : April 12, 2017
#
# -----------------------------------------------

# Executable/Project Name
TARGET          = ./GLS3D

# Compilers
GCC              = gcc -fopenmp
ICC              = icc -fopenmp

# CompileFlags
CFLAGS           = -g3 -O2
IFLAGS           = -axW -ipo -ansi-alias
DFLAGS           = -g3 -O0

LINKER           = gcc -o -fopenmp

# LinkingFlags
LFLAGS           = -Wall -lm

# ProjectDirectories
SRCDIR           = src
OBJDIR           = obj
INCDIR           = include

# Scan for the files
SOURCES         := $(wildcard $(SRCDIR)/*.c)
INCLUDES        := $(wildcard $(INCDIR)/*.h)
OBJECTS         := $(SOURCES:$(SRCDIR)/%.c=$(OBJDIR)/%.o)
rm               = rm -f

# Recipes
$(TARGET): $(OBJECTS)
	@$(GCC) $(OBJECTS) -o $@ $(LFLAGS)
	@echo "Linking complete!"

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.c
	@$(CC) $(DFLAGS) -c $< -o $@
	@echo "Compiled "$<" successfully!"

.PHONY: clean
clean:
	@$(rm) $(OBJECTS)
	@echo "Cleanup complete!"

.PHONY: remove
remove: clean
	@$(rm) $(TARGET)
	@echo "Executable removed!"
