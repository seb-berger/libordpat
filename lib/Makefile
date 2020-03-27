# libordpat - A function library for extracting and encoding
#             ordinal patterns from time series.
#
# For a detailed description of the algorithms used, see the article:
#
# [1] Berger S, Kravtsiv A, Schneider G, Jordan D.
#     Teaching Ordinal Patterns to a Computer: Efficient Encoding
#     Algorithms Based on the Lehmer Code. Entropy. 2019; 21(10):1023.
#     https://doi.org/10.3390/e21101023

# Copyright (c) 2019, Sebastian Berger.
#
# Klinikum rechts der Isar der
# Technischen Universitaet Muenchen
# Munich, Germany
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in
#       the documentation and/or other materials provided with the
#       distribution.
#     * Neither the names of the copyright holders nor the names of its
#       contributors may be used to endorse or promote products derived
#       from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR OR
# THE KLINIKUM RECHTS DER ISAR BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
# AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
# THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGE.

# Configuration ################################################################

CC        := gcc
PREFIX    ?= /usr/local
ARCH      ?= x64
DEBUG     ?= 0
CFLAGS    := -std=c99 -Wall -Wextra -pedantic
LDFLAGS   := # will be populated below
LIB_NAME  := libordpat
LIB_OBJ   := ordpat.o
MAKEFLAGS := rR # do not use built-in rules or variables
MAJ_VER   := 0
MIN_VER   := 1
REV_VER   := 0
LIB_VER   := $(MAJ_VER).$(MIN_VER).$(REV_VER)


# Prepare Build Environment ####################################################

$(LIB_OBJ): CFLAGS += -DBUILD_LIBORDPAT

# Architecture
ifeq ($(ARCH), x86)
    CFLAGS  += -m32
    LDFLAGS += -m32
else
    ifeq ($(ARCH), x64)
        CFLAGS  += -m64
        LDFLAGS += -m64
    else
        $(error Unsupported target architecture: '$(ARCH)')
    endif
endif

# Debug/Release
ifeq ($(DEBUG), 0)
    CFLAGS  += -O3
    LDFLAGS += -O3
else
    CFLAGS  += -g -O0
    LDFLAGS += -g -O0
endif

# OS detection
override undefine LINUX  # make sure that all three
override undefine MACOS  # of these are initially
override undefine WIN32  # undefined

ifeq ($(OS), Windows_NT)
    WIN32 := 1
else
    KERNEL := $(shell uname -s)
    ifeq ($(KERNEL), Darwin)
        MACOS := 1
    else
        ifeq ($(KERNEL), Linux)
            LINUX := 1
        else
            $(error Unsupported target platform: '$(KERNEL)')
        endif
    endif
endif

# Settings for Linux
ifdef LINUX
    TARGET := $(LIB_NAME).so.$(LIB_VER)
    SHORTNAME := $(LIB_NAME).so.$(MAJ_VER)
    $(LIB_OBJ): CFLAGS  += -fPIC
    $(TARGET):  LDFLAGS += -shared -Wl,-soname,$(LIB_NAME).so.$(MAJ_VER)

    ifneq ($(DEBUG), 1)
        $(LIB_OBJ): CFLAGS += -fvisibility=hidden
        $(TARGET): LDFLAGS += -s
    endif

    RM := rm -f
endif

# Settings for macOS
ifdef MACOS
    TARGET := $(LIB_NAME).dylib
    SHORTNAME := $(TARGET)
    $(LIB_OBJ): CFLAGS  += -fPIC
    $(TARGET):  LDFLAGS += -dynamiclib -current_version $(LIB_VER) \
                           -compatibility_version $(MAJ_VER)

    ifneq ($(DEBUG), 1)
        $(LIB_OBJ): CFLAGS += -fvisibility=hidden
    endif

    RM := rm -f
endif

# Settings for Windows
ifdef WIN32
    TARGET := $(LIB_NAME).dll
    SHORTNAME := $(TARGET)
    $(TARGET): LDFLAGS += -shared

    ifneq ($(DEBUG), 1)
        $(TARGET): LDFLAGS += -s
    endif

    RM := del /F 2> NUL
endif


# Targets ######################################################################

.DEFAULT_GOAL := $(TARGET)

$(TARGET): $(LIB_OBJ)
	@$(CC) $(LDFLAGS) -o $@ $?

ordpat.o: ordpat.c ordpat.h
	@$(CC) $(CFLAGS) -o $@ -c ordpat.c

clean:
	@$(RM) $(LIB_OBJ)

distclean: clean
	@$(RM) $(TARGET)

install: $(TARGET)
# TODO: Add macOS support
ifdef LINUX
	@install -d $(PREFIX)/include/
	@install --mode=644 ./ordpat.h $(PREFIX)/include/
	@install -d $(PREFIX)/lib/
	@install $(TARGET) $(PREFIX)/lib/
	@ln -fs $(PREFIX)/lib/$(TARGET) $(PREFIX)/lib/$(LIB_NAME).so
	@ldconfig
endif

uninstall:
# TODO: Add macOS support
ifdef LINUX
	@$(RM) $(PREFIX)/include/ordpat.h
	@$(RM) $(PREFIX)/lib/$(TARGET)
	@if ! [ -e $$(readlink $(PREFIX)/lib/$(LIB_NAME).so) ]; then \
	    $(RM) $(PREFIX)/lib/$(LIB_NAME).so; \
	fi
	@ldconfig
endif

.PHONY: clean distclean install uninstall
