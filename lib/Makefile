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

# Platform
ifeq ($(OS), Windows_NT)
WIN32 := 1
else
KERNEL := $(shell uname -s)
ifeq ($(KERNEL), Darwin)
OSX := 1
endif
ifeq ($(KERNEL), Linux)
LINUX := 1
endif
endif

ifdef WIN32
RM := del /F 2> NUL
else
RM := rm -f
endif

# Configuration
MAJ_VER         := 0
MIN_VER         := 1
REV_VER         := 0
LIB_VER         := $(MAJ_VER).$(MIN_VER).$(REV_VER)

CC              := gcc
PREFIX          ?= /usr/local
LIBRARY_NAME    := libordpat

ifdef WIN32
TARGET          := $(LIBRARY_NAME).dll
SHORTNAME       := $(TARGET)
endif

ifdef OSX
TARGET          := $(LIBRARY_NAME).dylib
SHORTNAME       := $(TARGET)
endif

ifdef LINUX
TARGET          := $(LIBRARY_NAME).so.$(LIB_VER)
SHORTNAME       := $(LIBRARY_NAME).so.$(MAJ_VER)
endif

ifndef TARGET
$(error Unsupported target platform)
endif

CFLAGS          := -std=c99 -fPIC -Wall -Wextra -pedantic -DBUILD_LIBORDPAT
CFLAGS_DEBUG    := -g
CFLAGS_RELEASE  := -O3 -fvisibility=hidden

ifdef OSX
LDFLAGS         := -dynamiclib -current_version $(LIB_VER) \
                   -compatibility_version $(MAJ_VER)
LDFLAGS_RELEASE := -O3
else
LDFLAGS         := -shared -Wl,-soname,$(LIBRARY_NAME).so.$(MAJ_VER)
LDFLAGS_RELEASE := -O3 -s
endif

LDFLAGS_DEBUG   := -g

# Architecture
ARCH ?= x64
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
DEBUG ?= 0
ifeq ($(DEBUG), 1)
CFLAGS  += $(CFLAGS_DEBUG)
LDFLAGS += $(LDFLAGS_DEBUG)
else
CFLAGS  += $(CFLAGS_RELEASE)
LDFLAGS += $(LDFLAGS_RELEASE)
endif

# Targets
$(TARGET): ordpat.o
	@$(CC) $(LDFLAGS) -o $@ ordpat.o

ordpat.o: ordpat.c ordpat.h
	@$(CC) $(CFLAGS) -c -o $@ ordpat.c

clean:
	@$(RM) *.o

distclean: clean
	@$(RM) $(LIBRARY_NAME).*

install: $(TARGET)
ifndef WIN32
	@install -d $(PREFIX)/include/
	@install --mode=644 ./ordpat.h $(PREFIX)/include/
	@install -d $(PREFIX)/lib/
	@install $(TARGET) $(PREFIX)/lib/
	@ln -fs $(PREFIX)/lib/$(TARGET) $(PREFIX)/lib/$(LIBRARY_NAME).so
	@ldconfig
endif

uninstall:
ifndef WIN32
	@$(RM) $(PREFIX)/include/ordpat.h
	@$(RM) $(PREFIX)/lib/$(TARGET)
	@if ! [ -e $$(readlink $(PREFIX)/lib/$(LIBRARY_NAME).so) ]; then \
	$(RM) $(PREFIX)/lib/$(LIBRARY_NAME).so; \
	fi
	@ldconfig
endif

.DEFAULT_GOAL := $(TARGET)
.PHONY: clean distclean install uninstall
