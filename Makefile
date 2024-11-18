# Copyright Cristina Iordache 314CAa 2023-2024
# compiler setup
CC=gcc
CFLAGS=-Wall -Wextra -std=c99

# define targets
TARGETS = my_octave

build: $(TARGETS)

my_octave: my_octave.c
	$(CC) $(CFLAGS) my_octave.c -lm -o my_octave

pack:
	zip -FSr 314CA_IordacheCristina_Tema2.zip README Makefile *.c *.h

clean:
	rm -f $(TARGETS)

.PHONY: pack clean
