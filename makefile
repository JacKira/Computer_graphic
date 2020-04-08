CC = gcc
CFLAGS = -g -Wall -pedantic
LFLAGS = -lm

.PHONY:  render main.o tga.o model.o  

render: main.o tga.o model.o 
	$(CC) $(LFLAGS) -o $@ $^ 

main.o: main.c tga.h model.h
	$(CC) -c $(CFLAGS) -o $@ $<

tga.o:tga.c tga.h
	$(CC) -c $(CFLAGS) -o $@ $<

model.o:model.c model.h tga.h
	$(CC) -c $(CFLAGS) -o $@ $<

