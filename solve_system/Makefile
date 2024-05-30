TARGET = main
CC = g++
OBJS = normal_form.o main.o
CFLAGS = -I/home/alemoine/Documents/pml/ntl-extras/include
LDFLAGS = /home/alemoine/Documents/pml/ntl-extras/lib/*.o -lflint -lntl -pthread -lgmp

$(TARGET): $(OBJS)
	$(CC) -g -Wall -o $@ $^ $(LDFLAGS)

normal_form.o: normal_form.c normal_form.h 
	$(CC) -c -o $@ $< $(CFLAGS)

main.o: main.cpp normal_form.c normal_form.h 
	$(CC) -c -o $@ $< $(CFLAGS)

clean:
	rm $(OBJS)