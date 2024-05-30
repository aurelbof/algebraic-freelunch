CC = g++ --std=c++17
TARGET = main
OBJS = fast_multivariate_multiplication.o generate_system.o normal_form.o product.o main.o
CFLAGS = -I/home/alemoine/Documents/pml/ntl-extras/include
LDFLAGS = /home/alemoine/Documents/pml/ntl-extras/lib/*.o -lflint -lntl -pthread -lgmp 

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) -g -Wall -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) -c -g -Wall -O3 -o $@ $< $(CFLAGS)

%.o: %.cpp
	$(CC) -c -g -Wall -O3 -o $@ $< $(CFLAGS)

normal_form.o: normal_form.c normal_form.h 

generate_system.o: generate_system.cpp generate_system.h normal_form.c normal_form.h 

product.o: product.c product.h 

fast_multivariate_multiplication.o: fast_multivariate_multiplication.cpp fast_multivariate_multiplication.h product.c product.h 

clean:
	rm $(OBJS)