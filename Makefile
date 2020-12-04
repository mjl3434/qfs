CC=g++
CFLAGS=-g -Wall -std=c++17 -fno-stack-protector
INCLUDES=-I.
TARGET = qfs
LINK_PATH = -L/usr/lib/x86_64-linux-gnu
LINK = -lgmpxx -lgmp

all: $(TARGET)

$(TARGET): main.o qfs.o
	$(CC) -o $@ $^ $(CFLAGS) $(LINK_PATH) $(LINK)

main.o: main.cpp
	$(CC) -o $@ -c $< $(INCLUDES) $(CFLAGS)

qfs.o: qfs.cpp
	$(CC) -o $@ -c $< $(INCLUDES) $(CFLAGS)

clean:
	rm -f $(TARGET) main.o qfs.o
