CC=g++
CFLAGS=-I.
#DEPS = hellomake.h
OBJS = main.o
TARGET = qfs
LINK = -lcrypto -lssl

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(TARGET): $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(LINK)

clean:
	rm -f $(TARGET) $(OBJS)
