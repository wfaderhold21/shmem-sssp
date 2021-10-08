CC=shmemcc
CFLAGS=-O2 -g -Wall -Werror -std=c99 -DBF_DEBUG=0
LDFLAGS=

all: shmem_bf io.o

io.o:	io.c

shmem_bf:	shmem-bf.c io.o common.h
	${CC} ${CFLAGS} shmem-bf.c io.o -o shmem-bf ${LDFLAGS}
	
clean:
	rm shmem_bf io.o
