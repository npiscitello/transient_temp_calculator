EXECNAME 	:= 	transient
THREADED 	:=	$(EXECNAME)_threaded

CC				:=	gcc
CFLAGS 		:=	
LD_FLAGS 	:=	-lpthread
INCLUDES 	:= 	

src 			:= 	$(wildcard *.c)
obj 			:=	$(src:.c=.o)

all: $(EXECNAME)
threaded: $(THREADED)

# link
#$(EXECNAME): $(obj)
$(EXECNAME): $(EXECNAME).o
		$(CC) -o $@ $^ $(LD_FLAGS)

$(THREADED): $(THREADED).o
		$(CC) -o $@ $^ $(LD_FLAGS)

# compile
%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $< 

.PHONY: clean
clean:
	rm $(obj) $(EXECNAME) $(THREADED)
