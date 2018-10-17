EXECNAME 	:= 	transient

CC				:=	gcc
CFLAGS 		:=	
LD_FLAGS 	:=	-lm
INCLUDES 	:= 	

src 			:= 	$(wildcard *.c)
obj 			:=	$(src:.c=.o)

all: $(EXECNAME)

# link
$(EXECNAME): $(obj)
		$(CC) -o $@ $^ $(LD_FLAGS)

# compile
%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $< 

.PHONY: clean
clean:
	rm $(obj) $(EXECNAME)
