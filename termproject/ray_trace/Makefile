
EXEC   = ray_trace.so

SRCS   = ray_trace.c 

OBJS   = $(SRCS:.c=.o)
INCL   = overhead.h

CFLAGS =  -shared -fPIC -O3 -g  #  -Wall
LNKCMD =  ld -L/usr/lib -L/usr/local/lib  -shared

LIBS   =  -lm 

CC     =  gcc 

$(EXEC): $(OBJS) 
	$(CC) -L/usr/lib -L/usr/local/lib  -shared $(OBJS) $(LIBS) $(CFLAGS)  -o $(EXEC)


$(OBJS): $(INCL) 



clean:
	rm $(OBJS) $(EXEC)










