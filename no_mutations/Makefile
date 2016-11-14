CC = g++ 
#-Wall
CFLAGS =-c -g -O3 -std=gnu++11 -Wno-unknown-pragmas

#extremely conserved element detection with single mutation
App=limes_no_mutations
all: $(App)
	
$(App):	main.o Serial.o IO.o Util.o Timing.o
	 		$(CC) main.o IO.o Util.o Serial.o Timing.o -o $(App)

build: 		main.o Serial.o Timing.o

install:	LIMES
			cp $(App) .
			
main.o:		main.cpp
			$(CC) $(CFLAGS) main.cpp

Serial.o:	Serial.cpp
			$(CC) $(CFLAGS) Serial.cpp

IO.o:		IO.cpp
			$(CC) $(CFLAGS) IO.cpp

Util.o:		Util.cpp
			$(CC) $(CFLAGS) Util.cpp

Timing.o:	Timing.cpp
			$(CC) $(CFLAGS) Timing.cpp

clean:
		rm	-rf	$(App)	*~	\#*\#	*.o		




