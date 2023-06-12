#-----------------
#   Compiler
#------------------

CC = g++

#------------------
#   Libraries
#------------------

ITINCLUDE3  =   #-I/... # when include files are requested
ITINCLUDE   = ${ITINCLUDE3}
ITLIB       = #-L/usr/lib/ # directory where the required libraries are

#----------------------------
#   Compilation options
#----------------------------  

INCLUDE_PATH = -I./ ${ITINCLUDE}
LDFLAGS      = #--enable-auto-import #-Wl#,--enable-auto-import 
LIBRARY      =  
WARN         = #-Wno-deprecated 
CPPFLAGS     = -ffast-math -funroll-loops ${WARN} ${INCLUDE_PATH}
DEFINES      = -Wextra -Wconversion

.cpp.o: ; ${CC} ${DEFINES}  -std=c++11  ${INCLUDE_PATH} ${WARN} -g -c $*.cpp

#----------------------------
#   Files to be compiled
#---------------------------- 
COMOBJ = input.o matrix.o hamiltonian.o rkr.o clist.o algorithm.o main.o 
OBJlangavEM = ${COMOBJ} sampling1order.o langavEM.o  
OBJlangavSEB = ${COMOBJ} sampling1order.o langavSEB.o
OBJlangavSEAC = ${COMOBJ} sampling1order.o langavSEAC.o
OBJlangavABAO = ${COMOBJ} sampling1order.o langavABAO.o
OBJlangavSOILEA = ${COMOBJ} sampling2order.o langavSOILEA.o
OBJlangavSOILEB = ${COMOBJ} sampling2order.o langavSOILEB.o

default:  ${OBJlangavEM}
	@echo 'Targets:' 
	@echo '  em     - Euler    R-KR (default)' 
	@echo '  seb    - SEB      R-KR' 
	@echo '  seac   - SEAC     R-KR' 
	@echo '  abao   - ABAO     R-KR' 
	@echo '  soilea - SOILE-A  R-KR' 
	@echo '  soileb - SOILE-B  R-KR' 
	@echo '  clean - remove object files' 
	@echo ' '
	@echo ' '
	${CC} ${CPPFLAGS} ${DEFINES} ${OBJlangavEM}  -o em ${LDFLAGS} ${LIBRARY} 
seb:  ${OBJlangavSEB}
	${CC} ${CPPFLAGS} ${DEFINES} ${OBJlangavSEB}  -o seb ${LDFLAGS} ${LIBRARY}
seac:  ${OBJlangavSEAC}
	${CC} ${CPPFLAGS} ${DEFINES} ${OBJlangavSEAC}  -o seac ${LDFLAGS} ${LIBRARY}
abao:  ${OBJlangavABAO}
	${CC} ${CPPFLAGS} ${DEFINES} ${OBJlangavABAO}  -o abao ${LDFLAGS} ${LIBRARY}
soilea:  ${OBJlangavSOILEA}
	${CC} ${CPPFLAGS} ${DEFINES} ${OBJlangavSOILEA}  -o soilea ${LDFLAGS} ${LIBRARY}
soileb:  ${OBJlangavSOILEB}
	${CC} ${CPPFLAGS} ${DEFINES} ${OBJlangavSOILEB}  -o soileb ${LDFLAGS} ${LIBRARY}


clean: 
	rm -f *.o em seb saec abao soilea soileb


