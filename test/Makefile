#
#	Makefile do projeto Ray Tracing
#

#
# A variável CC indica que compilador será utilizado
#
CC= g++

#
# A variável CFLAGS indica que opções de compilação queremos
#
CFLAGS=	-Wall -pedantic -ansi -g -c

#
# A variável LFLAGS indica que opções de compilação queremos
# MACOS: -openmp  Unix: -fopenmp 
LFLAGS=	-Wall -g #-openmp #-fopenmp 

#
# A variável DIR indica o diretório dos arquivos de cabeçalho do projeto
#
DIR= ~/raytracing

#
# A variável INCS indica o caminho dos arquivos de cabeçalho
#
INCS= -I../ -I/usr/local/include/opencv2 -I/usr/local/include

#
# A variável LIBS indica o caminho das bibliotecas e as bibliotecas usadas na ligação
#
#LIBS= -lm -lGL -lglut -lGLU
LIBS= /System/Library/Frameworks/GLUT.framework/GLUT /System/Library/Frameworks/OpenGL.framework/OpenGL
#LIBS= -lopengl32 -lglut32 -lglu32

#
# A variável OBJS indica os arquivos objetos
#
OBJS= vetor.o raio.o luz.o objeto.o cena.o textura.o ray_tracing.o #main.o

#
# A variável ALL indica os arquivos executáveis
#
all: main

#
# Regra de compilação e ligação do executável
# 
main: $(OBJS) ../renderMMF.h main.cpp
	$(CC) $(LFLAGS) $(OBJS) main.cpp `pkg-config --cflags --libs opencv` -o main $(LIBS) -O3

#
# Regra de compilação e ligação do executáveis com debug
#
debug:	$(OBJS) ../renderMMF.h main.cpp
	$(CC) $(LFLAGS) $(OBJS) main.cpp `pkg-config --cflags --libs opencv` -o main $(LIBS) -DDEBUG=1 -O3

#
# Regra de compilação do arquivo objeto vetor.o
# 
vetor.o: $(DIR)/vetor.cpp $(DIR)/vetor.hpp
	$(CC) $(CFLAGS) $(DIR)/vetor.cpp -o vetor.o

#
# Regra de compilação do arquivo objeto raio.o
# 
raio.o: $(DIR)/raio.cpp $(DIR)/raio.hpp
	$(CC) $(CFLAGS) $(DIR)/raio.cpp -o raio.o
#
# Regra de compilação do arquivo objeto luz.o
# 
luz.o: $(DIR)/luz.cpp $(DIR)/luz.hpp
	$(CC) $(CFLAGS) $(DIR)/luz.cpp -o luz.o

#
# Regra de compilação do arquivo objeto raio.o
# 
objeto.o: $(DIR)/objeto.cpp $(DIR)/objeto.hpp $(DIR)/textura.cpp $(DIR)/textura.hpp
	$(CC) $(CFLAGS) $(DIR)/objeto.cpp -o objeto.o

#
# Regra de compilação do arquivo objeto cena.o
# 
cena.o: $(DIR)/cena.cpp $(DIR)/cena.hpp
	$(CC) $(CFLAGS) $(DIR)/cena.cpp -o cena.o

#
# Regra de compilação do arquivo objeto cena.o
# 
textura.o: $(DIR)/textura.cpp $(DIR)/textura.hpp
	$(CC) $(CFLAGS) $(DIR)/textura.cpp -o textura.o

#
# Regra de compilação do arquivo objeto raytracing.o
# 
ray_tracing.o: $(DIR)/ray_tracing.cpp $(DIR)/ray_tracing.hpp $(DIR)/vetor.cpp $(DIR)/vetor.hpp $(DIR)/raio.cpp $(DIR)/raio.hpp $(DIR)/luz.cpp $(DIR)/luz.hpp $(DIR)/objeto.cpp $(DIR)/objeto.hpp $(DIR)/cena.cpp $(DIR)/cena.hpp $(DIR)/textura.cpp $(DIR)/textura.hpp
	$(CC) $(CFLAGS) $(DIR)/ray_tracing.cpp -o ray_tracing.o $(LIBS)

#
# Regra de compilação do arquivo objeto main.o
# 
main.o: main.cpp $(DIR)/vetor.cpp $(DIR)/vetor.hpp $(DIR)/raio.cpp $(DIR)/raio.hpp $(DIR)/luz.cpp $(DIR)/luz.hpp $(DIR)/objeto.cpp $(DIR)/objeto.hpp $(DIR)/cena.cpp $(DIR)/cena.hpp $(DIR)/textura.cpp $(DIR)/textura.hpp $(DIR)/ray_tracing.cpp $(DIR)/ray_tracing.hpp ../renderMMF.h
	$(CC) $(CFLAGS) main.cpp `pkg-config --cflags --libs opencv` -o main.o

#
# Regra que elimina arquivos intermediários
# 
clean:
	rm -f *.o *.~

#
# Regra que elimina arquivos intermediários e executável
# 
realclean:
	rm -f *.o *.~ main

#
# Regra para executar - Abre processo mais nao exclui
# 
exec:
	./main
