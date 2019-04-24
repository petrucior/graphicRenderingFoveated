#
#	Makefile do projeto projectFoveaCuda
#

#
# A variável CC indica que compilador que será utilizado
#
CC= g++

#
# A variável CFLAGS indica que opção de compilação queremos
#
CFLAGS= -Wall -pedantic -ansi -g -c

#
# A variável LFLAGS indica que opção de compilação queremos
#
LFLAGS= -Wall -g -fopenmp

#
# A variável DIR indica o diretório dos arquivos de cabeçalho do projeto
#
DIR= ~/projectFoveaCuda

#
# A variável INCS indica o caminho dos arquivos de cabeçalho
#
INCS= -I../ -I/usr/local/include/opencv2 -I/usr/local/include

#
# A variável LIBS indica o caminho das bibliotecas e as bibliotecas usadas na ligação
#
LIBS=

#
# A variável OBJS indica os arquivos objetos comuns a todos os executáveis
#
OBJS= 

#
# A variável ALL indica os arquivos executáveis
#
all: mmf

#
# Regra de compilação e ligação do executáveis
# 
mmf:	mmf.cpp $(DIR)/mmf.h
	$(CC) $(LFLAGS) mmf.cpp `pkg-config --cflags --libs opencv` -o mmf $(LIBS) -O3

#
# Regra de compilação e ligação do executáveis com debug
# 
debug:	mmf.cpp $(DIR)/mmf.h
	$(CC) $(LFLAGS) mmf.cpp `pkg-config --cflags --libs opencv` -o mmf $(LIBS) -DDEBUG=1 -O3


#
# Regra que elimina arquivos intermediários
# 
clean:
	rm -f *.o *.cpp~

#
# Regra que elimina arquivos intermediários e executável
# 
realclean:
	rm -f *.o *.cpp~ mmf
