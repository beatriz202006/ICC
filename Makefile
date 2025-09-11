# PROGRAMA
    PROG = resolveEDO
    OBJS = $(PROG).o edo.o utils.o gauss_seidel.o # mod1.o mod2.o
     VERIF = verificaEP02

# Compilador
    CC     = gcc

# Acrescentar onde apropriado as opções para incluir uso da biblioteca LIKWID
    CFLAGS = -O0
    LFLAGS = -lm

# Lista de arquivos para distribuição. Acrescentar mais arquivos se necessário.
DISTFILES = *.c *.h LEIAME* Makefile *.dat
DISTDIR = ${USER}

.PHONY: clean purge dist all

%.o: %.c %.h utils.h
	$(CC) -c $(CFLAGS) -o $@ $<

$(PROG): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

$(VERIF): $(VERIF).c
	$(CC) -Wno-format -o $@ $<

testeFormato: $(PROG) $(VERIF)
	@cat teste.dat | ./$(PROG) | ./$(VERIF)

clean:
	@echo "Limpando sujeira ..."
	@rm -f *~ *.bak

purge:  clean
	@echo "Limpando tudo ..."
	@rm -f core a.out $(OBJS)
	@rm -f $(PROG) $(VERIF) $(DISTDIR) $(DISTDIR).tgz

dist: purge
	@echo "Gerando arquivo de distribuição ($(DISTDIR).tgz) ..."
	@ln -s . $(DISTDIR)
	@tar -chzvf $(DISTDIR).tgz $(addprefix ./$(DISTDIR)/, $(DISTFILES))
	@rm -f $(DISTDIR)
