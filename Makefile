# Beatriz Pontes Camargo 
# GRR 20242966

# PROGRAMA
    PROG = resolveEDO
    OBJS = $(PROG).o edo.o utils.o gauss_seidel.o # mod1.o mod2.o
    VERIF = verificaEP02

# Compilador
    CC     = gcc

# Para incluir uso da biblioteca LIKWID
CFLAGS = -O0 -DLIKWID_PERFMON -I${LIKWID_INCLUDE}
LFLAGS = -L${LIKWID_LIB} -llikwid -lm


# Lista de arquivos para distribuição
DISTFILES = *.c *.h LEIAME* Makefile *.dat likwid_script.sh
DISTDIR = ${USER}

.PHONY: clean purge dist all perm

%.o: %.c %.h utils.h
	$(CC) -c $(CFLAGS) -o $@ $<

$(PROG): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

$(VERIF): $(VERIF).c
	$(CC) -Wno-format -o $@ $<

testeFormato: $(PROG) $(VERIF)
	@cat teste.dat | ./$(PROG) | ./$(VERIF)

all: $(PROG)
	chmod +x likwid.sh

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

perm:
	chmod +x likwid_script.sh
