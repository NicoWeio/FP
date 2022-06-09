all: build/main.pdf

PLOTS = build/plt/waermekapazitaet.pdf
TABLES = build/tab/mess.tex

build/plt/waermekapazitaet.pdf build/tab/mess.tex: auswertung.py dat/messung.txt | build/plt build/tab
	BUILD=true MATPLOTLIBRC=../matplotlibrc TEXINPUTS=$$(pwd)/..: python3 auswertung.py

build/main.pdf: main.tex lit.bib $(wildcard content/*.tex) $(PLOTS) $(TABLES) | build
	  TEXINPUTS=build: \
	  BIBINPUTS=build: \
	  max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=nonstopmode \
	  --halt-on-error \
	main.tex

build build/plt build/tab &:
	mkdir -p build
	mkdir -p build/plt
	mkdir -p build/tab

clean:
	rm -rf build

.PHONY: all clean plots tables

#NOTE: Vollständiges Rekompilieren lässt sich mit `make -B` erzwingen
