all: build/main.pdf

PLOTS = \
	build/plt/g_F.pdf \

TABLES = \
	build/tab/messwerte.tex \

PYTHONENV = MATPLOTLIBRC=../matplotlibrc TEXINPUTS=$$(pwd)/..: BUILD=true

build/plt/g_F.pdf build/tab/messwerte.tex &: auswertung.py | build/plt build/tab
	$(PYTHONENV) python auswertung.py

build/main.pdf: ../header_common.tex main.tex $(wildcard content/*.tex) lit.bib $(PLOTS) $(TABLES) | build
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
	mkdir -p $@

clean:
	rm -rf build

.PHONY: all clean plots tables

#NOTE: Vollständiges Rekompilieren lässt sich mit `make -B` erzwingen
