all: build/main.pdf

PLOTS = \
	build/plt/1_koinzidenz.pdf \
	build/plt/2_mca.pdf \
	build/plt/3_lebensdauer_linear.pdf \
	build/plt/3_lebensdauer_log.pdf \

TABLES = \
	build/tab/1_koinzidenz.tex \
	build/tab/2_mca.tex \
	# build/tab/3_lebensdauer.tex \

PYTHONENV = MATPLOTLIBRC=../matplotlibrc TEXINPUTS=$$(pwd)/..: BUILD=true

build/plt/1_koinzidenz.pdf build/tab/1_koinzidenz.tex &: 1_koinzidenz.py 1_koinzidenz.csv | build/plt build/tab
	$(PYTHONENV) python 1_koinzidenz.py

build/plt/2_mca.pdf build/tab/2_mca.tex &: 2_mca.py 2_mca.csv | build/plt build/tab
	$(PYTHONENV) python 2_mca.py

build/plt/3_lebensdauer_linear.pdf build/plt/3_lebensdauer_log.pdf build/tab/3_lebensdauer.tex &: 3_lebensdauer.py 3_lebensdauer.csv | build/plt build/tab
	$(PYTHONENV) python 3_lebensdauer.py

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
	mkdir -p build
	mkdir -p build/plt
	mkdir -p build/tab

clean:
	rm -rf build

.PHONY: all clean plots tables

#NOTE: Vollständiges Rekompilieren lässt sich mit `make -B` erzwingen
