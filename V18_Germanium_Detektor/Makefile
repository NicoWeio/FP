all: build/main.pdf

PLOTS = build/plt/effizienz.pdf

# TABLES = build/tab/mu.tex \
# 	build/tab/wuerfel1.tex \
# 	build/tab/wuerfel2.tex \
# 	build/tab/wuerfel3.tex \
# 	build/tab/wuerfel4.tex \

PYTHONENV = MATPLOTLIBRC=../matplotlibrc TEXINPUTS=$$(pwd)/..: BUILD=true PLOTS=true

build/plt/effizienz.pdf &: auswertung.py | build/plt/
	$(PYTHONENV) python3 auswertung.py

# build/tab/wuerfel%.tex &: auswertung2.py | build/tab/
# 	python3 auswertung2.py

# build/tab/mu.tex: calc_mu.py
# 	python3 calc_mu.py

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
