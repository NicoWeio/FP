all: build/main.pdf

PLOTS = build/plt/spektrum.pdf

TABLES = build/tab/mu.tex

build/plt/spektrum.pdf: 1_spektrum.py | build/plt/
	python3 1_spektrum.py

build/tab/mu.tex: calc_mu.py
	python3 calc_mu.py

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
