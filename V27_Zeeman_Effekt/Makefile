all: build/main.pdf

plots: \
		build/plt/1_magnet.pdf \
		build/plt/rot_1.pdf \
		build/plt/rot_2.pdf \


plots &: 1_auswertung.py 1a_bildanalyse.py | build/plt
		BUILD=true MATPLOTLIBRC=../matplotlibrc TEXINPUTS=$$(pwd)/..: python3 1_auswertung.py


build/main.pdf: main.tex lit.bib $(wildcard content/*.tex) $(wildcard build/tab/*.tex) plots | build
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
