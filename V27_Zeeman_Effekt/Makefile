all: build/main.pdf

plots: \
		build/plt/1_magnet.pdf \


build/plt/1_magnet.pdf: 1_magnet.py | build/plt
		MATPLOTLIBRC=../matplotlibrc TEXINPUTS=$$(pwd)/..: python3 1_magnet.py


build/main.pdf: main.tex lit.bib plots | build
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
