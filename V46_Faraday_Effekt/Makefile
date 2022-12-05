all: build/main.pdf

PYTHONENV = MATPLOTLIBRC=../matplotlibrc TEXINPUTS=$$(pwd)/..: BUILD=true

build/plt/magnet.pdf &: 1_magnet.py data/magnet.csv | build/plt
	$(PYTHONENV) python 1_magnet.py

build/plt/rotation.pdf &: 2_rotation.py data/probe_1.csv | build/plt
	$(PYTHONENV) python 2_rotation.py

#build/plt/3_lebensdauer_linear.pdf build/plt/3_lebensdauer_log.pdf build/tab/3_lebensdauer.tex &: 3_lebensdauer.py 3_lebensdauer.csv | build/plt build/tab
#	$(PYTHONENV) python 3_lebensdauer.py

build/main.pdf: main.tex lit.bib build/plt/magnet.pdf build/plt/rotation.pdf | build
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
