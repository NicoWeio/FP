all: build/main.pdf

PLOTS = \
	build/plt/1_detektor.pdf \
	build/plt/2_z1.pdf \
	build/plt/4_rocking1_linear.pdf \
	build/plt/4_rocking1_log.pdf \
	build/plt/schichtdicke_a.pdf \
	build/plt/schichtdicke_b.pdf \

TABLES = \
	# build/tab/1_detektor.tex \
	# build/tab/2_z1.tex \
	# build/tab/4_rocking1.tex \

PYTHONENV = MATPLOTLIBRC=../matplotlibrc TEXINPUTS=$$(pwd)/..: BUILD=true PLOTS=true

$(PLOTS) &: main.py $(wildcard code/*.py) data_Mampfzwerg/1_detektor.txt | build/plt build/tab
	$(PYTHONENV) python main.py

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
