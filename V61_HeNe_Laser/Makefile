all: build/main.pdf

plots: \
	build/plt/tem_00.png \
	build/plt/tem_10.png \
	build/plt/2_stabilitaetsbedingung.pdf \
	build/plt/2_stabilitaetsbedingung_theorie.pdf \
	build/plt/3_tem_00.pdf \
	build/plt/3_tem_01.pdf \
	build/plt/4_polarisation.pdf \

tables: \
	build/tab/3_tem_00.tex \
	build/tab/3_tem_01.tex \

build/plt/tem_% &: 0_tem_moden_theorie.py | build/plt
	MATPLOTLIBRC=../matplotlibrc TEXINPUTS=$$(pwd)/..: BUILD=true python 0_vorbereitung.py

build/plt/2_stabilitaetsbedingung_theorie.pdf: 0_vorbereitung.py | build/plt
	MATPLOTLIBRC=../matplotlibrc TEXINPUTS=$$(pwd)/..: BUILD=true python 0_vorbereitung.py

build/plt/2_stabilitaetsbedingung.pdf: 2_stabilitaetsbedingung.py $(wildcard dat/2_stabilitaetsbedingung_*.csv) | build/plt
	MATPLOTLIBRC=../matplotlibrc TEXINPUTS=$$(pwd)/..: BUILD=true python 2_stabilitaetsbedingung.py

build/plt/3_tem_%.pdf build/tab/3_tem_%.tex &: 3_tem_moden.py $(wildcard dat/3_tem_*.csv) | build/plt
	MATPLOTLIBRC=../matplotlibrc TEXINPUTS=$$(pwd)/..: BUILD=true python 3_tem_moden.py

build/plt/4_polarisation.pdf: 4_polarisation.py dat/4_polarisation.csv | build/plt
	MATPLOTLIBRC=../matplotlibrc TEXINPUTS=$$(pwd)/..: BUILD=true python 4_polarisation.py

build/main.pdf: main.tex lit.bib plots tables | build
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
