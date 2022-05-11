all: build/main.pdf

PLOTS = \
	build/plt/tem_00.png \
	build/plt/tem_10.png \
	build/plt/2_stabilitaetsbedingung.pdf \
	build/plt/2_stabilitaetsbedingung_theorie.pdf \
	build/plt/3_tem_00.pdf \
	build/plt/3_tem_01.pdf \
	build/plt/4_polarisation.pdf \

TABLES = \
	build/tab/2_stabilitaetsbedingung_konkav_konkav.tex \
	build/tab/2_stabilitaetsbedingung_plan_konkav.tex \
	build/tab/3_tem_00.tex \
	build/tab/3_tem_01.tex \

build/plt/tem_% &: 0_tem_moden_theorie.py | build/plt
	MATPLOTLIBRC=../matplotlibrc TEXINPUTS=$$(pwd)/..: BUILD=true python 0_tem_moden_theorie.py

build/plt/2_stabilitaetsbedingung_theorie.pdf: 0_vorbereitung.py | build/plt
	MATPLOTLIBRC=../matplotlibrc TEXINPUTS=$$(pwd)/..: BUILD=true python 0_vorbereitung.py

build/plt/2_stabilitaetsbedingung.pdf \
build/tab/2_stabilitaetsbedingung_konkav_konkav.tex \
build/tab/2_stabilitaetsbedingung_plan_konkav.tex \
: 2_stabilitaetsbedingung.py $(wildcard dat/2_stabilitaetsbedingung_*.csv) | build/plt
	MATPLOTLIBRC=../matplotlibrc TEXINPUTS=$$(pwd)/..: BUILD=true python 2_stabilitaetsbedingung.py

build/plt/3_tem_%.pdf build/tab/3_tem_%.tex &: 3_tem_moden.py $(wildcard dat/3_tem_*.csv) | build/plt
	MATPLOTLIBRC=../matplotlibrc TEXINPUTS=$$(pwd)/..: BUILD=true python 3_tem_moden.py

build/plt/4_polarisation.pdf: 4_polarisation.py dat/4_polarisation.csv | build/plt
	MATPLOTLIBRC=../matplotlibrc TEXINPUTS=$$(pwd)/..: BUILD=true python 4_polarisation.py

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
