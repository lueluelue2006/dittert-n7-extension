PYTHON ?= python3
PDFLATEX ?= pdflatex

.PHONY: test verify verify-quick pdf hashes all clean

test:
	$(PYTHON) -m unittest discover -s tests -v

verify:
	$(PYTHON) run_all_certificates.py

verify-quick:
	$(PYTHON) run_all_certificates.py --skip-independent --skip-optimized

pdf:
	$(PDFLATEX) -interaction=nonstopmode -halt-on-error main.tex
	$(PDFLATEX) -interaction=nonstopmode -halt-on-error main.tex

hashes:
	$(PYTHON) run_all_certificates.py --hashes-only

all: test pdf verify

clean:
	rm -f main.aux main.log main.out main.toc
	rm -rf __pycache__ tests/__pycache__
