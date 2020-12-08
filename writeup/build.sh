#!/bin/bash
latexrun project.tex
latexrun supplement.tex
pdfunite project.pdf supplement.pdf final-project.pdf
