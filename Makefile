objects := $(wildcard R/*.R) DESCRIPTION
# man := $(wildcard man/*.Rd) NAMESPACE
dir := $(shell pwd)
version := $(shell grep "Version" DESCRIPTION | sed "s/Version: //")
pkg := $(shell grep "Package" DESCRIPTION | sed "s/Package: //")
tar := $(pkg)_$(version).tar.gz
checkLog := $(pkg).Rcheck/00check.log

rmd := vignettes/$(pkg)-intro.Rmd
vignettes := vignettes/$(pkg)-intro.html
cprt := COPYRIGHT


.PHONY: check
check: $(checkLog)

.PHONY: build
all: $(tar)

.PHONY: preview
preview: $(vignettes)

$(tar): $(objects)
	Rscript -e "library(methods); devtools::document();";
	R CMD build $(dir)

$(checkLog): $(tar)
	R CMD check --as-cran $(tar)

$(vignettes): $(rmd)
	Rscript -e "rmarkdown::render('$(rmd)')"

.PHONY: install
install: $(tar)
	R CMD INSTALL $(tar)

## update copyright year in HEADER, R script and date in DESCRIPTION
.PHONY: updateDate
updateDate: $(cprt)
	yr=$$(date +"%Y");\
	sed -i "s/Copyright (C) 2016[-]*[0-9]*/Copyright (C) 2016-$$yr/" $(cprt);\
# add HEADER file if there is no header
	for Rfile in R/*.R; do \
	if ! grep -e 'Copyright (C)' $$Rfile ;\
	then cat $(cprt) $$Rfile > tmp ;\
	mv tmp $$Rfile;\
	fi;\
	yr=$$(date +"%Y");\
	sed -i "s/Copyright (C) 2016[-]*[0-9]*/Copyright (C) 2016-$$yr/" $$Rfile;\
	done;\
	dt=$$(date +"%Y-%m-%d");\
	sed -i "s/Date: [0-9]\{4\}-[0-9]\{1,2\}-[0-9]\{1,2\}/Date: $$dt/" DESCRIPTION;

.PHONY: clean
clean:
	rm -rf *~ */*~ *.Rhistroy *.tar.gz *.Rcheck/ .\#*
