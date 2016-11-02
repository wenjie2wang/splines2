objects := $(wildcard R/*.R) DESCRIPTION
man := $(wildcard man/*.Rd) NAMESPACE
dir := $(shell pwd)
version := $(shell grep "Version" DESCRIPTION | sed "s/Version: //")
pkg := $(shell grep "Package" DESCRIPTION | sed "s/Package: //")
tar := $(pkg)_$(version).tar.gz
checkLog := $(pkg).Rcheck/00check.log
rmd := vignettes/$(pkg)-intro.Rmd
vignettes := vignettes/$(pkg)-intro.html
cprt := COPYRIGHT


$(tar): $(objects)
	Rscript -e "library(methods); devtools::document();";
	R CMD build $(dir)

.PHONY: check
check: $(checkLog)

.PHONY: preview
preview: $(vignettes)

$(checkLog): $(tar)
	R CMD check --as-cran $(tar)

$(vignettes): $(rmd)
	Rscript -e "rmarkdown::render('$(rmd)')"

.PHONY: INSTALL
INSTALL: $(tar)
	R CMD INSTALL --build $(tar)

## update copyright year in HEADER, R script and date in DESCRIPTION
.PHONY: updateHeader
updateHeader: $(cprt)
	yr=$$(date +"%Y");\
	sed -i "s/Copyright (C) 2016-[0-9]\{4\}/Copyright (C) 2016-$$yr/" $(cprt);\
# add HEADER file if there is no header
	for Rfile in R/*.R; do \
	if ! grep -e 'Copyright (C)' $$Rfile ;\
	then cat $(cprt) $$Rfile > tmp ;\
	mv tmp $$Rfile;\
	fi;\
	yr=$$(date +"%Y");\
	sed -i "s/Copyright (C) 2016-[0-9]*/Copyright (C) 2016-$$yr/" $$Rfile;\
	done;\
	dt=$$(date +"%Y-%m-%d");\
	sed -i "s/Date: [0-9]\{4\}-[0-9]\{2\}-[0-9]\{2\}/Date: $$dt/" DESCRIPTION;

.PHONY: clean
clean:
	rm -rf *~ */*~ *.Rhistroy *.tar.gz *.Rcheck/ .\#*
