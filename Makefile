objects := $(wildcard R/*.R) DESCRIPTION
version := $(shell grep "Version" DESCRIPTION | awk '{print $$NF}')
pkg := $(shell grep "Package" DESCRIPTION | awk '{print $$NF}')
tar := $(pkg)_$(version).tar.gz
tinytest := $(wildcard tests/testthat/*.R)
checkLog := $(pkg).Rcheck/00check.log
rmd := $(wildcard vignettes/*.Rmd)
vignettes := $(patsubst %.Rmd,%.html,$(rmd))


.PHONY: check
check: $(checkLog)

.PHONY: build
build: $(tar)

.PHONY: install
install: $(tar)
	R CMD INSTALL $(tar)

.PHONY: preview
preview: $(vignettes)

.PHONY: pkgdown
pkgdown:
	Rscript -e "library(methods); pkgdown::build_site();"

$(tar): $(objects)
	@$(RM) -rf src/RcppExports.cpp R/RcppExports.R
	@Rscript -e "library(methods);" \
	-e "Rcpp::compileAttributes()" \
	-e "devtools::document();";
	@$(MAKE) updateTimestamp
	R CMD build .

$(checkLog): $(tar) $(tinytest)
	R CMD check --as-cran $(tar)

.PHONY: check-as-cran
check-as-cran: $(tar)
	R CMD check --as-cran $(tar)

vignettes/%.html: vignettes/%.Rmd
	Rscript -e "library(methods); rmarkdown::render('$?')"

## update copyright year in HEADER, R script and date in DESCRIPTION
.PHONY: updateTimestamp
updateTimestamp:
	@bash misc/update_timestamp.sh

## make tags
.PHONY: TAGS
TAGS:
	Rscript -e "utils::rtags(path = 'R', ofile = 'TAGS')"
	gtags

.PHONY: clean
clean:
	@$(RM) -rf *~ */*~ *.Rhistroy *.tar.gz src/*.so src/*.o *.Rcheck/ .\#*
