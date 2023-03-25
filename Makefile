objects := $(wildcard R/*.R) DESCRIPTION
version := $(shell grep -E "^Version:" DESCRIPTION | awk '{print $$NF}')
pkg := $(shell  grep -E "^Package:" DESCRIPTION | awk '{print $$NF}')
tar := $(pkg)_$(version).tar.gz
tinytest := $(wildcard inst/tinytest/*.R)
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

.PHONY: check-rcpp
check-rcpp: $(tar)
	R CMD INSTALL $(tar)
	Rscript inst/run_rcpp_test.R > check-rcpp.Rout &

.PHONY: check-revdep
check-revdep: $(tar)
	@mkdir -p revdep
	@rm -rf revdep/{*.Rcheck,*.tar.gz}
	@cp $(tar) revdep
	R CMD BATCH --no-save --no-restore misc/revdep_check.R &

vignettes/%.html: vignettes/%.Rmd
	Rscript -e "library(methods); rmarkdown::render('$?')"

.PHONY: readme
readme: README.md
README.md: README.Rmd
	@Rscript -e "rmarkdown::render('$<')"

## update copyright year in HEADER, R script and date in DESCRIPTION
.PHONY: updateTimestamp
updateTimestamp:
	@bash misc/update_timestamp.sh

## make tags
.PHONY: tags
tags:
	Rscript -e "utils::rtags(path = 'R', ofile = 'TAGS')"
	gtags

.PHONY: clean
clean:
	@$(RM) -r *~ */*~ *.Rhistroy *.tar.gz src/*.so src/*.o \
	*.Rcheck/ *.Rout .\#* *_cache
