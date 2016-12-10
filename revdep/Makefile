R := R CMD BATCH --vanilla

.PHONY: check
check:
	$(R) check.R &

.PHONY: clean
clean:
	rm -rf *.rds *.Rout *.Rhistory
