all:  docs install

docs:
	R -e "devtools::document()"
build:
	(cd ..; R CMD build --no-build-vignettes TrenaProjectIGAP)

install:
	(cd ..; R CMD INSTALL TrenaProjectIGAP)

check:
	(cd ..; R CMD check `ls -t TrenaProjectIGAP) | head -1`)

tests:
	for x in inst/unitTests/test_*.R; do echo $$x; R -f $$x; done
