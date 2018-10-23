all:  docs install

docs:
	R -e "devtools::document()"
build:
	(cd ..; R CMD build --no-build-vignettes trenaProjectIGAP)

install:
	(cd ..; R CMD INSTALL trenaProjectIGAP)

check:
	(cd ..; R CMD check `ls -t trenaProjectIGAP) | head -1`)

