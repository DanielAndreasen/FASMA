SHELL = /bin/bash

install:
	@./check_moog.sh
	@echo "Installing the model atmospheres..."
	@mkdir -p models
	@mkdir -p results
	@mkdir -p spectra
	@tar zxf models/kurucz95.tar.gz
	@mv -f kurucz95 models
	@tar zxf models/marcs.tar.gz
	@mv -f marcs models
	@tar zxf models/APOGEE_kurucz.tar.gz
	@mv -f APOGEE_kurucz models/apogee_kurucz
	@echo "Atmosphere models installed in dir: models"
	@echo "Installing dependencies..."
	@pip install -r requirements.txt
#	@conda install -c anaconda wxpython=3.0.0.0
	@echo "Installing ARES"
	@cd ARES; make install; cd ..
	@echo "Installing TMCALC"
	@cd TMCALC; make; cd ..
	@echo "Dependencies installed"
	@echo ""
	@echo "MOOGme is successfully installed!"
	@echo "Type"
	@echo "    python MOOGme.py"
	@echo "to start. Happy spectroscopying :)"
