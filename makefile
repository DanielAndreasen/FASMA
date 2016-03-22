SHELL = /bin/bash

install:
	@./check_moog.sh
	@echo "Downloading the model atmosphere..."
	@mkdir -p models
	@mkdir -p linelist
	@mkdir -p spectra
	@wget www.astro.up.pt/~dandreasen/kurucz95.tar.gz
	@tar zxf kurucz95.tar.gz
	@rm -rf models/kurucz95*
	@mv -f kurucz95.tar.gz models
	@mv -f kurucz95 models
	@echo "Atmosphere models installed in dir: models"
	@echo "Installing dependencies..."
	@pip install -r requirements.txt
#	@pip install -U --trusted-host wxpython.org --pre -f http://wxpython.org/Phoenix/snapshot-builds/ wxPython_Phoenix
	@echo "Dependencies installed"
	@echo ""
	@echo "MOOGme is successfully installed!"
	@echo "Type"
	@echo "    python MOOGme.py"
	@echo "to start. Happy spectroscopying :)"
