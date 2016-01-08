SHELL = /bin/bash

install:
	@./check_moog.sh
	@echo "Downloading the model atmosphere..."
	@mkdir -p models
	@mkdir -p linelist
	@wget www.astro.up.pt/~dandreasen/kurucz95.tar.gz
	@tar zxf kurucz95.tar.gz
	@rm -rf models/kurucz95*
	@mv -f kurucz95.tar.gz models
	@mv -f kurucz95 models
	@echo "Atmosphere models installed in dir: models"
	@echo "Installing dependencies..."
	@pip install -r requirements.txt
	@echo "Dependencies installed"
	@echo "\nMOOGme is ready to go.\nType python MOOGme.py to start"
