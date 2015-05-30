
# Need to install profiling.
# pip install git+https://github.com/what-studio/profiling.git
profiling:
	python -m profiling profile pymoog.py --dump=pymoog.prf
	python -m profiling view pymoog.prf
