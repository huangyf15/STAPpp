@cd ..
@if not exist build goto makenew
:resume
@cd build
@goto test

:makenew
@mkdir build
@cd build
@cmake -G "MinGW Makefiles" ../src
@goto resume

:test
@make rebuild_cache
@make
@stap++ ../data/test_truss_22.dat > nul
@if not %errorlevel%==0 goto fail
@stap++ ../data/truss > nul
@if not %errorlevel%==0 goto fail
@stap++ ../data/bar-6 > nul
@if not %errorlevel%==0 goto fail
@cd ../data/4Q
@python run-patch.py
@if not %errorlevel%==0 goto fail
@cd ../3T
@python run-patch.py
@if not %errorlevel%==0 goto fail
@cd ../Beam
@python run-patch.py
@if not %errorlevel%==0 goto fail
@cd ../8H
@python run-patch.py
@if not %errorlevel%==0 goto fail
@cd ../shell
@python run-patch.py
@if not %errorlevel%==0 goto fail
@cd ../plate
@python run-patch.py
@if not %errorlevel%==0 goto fail
@echo test success
@goto quit

:fail
@echo test failed
:quit