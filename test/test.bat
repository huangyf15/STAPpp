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
@if not %errorlevel%==0 (echo run test failed. & exit)
@stap++ ../data/truss > nul
@if not %errorlevel%==0 (echo run test failed. & exit)
@stap++ ../data/bar-6 > nul
@if not %errorlevel%==0 (echo run test failed. & exit)
@cd ../data/4Q
@python run-patch.py
@if not %errorlevel%==0 (echo test failed. & exit)
@cd ../3T
@python run-patch.py
@if not %errorlevel%==0 (echo test failed. & exit)
