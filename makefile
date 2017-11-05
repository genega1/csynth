csynth: csynth.cpp
	g++ csynth.cpp

run:
	@$(MAKE) && ./$(csynth)
