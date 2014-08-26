all: main.cpp
	g++ -std=c++11 main.cpp read_config.cpp read_data.cpp edgelist.cpp output.cpp dynamics.cpp -fopenmp -lyaml-cpp -pthread -o main	

clean:
	rm main
