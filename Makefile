project:
	nvcc -o out/project -std=c++11 src/main.cu

clean:
	rm -rf out/*.txt

remove_project:
	rm -rf out/project
