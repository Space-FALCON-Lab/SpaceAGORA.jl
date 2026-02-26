CC = g++

./lib/libgtest.a: gtest-all.o
	mkdir ./lib
	ar -rv ./lib/libgtest.a gtest-all.o
	rm gtest-all.o

gtest-all.o: ./src/gtest-all.cc
	${CC} -isystem ./include -I. -pthread -c ./src/gtest-all.cc

clean:
	rm -Rf ./*.o ./lib
