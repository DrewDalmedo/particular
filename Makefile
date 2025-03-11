particular: clean
	mkdir -p build
	gcc  $@.c -o build/$@ -g -lSDL2 -lSDL2_image -lm

clean:
	rm -rf build

run: clean particular
	./build/particular
