particular: clean
	gcc  $@.c -o $@ -g -lSDL2 -lSDL2_image -lm

clean:
	rm -f particular

run: clean particular
	./particular
