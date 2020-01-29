CC = gcc

segmentation: segmentation.c
	$(CC) -O3 -o segmentation segmentation.c -lm

run-mcm:
	./segmentation mcm < ./params-mcm.txt
	animate -delay 20 *_head-out.pgm

run-cv:
	./segmentation cv < ./params-cv.txt