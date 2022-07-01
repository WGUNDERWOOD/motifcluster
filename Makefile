all: python R

.PHONY: python R

python:
	@cd python/ && make all

R:
	@cd R/ && make all
