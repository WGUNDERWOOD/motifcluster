# Makefile: Julia

all: test format doc tutorial

.PHONY: test format doc tutorial clean

test:
	@echo -e "\e[0;35m\033[1mTesting Julia package...\e[0;30m\033[0m"
	@cd MotifCluster && julia --depwarn=yes --project -e 'using Pkg; Pkg.test(coverage=true)'

format:
	@echo -e "\e[0;35m\033[1mFormatting Julia package...\e[0;30m\033[0m"
	@cd MotifCluster && julia -e 'using JuliaFormatter; format(".")'

doc:
	@echo -e "\e[0;35m\033[1mMaking Julia documentation...\e[0;30m\033[0m"
	@cd MotifCluster/docs && julia make.jl

tutorial:
	@echo -e "\e[0;35m\033[1mMaking Julia tutorial...\e[0;30m\033[0m"
	@cd MotifCluster/tutorial && make

clean:
	@echo -e "\e[0;35m\033[1mCleaning up Julia package...\e[0;30m\033[0m"
	@cd MotifCluster/tutorial && texclean
	@rm -r MotifCluster/docs/build
