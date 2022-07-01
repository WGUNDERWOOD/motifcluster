# Makefile: python

all: test lint doc tutorial dist

.PHONY: test lint doc tutorial dist uninstall install upload_test_pypi upload_pypi clean

test:
	@echo -e "\e[0;35m\033[1mTesting python package...\e[0;30m\033[0m"
	@pytest --cov=motifcluster --profile-svg tests/

lint:
	@echo -e "\e[0;35m\033[1mLinting python package...\e[0;30m\033[0m"
	@pylint -j 8 --rcfile=.pylintrc motifcluster

doc:
	@echo -e "\e[0;35m\033[1mMaking python documentation...\e[0;30m\033[0m"
	@cd doc/ && make html

tutorial:
	@echo -e "\e[0;35m\033[1mMaking python tutorial...\e[0;30m\033[0m"
	@cd tutorial && latexmk -pdf -quiet -rc-report- motifcluster_tutorial.tex
	@cd tutorial && pythontex motifcluster_tutorial.tex
	@cd tutorial && latexmk -pdf -quiet -rc-report- motifcluster_tutorial.tex
	@optpdf tutorial/motifcluster_tutorial.pdf
	@optpdf tutorial/motifcluster_tutorial.pdf

dist:
	@echo -e "\e[0;35m\033[1mMaking python distribution...\e[0;30m\033[0m"
	@python setup.py sdist bdist_wheel
	@python -m twine check dist/*

uninstall:
	@echo -e "\e[0;35m\033[1mUninstalling python package...\e[0;30m\033[0m"
	@pip uninstall motifcluster

install: clean dist
	@echo -e "\e[0;35m\033[1mInstalling python package...\e[0;30m\033[0m"
	@pip install --user dist/motifcluster-*.tar.gz

upload_test_pypi:
	@echo -e "\e[0;35m\033[1mUploading python package to test PyPI...\e[0;30m\033[0m"
	@python -m twine upload --repository testpypi -u __token__ dist/*

upload_pypi:
	@echo -e "\e[0;35m\033[1mUploading python package to PyPI...\e[0;30m\033[0m"
	@python -m twine upload -u __token__ dist/*

clean:
	@echo -e "\e[0;35m\033[1mCleaning up python package...\e[0;30m\033[0m"
	@rm -f .coverage
	@rm -rf dist/
	@rm -rf build
	@rm -rf motifcluster.egg-info/
	@rm -rf doc/_build/
	@rm -rf motifcluster/__pycache__/
	@rm -rf tests/__pycache__/
	@rm -rf prof/
	@rm -rf .pytest_cache/
	@cd tutorial && texclean && rm -rf pythontex-files-motifcluster_tutorial/
