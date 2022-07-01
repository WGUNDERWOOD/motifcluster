# Makefile: python

all: test lint doc tutorial dist

.PHONY: test lint doc tutorial dist uninstall install upload_test_pypi upload_pypi clean

test:
	@pytest --cov=motifcluster --profile-svg tests/

lint:
	@pylint -j 8 --rcfile=.pylintrc motifcluster

doc:
	@cd doc/ && make html

tutorial:
	@cd tutorial && latexmk -pdf -quiet -rc-report- motifcluster_tutorial.tex
	@cd tutorial && pythontex motifcluster_tutorial.tex
	@cd tutorial && latexmk -pdf -quiet -rc-report- motifcluster_tutorial.tex

dist:
	@python setup.py sdist bdist_wheel
	@python -m twine check dist/*

uninstall:
	@pip uninstall motifcluster

install: clean dist
	@pip install --user dist/motifcluster-*.tar.gz

upload_test_pypi:
	@python -m twine upload --repository testpypi -u __token__ dist/*

upload_pypi:
	@python -m twine upload -u __token__ dist/*

clean:
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
