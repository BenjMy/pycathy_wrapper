# Build, package, test, and clean

PROJECT=pyCATHY
STYLE_CHECK_FILES=lib


format: isort black

check: isort-check black-check license-check flake8

black:
	black $(STYLE_CHECK_FILES)

black-check:
	black --check $(STYLE_CHECK_FILES)

isort:
	isort $(STYLE_CHECK_FILES)

isort-check:
	isort --check $(STYLE_CHECK_FILES)
	
	
