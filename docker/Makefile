
version ?= latest

build:
	docker build --tag crukcibioinformatics/mga2:${version} .

release: build
	docker push crukcibioinformatics/mga2:${version}

