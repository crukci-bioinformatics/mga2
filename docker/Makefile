
version ?= latest

build:
	docker build --tag crukcibioinformatics/mga2:${version} .

release: build
	docker push crukcibioinformatics/mga2:${version}

singularity: build
	singularity build mga2-${version}.sif docker-daemon://crukcibioinformatics/mga2:${version}

